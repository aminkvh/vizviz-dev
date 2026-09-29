//! 3D-SNFG glycan detection and per-residue placement geometry.
//!
//! Ports the residue-recognition table, shape/color assignment, and
//! ring-orientation math of the `3D-SNFG.tcl` script (version 1), written by
//! David F. Thieker and Jodi A. Hadden (<http://glycam.org/3d-snfg>).
//! Cite: Thieker, Hadden, Schulten & Woods (2016) "3D Implementation of
//! the Symbol Nomenclature for Graphical Representation of Glycans",
//! Glycobiology 26(8):786-787, doi:10.1093/glycob/cww076. The symbol
//! system itself: Varki et al. (2015) Glycobiology 25(12):1323-1324,
//! doi:10.1093/glycob/cwv091.
//!
//! # Detection
//!
//! A residue is a glycan ring candidate when its resname looks up a
//! monosaccharide in [`table::MONO_TABLE`] (or, unrecognized, when a ring
//! is still found by trying both patterns below) and at least 5 of its
//! atoms are named `C1,C2,C3,C4,C5,O5` (aldose numbering) or
//! `C2,C3,C4,C5,C6,O6` (2-ketose numbering: the anomeric carbon is C2,
//! e.g. the sialic acids). Which pattern applies is fixed by the sugar's
//! chemistry ([`is_ketose`]), not guessed from which names are present:
//! a 2-ketose like Neu5Ac still has an (exocyclic) atom named `C1`, its
//! carboxyl carbon, so counting name hits without knowing the sugar
//! would misdetect its ring as C1-anomeric and its missing `O5` as a
//! hole. The source script avoids this by running ring
//! perception (true bond connectivity) instead of a name count; fixing
//! the pattern from the looked-up sugar is the equivalent for every
//! recognized monosaccharide, and is O(residues) rather than O(atoms)
//! graph search. An unrecognized name tries both patterns (aldose
//! first) and still becomes a (white hexagon) glycan residue, matching
//! the source script's fallback for a carbohydrate ring it cannot
//! identify.
//!
//! # Orientation
//!
//! Each shape is placed at the ring centroid and built in an orthonormal
//! frame (see [`frame`]) so it faces the residue it is glycosidically
//! linked to, exactly as the source script's `linked_carb.oriented.*`
//! procs do: `u` points from the ring centroid toward the attachment
//! point, `v` is `u` crossed with the vector from the ring oxygen to the
//! centroid (so it tracks the ring's pucker), and `w = u x v`.
//!
//! # Linkage
//!
//! The source script finds an attachment by searching for an oxygen or
//! nitrogen within 1.6 A of the residue's anomeric carbon, each frame.
//! This instead walks the structure's already-perceived [`BondTable`]
//! once, at topology-build time: bonds don't change frame to frame, only
//! positions do, so this is the one-time "precompute topology" step the
//! `Glycan` representation needs to stay under a millisecond per frame.

mod mesh;
mod table;

pub use mesh::{build_linkage_mesh, build_mesh, PolytopeMesh};

use std::collections::HashMap;
use std::sync::OnceLock;

use glam::Vec3;

use crate::topology::flags;
use crate::{BondTable, Element, Topology};

/// One SNFG glyph shape. `size_factor` is the source script's per-shape
/// constant relative to its global `size` variable (`sphere_size =
/// size*0.5`, ...), so a rep's `size` option scales every shape the same
/// way the script's `snfg-enable` full/icon presets do. `FlatDiamond` is
/// the spec's separate "Flattened Diamond" (di-deoxynonulosonates: Pse,
/// Leg, Aci, 4eLeg), visually thinner than the sialic-acid `Diamond`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Shape {
    Sphere,
    Cube,
    Diamond,
    FlatDiamond,
    Cone,
    Rectangle,
    Star,
    Hexagon,
    Pentagon,
}

impl Shape {
    pub const fn size_factor(self) -> f32 {
        match self {
            Shape::Sphere => 0.5,
            Shape::Cube => 0.806,
            Shape::Diamond => 1.3,
            Shape::FlatDiamond => 1.1,
            Shape::Cone | Shape::Rectangle | Shape::Star | Shape::Pentagon => 1.0,
            Shape::Hexagon => 1.15,
        }
    }
}

/// The official SNFG palette (Table 2, Neelamegham et al. 2019
/// Glycobiology 29:620-624; also NCBI's `glycans/docs/notes.pdf`), as
/// the RGB values the spec itself publishes -- not re-derived from its
/// CMYK column, whose naive `R=(1-C)(1-K)` conversion drifts from these
/// (e.g. Blue's CMYK gives (0,128,255), not the spec's (0,144,188)).
/// Independent of `vv_scene`'s `color NAME` table.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SnfgColor {
    White,
    Blue,
    Green,
    Yellow,
    LightBlue,
    Pink,
    Purple,
    Brown,
    Orange,
    Red,
}

impl SnfgColor {
    pub const fn rgb(self) -> [u8; 3] {
        match self {
            SnfgColor::White => [255, 255, 255],
            SnfgColor::Blue => [0, 144, 188],
            SnfgColor::Green => [0, 166, 81],
            SnfgColor::Yellow => [255, 212, 0],
            SnfgColor::LightBlue => [143, 204, 233],
            SnfgColor::Pink => [246, 158, 161],
            SnfgColor::Purple => [165, 67, 153],
            SnfgColor::Brown => [161, 122, 77],
            SnfgColor::Orange => [244, 121, 32],
            SnfgColor::Red => [237, 28, 36],
        }
    }
}

/// One recognized monosaccharide: its shape, its one or two colors (two
/// for the "crossed"/"divided" shapes, e.g. white/blue for a hexosamine
/// cube), and every residue name the source script accepts for it across
/// the common/PDB, CHARMM and GLYCAM conventions.
pub struct MonoEntry {
    pub key: &'static str,
    pub label: &'static str,
    pub shape: Shape,
    pub color1: SnfgColor,
    pub color2: SnfgColor,
    pub names: &'static [&'static str],
}

/// The fallback for a ring-detected residue absent from
/// [`table::MONO_TABLE`], matching the source script's own `else` branch.
const UNKNOWN: MonoEntry = MonoEntry {
    key: "?",
    label: "Unknown glycan",
    shape: Shape::Hexagon,
    color1: SnfgColor::White,
    color2: SnfgColor::White,
    names: &[],
};

/// Case-sensitive: GLYCAM names distinguish sugars by letter case alone
/// (e.g. `YGA` is glucose, `YGa` is abequose), and the source script's
/// Tcl `lsearch` is case-sensitive too.
pub fn lookup(resname: &str) -> Option<&'static MonoEntry> {
    fn index() -> &'static HashMap<&'static str, &'static MonoEntry> {
        static INDEX: OnceLock<HashMap<&'static str, &'static MonoEntry>> = OnceLock::new();
        INDEX.get_or_init(|| {
            let mut m = HashMap::new();
            for entry in table::MONO_TABLE {
                for &name in entry.names {
                    // First match wins, as the source's elseif chain does
                    // for the one name two groups share (see table.rs).
                    m.entry(name).or_insert(entry);
                }
            }
            m
        })
    }
    index().get(resname).copied()
}

/// Whether `resname` is one of [`table::MONO_TABLE`]'s recognized names,
/// for a `select glycan` keyword: a resname-only test, unlike
/// [`GlycanPlan`], which also requires the residue to carry ring atoms.
pub fn is_glycan_name(resname: &str) -> bool {
    lookup(resname).is_some()
}

/// Ring atom names in `[C(k+1)..C(k+5), O(k+5)]` order: aldose numbering
/// (anomeric carbon C1) and 2-ketose numbering (anomeric carbon C2).
const RING_NAMES: [[&str; 6]; 2] = [
    ["C1", "C2", "C3", "C4", "C5", "O5"],
    ["C2", "C3", "C4", "C5", "C6", "O6"],
];

/// Monosaccharide keys whose anomeric carbon is C2, not C1: the
/// sialic-acid-like nonulosonic/ulosonic acids (deoxynonulosonate and
/// di-deoxynonulosonate rows: all "non-2-ulopyranosonic acid") and the
/// ketohexoses/ketopentoses (SNFG's diamond and pentagon shapes;
/// `snfg-update`'s `shiftbool`, keyed here by [`MonoEntry::key`] since
/// `table::MONO_TABLE` has no separate enum).
const KETOSE_KEYS: &[&str] = &[
    "Kdn", "Neu5Ac", "Neu5Gc", "Neu", "Sia", "Pse", "Leg", "Aci", "4eLeg", "Kdo", "Dha", "Fruc",
    "Tag", "Sor", "Psi",
];

fn is_ketose(entry: &MonoEntry) -> bool {
    KETOSE_KEYS.contains(&entry.key)
}
/// The source script's `>=5 of the 6` threshold for a partially-resolved
/// ring (missing an atom to disorder or an incomplete model) to still count.
const MIN_RING_ATOMS: usize = 5;
/// Absent ring atom (of the 6 named), for `GlycanResidue::ring_atoms`.
const NO_ATOM: u32 = u32::MAX;
/// GLYCAM's fixed Cg-Oh bond length, used to place a synthetic
/// attachment point for a residue with no resolved linkage at all.
const CG_OH_BOND: f32 = 1.43;

/// How a residue's anomeric carbon (its own glycosidic bond, going
/// toward the reducing end) connects onward.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Attachment {
    /// Another detected glycan residue, by index into `GlycanPlan::residues`.
    Residue(u32),
    /// A protein/other polymer residue's alpha carbon (N- or O-glycosidic bond).
    ProteinCa(u32),
    /// A terminal oxygen: the reducing end's free hydroxyl, or a cap
    /// (GLYCAM `OME`/`TBT`) with no further neighbor.
    Terminal(u32),
    /// No bonded neighbor at all on the anomeric carbon.
    None,
}

/// One detected glycan ring: its atoms and its place in the identified
/// monosaccharide table, precomputed once per topology.
#[derive(Clone, Copy)]
pub struct GlycanResidue {
    pub residue: u32,
    /// `[C(k+1), C(k+2), C(k+3), C(k+4), C(k+5), O(k+5)]` atom indices;
    /// `NO_ATOM` for one missing under `MIN_RING_ATOMS`'s slack.
    ring_atoms: [u32; 6],
    pub shape: Shape,
    pub color1: SnfgColor,
    pub color2: SnfgColor,
    pub label: &'static str,
    pub attachment: Attachment,
}

impl GlycanResidue {
    /// The anomeric carbon: this residue's own glycosidic bond atom.
    pub fn anomeric(&self) -> u32 {
        self.ring_atoms[0]
    }

    /// The ring oxygen (`p6` in the source script): orients the shape
    /// with the ring pucker.
    pub fn ring_oxygen(&self) -> u32 {
        self.ring_atoms[5]
    }

    pub fn ring_atom_indices(&self) -> impl Iterator<Item = u32> + '_ {
        self.ring_atoms.iter().copied().filter(|&a| a != NO_ATOM)
    }
}

/// Every glycan residue a topology contains, with its shape/color and
/// linkage resolved once. Cheap to recompute in full: rebuild whenever
/// the topology (not just coordinates) changes.
pub struct GlycanPlan {
    pub residues: Vec<GlycanResidue>,
}

/// The atoms of `residue_atoms` named by `names` (`RING_NAMES`'s aldose
/// or 2-ketose pattern), by slot; `Some` when at least [`MIN_RING_ATOMS`]
/// of the 6 are present (a little slack for a disordered or incomplete
/// model).
fn ring_atoms_named(
    topology: &Topology,
    residue_atoms: std::ops::Range<u32>,
    names: &[&str; 6],
) -> Option<[u32; 6]> {
    let mut found = [NO_ATOM; 6];
    let mut count = 0;
    for atom in residue_atoms {
        let name = topology.atom_name(atom as usize);
        if let Some(slot) = names.iter().position(|&n| n == name) {
            found[slot] = atom;
            count += 1;
        }
    }
    (count >= MIN_RING_ATOMS).then_some(found)
}

/// The residue's ring atoms: for a recognized sugar, its chemistry fixes
/// which of `RING_NAMES`'s two patterns applies ([`is_ketose`]), so an
/// exocyclic atom that happens to be named like the other pattern's
/// anomeric carbon (Neu5Ac's carboxyl `C1`) is never mistaken for it. An
/// unrecognized ring tries both, aldose first, as the source script's own
/// `C1`-before-`C2` check does when it has no name to go on either.
fn ring_atoms_of(
    topology: &Topology,
    residue_atoms: std::ops::Range<u32>,
    mono: Option<&MonoEntry>,
) -> Option<[u32; 6]> {
    match mono {
        Some(m) => ring_atoms_named(topology, residue_atoms, &RING_NAMES[is_ketose(m) as usize]),
        None => RING_NAMES
            .iter()
            .find_map(|names| ring_atoms_named(topology, residue_atoms.clone(), names)),
    }
}

fn is_hetero(topology: &Topology, atom: u32) -> bool {
    topology
        .flags
        .get(atom as usize)
        .is_some_and(|f| f & flags::HETERO != 0)
}

/// Adjacency built once from the bond table, for the linkage walk below.
fn adjacency(bonds: &BondTable, atom_count: usize) -> Vec<Vec<u32>> {
    let mut adj = vec![Vec::new(); atom_count];
    for &[a, b] in &bonds.pairs {
        adj[a as usize].push(b);
        adj[b as usize].push(a);
    }
    adj
}

fn find_atom_named(
    topology: &Topology,
    residue_atoms: std::ops::Range<u32>,
    name: &str,
) -> Option<u32> {
    residue_atoms
        .into_iter()
        .find(|&a| topology.atom_name(a as usize) == name)
}

/// The attachment reached by walking outward from `anomeric`, per the
/// source script's `snfg-update`: an oxygen bonded to it belongs to
/// either another detected ring (a glycosidic bond to that residue) or a
/// residue with a `CA` (an O-linked glycoprotein anchor); a nitrogen
/// bonded to it is an N-linked anchor (its residue's `CA`); anything
/// else is a terminal cap or the free reducing end; nothing bonded at
/// all resolves to `None`.
fn resolve_attachment(
    topology: &Topology,
    adj: &[Vec<u32>],
    anomeric: u32,
    own_residue: u32,
    ring_of_residue: &HashMap<u32, u32>,
) -> Attachment {
    let residue_of = |atom: u32| topology.residue_index[atom as usize];
    let neighbors: Vec<u32> = adj[anomeric as usize]
        .iter()
        .copied()
        .filter(|&a| residue_of(a) != own_residue)
        .collect();

    if let Some(&o_att) = neighbors
        .iter()
        .find(|&&a| topology.element[a as usize] == Element::OXYGEN)
    {
        let o_residue = residue_of(o_att);
        if let Some(&idx) = ring_of_residue.get(&o_residue) {
            return Attachment::Residue(idx);
        }
        let res = &topology.residues[o_residue as usize];
        if let Some(ca) = find_atom_named(topology, res.atoms.clone(), "CA") {
            return Attachment::ProteinCa(ca);
        }
        return Attachment::Terminal(o_att);
    }
    if let Some(&n_att) = neighbors
        .iter()
        .find(|&&a| topology.element[a as usize] == Element::NITROGEN)
    {
        let res = &topology.residues[residue_of(n_att) as usize];
        if let Some(ca) = find_atom_named(topology, res.atoms.clone(), "CA") {
            return Attachment::ProteinCa(ca);
        }
    }
    Attachment::None
}

impl GlycanPlan {
    /// Detects every glycan ring in `topology` and resolves its shape,
    /// color and linkage from `bonds` (`vv_core::bonds::perceive`, or the
    /// structure's own merged bond table). O(atoms + bonds), meant to run
    /// once per topology, not per frame.
    pub fn build(topology: &Topology, bonds: &BondTable) -> GlycanPlan {
        let mut residues = Vec::new();
        let mut ring_of_residue = HashMap::new();
        for (r, rec) in topology.residues.iter().enumerate() {
            if !is_hetero(topology, rec.atoms.start) {
                continue;
            }
            let mono = lookup(topology.residue_name(r));
            let Some(ring_atoms) = ring_atoms_of(topology, rec.atoms.clone(), mono) else {
                continue;
            };
            ring_of_residue.insert(r as u32, residues.len() as u32);
            let mono = mono.unwrap_or(&UNKNOWN);
            residues.push(GlycanResidue {
                residue: r as u32,
                ring_atoms,
                shape: mono.shape,
                color1: mono.color1,
                color2: mono.color2,
                label: mono.label,
                attachment: Attachment::None, // resolved below, once every ring is known
            });
        }

        let adj = adjacency(bonds, topology.atom_count());
        for res in &mut residues {
            res.attachment = resolve_attachment(
                topology,
                &adj,
                res.anomeric(),
                res.residue,
                &ring_of_residue,
            );
        }
        GlycanPlan { residues }
    }
}

/// Per-frame placement: ring centroid, ring-oxygen position (orients the
/// shape with the ring pucker) and attachment point (orients it toward
/// its glycosidic bond), one entry per `GlycanPlan::residues`.
#[derive(Default)]
pub struct GlycanFrame {
    pub centroid: Vec<Vec3>,
    pub ring_oxygen: Vec<Vec3>,
    pub attach_point: Vec<Vec3>,
}

impl GlycanPlan {
    /// Recomputes `out` for `positions` (the structure's current frame).
    /// Reuses `out`'s buffers instead of allocating: after the plan's
    /// residue count first sizes them, replaying a trajectory allocates
    /// nothing here.
    pub fn update_into(&self, positions: &[Vec3], out: &mut GlycanFrame) {
        let n = self.residues.len();
        out.centroid.resize(n, Vec3::ZERO);
        out.ring_oxygen.resize(n, Vec3::ZERO);
        out.attach_point.resize(n, Vec3::ZERO);

        for (i, res) in self.residues.iter().enumerate() {
            let mut sum = Vec3::ZERO;
            let mut count = 0.0f32;
            for a in res.ring_atom_indices() {
                sum += positions[a as usize];
                count += 1.0;
            }
            out.centroid[i] = sum / count.max(1.0);
            out.ring_oxygen[i] = positions[res.ring_oxygen() as usize];
        }
        // A `Residue` attachment reads another ring's centroid, so this
        // pass runs after every centroid above is filled.
        for (i, res) in self.residues.iter().enumerate() {
            out.attach_point[i] = match res.attachment {
                Attachment::Residue(j) => out.centroid[j as usize],
                Attachment::ProteinCa(a) | Attachment::Terminal(a) => positions[a as usize],
                Attachment::None => {
                    let p1 = positions[res.anomeric() as usize];
                    let dir = (p1 - out.centroid[i]).normalize_or_zero();
                    p1 + dir * CG_OH_BOND
                }
            };
        }
    }
}

/// The orthonormal frame a shape is drawn in: `u` toward the attachment
/// (the source script's `vec_AB`, normalized), `v` derived from the ring
/// oxygen so the shape tracks the ring pucker, `w = u x v`. Degenerate
/// input (a zero-length `u`, or `v` colinear with it) falls back to an
/// arbitrary right-handed frame rather than producing NaNs.
pub fn frame(centroid: Vec3, attach_point: Vec3, ring_oxygen: Vec3) -> (Vec3, Vec3, Vec3) {
    let u = (attach_point - centroid).normalize_or(Vec3::X);
    let mut v = u.cross(centroid - ring_oxygen);
    if v.length_squared() < 1e-8 {
        v = u.cross(Vec3::Y);
        if v.length_squared() < 1e-8 {
            v = u.cross(Vec3::Z);
        }
    }
    let v = v.normalize();
    let w = u.cross(v);
    (u, v, w)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bonds::perceive;
    use crate::intern::Interner;
    use crate::topology::{ChainRec, ResidueRec, SecondaryStructure};

    #[test]
    fn looks_up_every_naming_convention_for_glcnac() {
        for name in ["NAG", "NDG", "AGLCNA", "BGLCNA", "0YA", "0yB"] {
            let e = lookup(name).unwrap_or_else(|| panic!("{name} not recognized"));
            assert_eq!(e.shape, Shape::Cube);
            assert_eq!(e.color1, SnfgColor::Blue);
        }
    }

    #[test]
    fn man_and_bma_are_the_green_sphere() {
        for name in ["MAN", "BMA", "AMAN"] {
            let e = lookup(name).unwrap();
            assert_eq!(e.shape, Shape::Sphere);
            assert_eq!(e.color1, SnfgColor::Green);
        }
    }

    #[test]
    fn case_sensitive_glycam_codes_pick_different_sugars() {
        // Glucose (YGA) and abequose (YGa) differ only in case, as GLYCAM
        // intends; the source script's Tcl `lsearch` is case-sensitive.
        assert_eq!(lookup("YGA").unwrap().shape, Shape::Sphere); // Glc
        assert_eq!(lookup("YGa").unwrap().shape, Shape::Rectangle); // Abe
    }

    #[test]
    fn shared_code_prefers_the_first_matching_group() {
        // 4YS is listed under both GlcNAc (tested first in the source's
        // elseif chain) and GlcN; GlcNAc must win.
        let e = lookup("4YS").unwrap();
        assert_eq!(e.shape, Shape::Cube);
        assert_eq!(e.color1, SnfgColor::Blue);
        assert_eq!(e.color2, SnfgColor::Blue);
    }

    #[test]
    fn unknown_name_is_not_found() {
        assert!(lookup("XYZZY").is_none());
    }

    /// Two linked NAG rings: NAG1-O4-C1-NAG2's C1 anomeric carbon, so
    /// residue 0 (the reducing end) has no outgoing attachment and
    /// residue 1's anomeric carbon links back to residue 0's ring.
    /// Coordinates are a rough hexagon, not a real pyranose geometry:
    /// only bond perception and the ring-atom-name detection are under
    /// test here, not real ring chemistry.
    fn two_linked_rings() -> (Topology, Vec<Vec3>) {
        let mut names = Interner::new();
        let nag = names.intern("NAG");
        let a = names.intern("A");
        let ring_names = ["C1", "C2", "C3", "C4", "C5", "O5", "O4"];
        let mut name = Vec::new();
        let mut element = Vec::new();
        let mut flags_ = Vec::new();
        let mut positions = Vec::new();
        // Ring detection is name-only, so these atoms need not really be
        // bonded into a ring: keep every atom far from every other one
        // (no spurious bonds) except the one pair that must perceive as
        // a bond -- ring 0's O4 to ring 1's C1, the glycosidic linkage.
        for ring in 0..2 {
            for (i, n) in ring_names.iter().enumerate() {
                let mut buf = [b' '; 4];
                buf[..n.len()].copy_from_slice(n.as_bytes());
                name.push(buf);
                element.push(if n.starts_with('O') {
                    Element::OXYGEN
                } else {
                    Element::CARBON
                });
                flags_.push(flags::HETERO);
                let slot = (ring * ring_names.len() + i) as f32;
                positions.push(Vec3::new(slot * 20.0, 0.0, 0.0));
            }
        }
        let atoms_per_ring = ring_names.len() as u32;
        let topology = Topology {
            element,
            name,
            flags: flags_,
            residue_index: (0..2u32)
                .flat_map(|r| std::iter::repeat_n(r, ring_names.len()))
                .collect(),
            residues: (0..2)
                .map(|r| ResidueRec {
                    atoms: r as u32 * atoms_per_ring..(r as u32 + 1) * atoms_per_ring,
                    chain: 0,
                    comp: nag,
                    seq_id: r + 1,
                    auth_seq_id: r + 1,
                    ins_code: 0,
                    ss: SecondaryStructure::Unknown,
                })
                .collect(),
            chains: vec![ChainRec {
                residues: 0..2,
                label_asym: a,
                auth_asym: a,
                entity: 1,
            }],
            names,
            ..Default::default()
        };
        // Move ring 1's C1 (index 7) next to ring 0's O4 (index 6) so
        // distance-based bond perception links them 1->4.
        positions[7] = positions[6] + Vec3::new(1.3, 0.0, 0.0);
        (topology, positions)
    }

    #[test]
    fn detects_both_rings_and_links_the_second_to_the_first() {
        let (topology, positions) = two_linked_rings();
        let bonds = perceive(&topology, &positions);
        let plan = GlycanPlan::build(&topology, &bonds);
        assert_eq!(plan.residues.len(), 2);
        assert_eq!(plan.residues[0].attachment, Attachment::None);
        assert_eq!(plan.residues[1].attachment, Attachment::Residue(0));
    }

    #[test]
    fn per_frame_update_reuses_its_buffers() {
        let (topology, positions) = two_linked_rings();
        let bonds = perceive(&topology, &positions);
        let plan = GlycanPlan::build(&topology, &bonds);
        let mut frame_out = GlycanFrame::default();
        plan.update_into(&positions, &mut frame_out);
        let cap = frame_out.centroid.capacity();
        plan.update_into(&positions, &mut frame_out);
        assert_eq!(
            frame_out.centroid.capacity(),
            cap,
            "resize should be a no-op the 2nd time"
        );
        assert_eq!(frame_out.centroid.len(), 2);
    }

    /// A single SIA (Neu5Ac) residue: ring atoms C2-C6,O6 plus an
    /// exocyclic C1 (its carboxyl carbon, not in the ring) and O5, which
    /// is NOT one of its atoms (regression test for the ring-shift bug:
    /// counting raw name hits without knowing the sugar finds C1..C5 -- 5
    /// names -- and wrongly picks aldose numbering with a missing O5).
    fn lone_sialic_acid_residue() -> (Topology, Vec<Vec3>) {
        let mut names = Interner::new();
        let sia = names.intern("SIA");
        let a = names.intern("A");
        let atom_names = ["C1", "C2", "C3", "C4", "C5", "C6", "O6"];
        let mut name = Vec::new();
        let mut element = Vec::new();
        let mut flags_ = Vec::new();
        let mut positions = Vec::new();
        for (i, n) in atom_names.iter().enumerate() {
            let mut buf = [b' '; 4];
            buf[..n.len()].copy_from_slice(n.as_bytes());
            name.push(buf);
            element.push(if n.starts_with('O') {
                Element::OXYGEN
            } else {
                Element::CARBON
            });
            flags_.push(flags::HETERO);
            positions.push(Vec3::new(i as f32 * 20.0, 0.0, 0.0));
        }
        let topology = Topology {
            element,
            name,
            flags: flags_,
            residue_index: vec![0; atom_names.len()],
            residues: vec![ResidueRec {
                atoms: 0..atom_names.len() as u32,
                chain: 0,
                comp: sia,
                seq_id: 1,
                auth_seq_id: 1,
                ins_code: 0,
                ss: SecondaryStructure::Unknown,
            }],
            chains: vec![ChainRec {
                residues: 0..1,
                label_asym: a,
                auth_asym: a,
                entity: 1,
            }],
            names,
            ..Default::default()
        };
        (topology, positions)
    }

    #[test]
    fn sialic_acid_ring_excludes_its_exocyclic_c1() {
        let (topology, positions) = lone_sialic_acid_residue();
        let bonds = perceive(&topology, &positions); // atoms are far apart: no bonds, and none needed
        let plan = GlycanPlan::build(&topology, &bonds);
        assert_eq!(plan.residues.len(), 1);
        let res = &plan.residues[0];
        assert_eq!(res.label, "N-Acetylneuraminic acid (purple diamond)");
        // Anomeric carbon is C2 (index 1), ring oxygen is O6 (index 6);
        // the exocyclic C1 (index 0) is in neither slot.
        assert_eq!(res.anomeric(), 1);
        assert_eq!(res.ring_oxygen(), 6);
        assert!(
            res.ring_atom_indices().all(|a| a != 0),
            "C1 must not be a ring atom"
        );

        let mut out = GlycanFrame::default();
        plan.update_into(&positions, &mut out);
        // Centroid over C2..C6,O6 (indices 1..=6), excluding C1 (index 0).
        let expected: Vec3 = (1..=6).map(|i| positions[i]).sum::<Vec3>() / 6.0;
        assert!((out.centroid[0] - expected).length() < 1e-4);
    }

    #[test]
    fn frame_is_orthonormal() {
        let (u, v, w) = frame(Vec3::ZERO, Vec3::X * 2.0, Vec3::new(0.0, 1.0, 0.0));
        assert!((u.length() - 1.0).abs() < 1e-5);
        assert!((v.length() - 1.0).abs() < 1e-5);
        assert!((w.length() - 1.0).abs() < 1e-5);
        assert!(u.dot(v).abs() < 1e-5);
        assert!(u.dot(w).abs() < 1e-5);
        assert!(v.dot(w).abs() < 1e-5);
    }
}
