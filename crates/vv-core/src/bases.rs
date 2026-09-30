//! Nucleotide bases for the cartoon: each base as a slab on its own ring
//! atoms, and base-pair detection for the one-rung-per-pair "ladder".
//!
//! A slab is the ring perimeter (purine: the nine atoms of the fused
//! rings; pyrimidine: six) fan-triangulated from the ring centroid and
//! extruded either side of the local plane, so a puckered or modified
//! base follows its actual atoms. Pairs are found from geometry alone:
//! a purine and a pyrimidine whose Watson-Crick edges (N1 of the purine,
//! N3 of the pyrimidine) are within hydrogen-bond distance, or that form
//! the wobble contacts (purine N1 to pyrimidine O2 and purine O6 to
//! pyrimidine N3), with roughly parallel, co-planar rings.

use glam::Vec3;

use crate::glycan::PolytopeMesh;
use crate::residue_class::ResidueClass;
use crate::topology::Topology;

/// Perimeter order, so consecutive atoms are bonded.
const PURINE_RING: [&str; 9] = ["N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"];
const PYRIMIDINE_RING: [&str; 6] = ["N1", "C2", "N3", "C4", "C5", "C6"];

/// Half the slab's thickness (Angstrom): thin enough to read as a plate,
/// thick enough to show an edge from the side.
pub const HALF_THICKNESS: f32 = 0.2;
/// How far the slab's outline sits beyond the ring atoms' centers, so the
/// plate covers the atoms' extent rather than stopping at their centers.
pub const MARGIN: f32 = 0.35;
/// Longest hydrogen-bond heavy-atom distance for a base-pair contact
/// (2.8 to 3.0 A in Watson-Crick pairs, with slack for crystal noise).
pub const PAIR_DISTANCE: f32 = 3.4;
/// Largest angle between the two rings' normals (propeller twist and
/// buckle stay well inside it; a stacked or perpendicular base does not).
pub const PAIR_MAX_TILT_DEG: f32 = 40.0;
/// Largest distance of one ring's centroid from the other's plane.
pub const PAIR_MAX_STAGGER: f32 = 2.0;
/// A glycosidic bond is at most this long (C1' to the base's ring atom).
const GLYCOSIDIC_MAX: f32 = 2.2;

/// One nucleotide's base.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Base {
    pub residue: u32,
    pub purine: bool,
    /// Ring atoms in perimeter order.
    pub ring: Vec<u32>,
    /// The sugar's `C1'` when it is bonded to the base; else the
    /// glycosidic atom itself.
    pub anchor: u32,
    /// The ring atom bonded to `anchor` (N9 / N1, or C5 in pseudouridine).
    pub glycosidic: u32,
    /// The backbone trace atom (`P`, else `C4'`) the base hangs from.
    pub trace: Option<u32>,
    /// N1 (purine) or N3 (pyrimidine): the Watson-Crick edge nitrogen.
    n_edge: u32,
    /// O6 (purine) or O2 (pyrimidine), when present: the wobble partner.
    o_edge: Option<u32>,
}

/// Which atoms may stand for a named base atom: the shown conformer's.
type Visible<'a> = &'a dyn Fn(u32) -> bool;

fn atom_named(topology: &Topology, residue: usize, name: &str, visible: Visible) -> Option<u32> {
    topology.residues[residue]
        .atoms
        .clone()
        .find(|&a| visible(a) && topology.atom_name(a as usize) == name)
}

fn ring_atoms(
    topology: &Topology,
    residue: usize,
    names: &[&str],
    visible: Visible,
) -> Option<Vec<u32>> {
    names
        .iter()
        .map(|n| atom_named(topology, residue, n, visible))
        .collect()
}

fn base_of(
    topology: &Topology,
    positions: &[Vec3],
    residue: usize,
    visible: Visible,
) -> Option<Base> {
    let named = |name| atom_named(topology, residue, name, visible);
    let purine = named("N9").is_some();
    let names: &[&str] = if purine {
        &PURINE_RING
    } else {
        &PYRIMIDINE_RING
    };
    let ring = ring_atoms(topology, residue, names, visible)?;
    let c1 = named("C1'");
    let glycosidic = match c1 {
        Some(c1) => *ring.iter().min_by(|&&a, &&b| {
            let d = |x: u32| positions[x as usize].distance(positions[c1 as usize]);
            d(a).total_cmp(&d(b))
        })?,
        None => ring[0],
    };
    let anchor = c1
        .filter(|&c| {
            positions[c as usize].distance(positions[glycosidic as usize]) <= GLYCOSIDIC_MAX
        })
        .unwrap_or(glycosidic);
    let (n_name, o_name) = if purine { ("N1", "O6") } else { ("N3", "O2") };
    Some(Base {
        residue: residue as u32,
        purine,
        ring,
        anchor,
        glycosidic,
        trace: crate::backbone::trace_atom(topology, residue).map(|(a, _)| a),
        n_edge: named(n_name)?,
        o_edge: named(o_name),
    })
}

/// Every nucleotide whose whole base ring is present among the `visible`
/// atoms, in residue order. `positions` picks the glycosidic atom (the
/// ring atom nearest `C1'`).
pub fn bases(topology: &Topology, positions: &[Vec3], visible: Visible) -> Vec<Base> {
    (0..topology.residues.len())
        .filter(|&r| topology.residue_class(r) == ResidueClass::Nucleic)
        .filter_map(|r| base_of(topology, positions, r, visible))
        .collect()
}

/// Unit normal of the ring's mean plane (Newell's method: right-handed
/// around the perimeter order), and the ring's centroid.
fn frame(ring: &[u32], positions: &[Vec3]) -> (Vec3, Vec3) {
    let points: Vec<Vec3> = ring.iter().map(|&a| positions[a as usize]).collect();
    let centroid = points.iter().copied().sum::<Vec3>() / points.len() as f32;
    let normal = (0..points.len())
        .map(|i| points[i].cross(points[(i + 1) % points.len()]))
        .sum::<Vec3>()
        .normalize_or_zero();
    (normal, centroid)
}

/// A slab's outline point for ring atom `p`: pushed out of the centroid
/// within the ring plane by [`MARGIN`].
fn outline(p: Vec3, centroid: Vec3, normal: Vec3) -> Vec3 {
    let radial = p - centroid;
    let flat = (radial - normal * radial.dot(normal)).normalize_or_zero();
    p + flat * MARGIN
}

/// Appends `base`'s slab to `out`, every vertex tagged with the ring atom
/// it sits on (the centroid with the glycosidic atom) so a pick names a
/// real atom of the base.
pub fn push_plate(base: &Base, positions: &[Vec3], color: [u8; 3], out: &mut PolytopeMesh) {
    let (normal, centroid) = frame(&base.ring, positions);
    let up = normal * HALF_THICKNESS;
    let rim: Vec<Vec3> = base
        .ring
        .iter()
        .map(|&a| outline(positions[a as usize], centroid, normal))
        .collect();
    let n = rim.len();
    for i in 0..n {
        let j = (i + 1) % n;
        let (a, b) = (base.ring[i], base.ring[j]);
        let (t0, t1, b0, b1) = (rim[i] + up, rim[j] + up, rim[i] - up, rim[j] - up);
        out.tri(centroid + up, t0, t1, color, base.glycosidic);
        out.tri(centroid - up, b1, b0, color, base.glycosidic);
        out.tri(b0, b1, t1, color, a);
        out.tri(b0, t1, t0, color, b);
    }
}

/// Each base's stick from the backbone trace atom to its sugar.
fn stems(bases: &[Base]) -> Vec<[u32; 2]> {
    bases
        .iter()
        .filter_map(|b| b.trace.map(|t| [t, b.anchor]))
        .filter(|[a, b]| a != b)
        .collect()
}

/// Each base's sticks from the backbone trace atom to its sugar and on
/// to the glycosidic atom.
pub fn connectors(bases: &[Base]) -> Vec<[u32; 2]> {
    let links = bases
        .iter()
        .filter(|b| b.anchor != b.glycosidic)
        .map(|b| [b.anchor, b.glycosidic]);
    stems(bases).into_iter().chain(links).collect()
}

/// Two paired bases, as indices into the `bases` they came from.
pub type BasePair = [usize; 2];

fn ring_frames(bases: &[Base], positions: &[Vec3]) -> Vec<(Vec3, Vec3)> {
    bases.iter().map(|b| frame(&b.ring, positions)).collect()
}

/// How closely `pu` (purine) and `py` (pyrimidine) meet at the edges, as
/// the longest of the contacts that make up their best pairing mode, or
/// `None` when neither Watson-Crick nor wobble reaches [`PAIR_DISTANCE`].
fn contact(pu: &Base, py: &Base, positions: &[Vec3]) -> Option<f32> {
    let d = |a: u32, b: u32| positions[a as usize].distance(positions[b as usize]);
    let watson = d(pu.n_edge, py.n_edge);
    let wobble = match (pu.o_edge, py.o_edge) {
        (Some(o6), Some(o2)) => Some(d(pu.n_edge, o2).max(d(o6, py.n_edge))),
        _ => None,
    };
    [Some(watson), wobble]
        .into_iter()
        .flatten()
        .filter(|&x| x <= PAIR_DISTANCE)
        .min_by(f32::total_cmp)
}

fn coplanar(fa: (Vec3, Vec3), fb: (Vec3, Vec3)) -> bool {
    let (na, ca) = fa;
    let (nb, cb) = fb;
    na.dot(nb).abs() >= PAIR_MAX_TILT_DEG.to_radians().cos()
        && (cb - ca).dot(na).abs() <= PAIR_MAX_STAGGER
        && (ca - cb).dot(nb).abs() <= PAIR_MAX_STAGGER
}

/// Base pairs among `bases`, each base in at most one: the closest
/// contacts claim their bases first. Both orders of a pair are not
/// repeated; the purine comes first.
pub fn pairs(bases: &[Base], positions: &[Vec3]) -> Vec<BasePair> {
    let frames = ring_frames(bases, positions);
    let edge_reach = PAIR_DISTANCE + 2.0;
    let mut found: Vec<(f32, BasePair)> = Vec::new();
    for (i, pu) in bases.iter().enumerate().filter(|(_, b)| b.purine) {
        for (j, py) in bases.iter().enumerate().filter(|(_, b)| !b.purine) {
            if positions[pu.n_edge as usize].distance(positions[py.n_edge as usize]) > edge_reach {
                continue;
            }
            if let Some(d) = contact(pu, py, positions).filter(|_| coplanar(frames[i], frames[j])) {
                found.push((d, [i, j]));
            }
        }
    }
    found.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut taken = vec![false; bases.len()];
    let mut out = Vec::new();
    for (_, [i, j]) in found {
        if !taken[i] && !taken[j] {
            taken[i] = true;
            taken[j] = true;
            out.push([i, j]);
        }
    }
    out.sort_unstable();
    out
}

/// A stem from the backbone for every base, one rung per pair (sugar to
/// sugar), and a stub to the glycosidic atom for every unpaired base.
pub fn ladder(bases: &[Base], pairs: &[BasePair]) -> Vec<[u32; 2]> {
    let mut paired = vec![false; bases.len()];
    let mut rungs = stems(bases);
    for &[i, j] in pairs {
        paired[i] = true;
        paired[j] = true;
        rungs.push([bases[i].anchor, bases[j].anchor]);
    }
    let stubs = bases
        .iter()
        .zip(&paired)
        .filter(|(b, &p)| !p && b.anchor != b.glycosidic)
        .map(|(b, _)| [b.anchor, b.glycosidic]);
    rungs.extend(stubs);
    rungs
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::{AtomRow, TopologyBuilder};

    /// A flat regular hexagon (a stand-in pyrimidine ring) at height `z`,
    /// 1.4 A from its center, with C1' beside the first atom.
    fn pyrimidine(z: f32) -> (Vec<(&'static str, Vec3)>, Vec3) {
        let center = Vec3::new(0.0, 0.0, z);
        let ring = PYRIMIDINE_RING
            .iter()
            .enumerate()
            .map(|(i, name)| {
                let t = std::f32::consts::TAU * i as f32 / 6.0;
                (*name, center + Vec3::new(1.4 * t.cos(), 1.4 * t.sin(), 0.0))
            })
            .collect();
        (ring, center)
    }

    fn topology(atoms: &[(&str, Vec3)]) -> (Topology, Vec<Vec3>) {
        let mut b = TopologyBuilder::new();
        for (i, (atom, position)) in atoms.iter().enumerate() {
            let mut name = [b' '; 4];
            name[..atom.len()].copy_from_slice(atom.as_bytes());
            b.push(&AtomRow {
                element: crate::Element::CARBON,
                name,
                serial: i as u32 + 1,
                alt_loc: 0,
                comp: "DC",
                asym: "A",
                auth_asym: "A",
                seq_id: 1,
                auth_seq_id: 1,
                ins_code: 0,
                entity: 1,
                position: *position,
                occupancy: 1.0,
                b_factor: 0.0,
                charge: 0,
                hetero: false,
            });
        }
        let structure = b.finish().unwrap();
        let positions = structure.frame(0).positions().to_vec();
        (
            std::sync::Arc::try_unwrap(structure.topology).unwrap(),
            positions,
        )
    }

    #[test]
    fn a_plate_sits_half_a_thickness_either_side_of_the_ring_plane() {
        let (mut atoms, _) = pyrimidine(3.0);
        atoms.push(("C1'", Vec3::new(1.4 + 1.5, 0.0, 3.0)));
        let (t, positions) = topology(&atoms);
        let found = bases(&t, &positions, &|_| true);
        assert_eq!(found.len(), 1);
        let mut mesh = PolytopeMesh::default();
        push_plate(&found[0], &positions, [1, 2, 3], &mut mesh);
        assert_eq!(mesh.triangle_count(), 6 * 4);
        for p in &mesh.positions {
            assert!(((p.z - 3.0).abs() - HALF_THICKNESS).abs() < 1e-4, "{p}");
        }
        for (n, atom) in mesh.normals.iter().zip(&mesh.source_atom) {
            assert!(n.z.abs() > 0.99 || n.z.abs() < 1e-4, "{n} {atom}");
        }
    }

    #[test]
    fn a_purine_has_nine_ring_atoms_and_a_pyrimidine_six() {
        assert_eq!(PURINE_RING.len(), 9);
        assert_eq!(PYRIMIDINE_RING.len(), 6);
        assert!(PURINE_RING.contains(&"N7") && !PYRIMIDINE_RING.contains(&"N7"));
    }

    #[test]
    fn a_base_missing_a_ring_atom_is_skipped() {
        let (mut atoms, _) = pyrimidine(0.0);
        atoms.pop();
        let (t, positions) = topology(&atoms);
        assert!(bases(&t, &positions, &|_| true).is_empty());
    }
}
