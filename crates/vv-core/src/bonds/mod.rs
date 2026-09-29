//! Bond table and how bonds are decided.
//!
//! [`perceive`] combines, in order:
//! 1. **Standard-residue templates** ([`templates`]): a residue whose name
//!    the PDB Chemical Component Dictionary covers (amino acids,
//!    nucleotides, water, common glycans, heme cofactors) bonds exactly
//!    per that residue's CCD atom graph, not by distance. This is what
//!    keeps a tight crystal-packing clash or a compressed MD frame from
//!    growing a false bond between two atoms of the same residue: real
//!    chemistry decides, geometry doesn't get a vote.
//! 2. **Named inter-residue links**: peptide C(i)-N(i+1), nucleic
//!    O3'(i)-P(i+1), disulfide SG-SG, and glycosidic links (glycan C1/C2
//!    to an acceptor oxygen or an N-glycosylation ASN's ND2), each within
//!    its own chemistry-specific distance window. This is the only way
//!    two different residues of a polymer end up bonded.
//! 3. **Distance-based perception** (the original algorithm) for
//!    everything a template doesn't cover -- ligands, unknown residues,
//!    modified residues -- now with a per-element valence cap and a
//!    clash guard, and never applied to a metal atom or a single-atom
//!    (ion) residue: those bond only through a cofactor template's own
//!    coordination bonds (heme's Fe-N) or an explicit file bond, never by
//!    proximity. That is what stops a metal from bonding a nearby water
//!    or an unrelated close atom.
//! 4. **Explicit file bonds** (`CONECT`, `_struct_conn`, `LINK`), always
//!    merged in verbatim and never removed by any of the above.
//!
//! If `topology.md_bonds` is `Some` (an MD topology file supplied its own
//! bond list -- PSF, PRMTOP), perception is skipped entirely and those
//! bonds are used as-is.

use std::ops::Range;
use std::sync::atomic::{AtomicU32, Ordering};

use glam::Vec3;
use rayon::prelude::*;

use crate::{Element, ExplicitBond, ExplicitBondKind, ResidueRec, Topology};

pub mod geometry;
pub mod templates;

pub use geometry::{
    adjacency, bond_strand_endpoints, bond_strands, strand_axis, Adjacency, BondStrand,
};

/// A bond's Kekulé order. Aromatic rings in the wwPDB Chemical Component
/// Dictionary are already Kekulized (alternating `SING`/`DOUB` with
/// `pdbx_aromatic_flag = Y`), so that alternation is what templates and
/// `_chem_comp_bond` parsing store directly as `Single`/`Double` -- no
/// separate aromatic geometry is needed for those. `Aromatic` exists only
/// for the rare literal `AROM`/`DELO` (undetermined, non-Kekulé)
/// `value_order` some depositor-submitted ligands use; absent ring
/// membership to place a single inner line correctly, it draws as a
/// `Double` (two parallel strands) -- the closest unambiguous
/// approximation without ring detection.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[repr(u8)]
pub enum BondOrder {
    #[default]
    Single = 0,
    Double = 1,
    Triple = 2,
    Aromatic = 3,
}

impl BondOrder {
    /// Parallel cylinders/lines this order draws as.
    pub fn strand_count(self) -> u32 {
        match self {
            BondOrder::Single => 1,
            BondOrder::Double | BondOrder::Aromatic => 2,
            BondOrder::Triple => 3,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BondTable {
    /// Each pair is stored once with `a < b`.
    pub pairs: Vec<[u32; 2]>,
    /// `(bond index into `pairs`, order)` for every bond that is not
    /// `Single`, sorted by index; a bond absent here is `Single`. Empty in
    /// the common all-single case (no MD topology carries order data at
    /// all, see `perceive`'s module doc), so it costs nothing there.
    pub orders: Vec<(u32, BondOrder)>,
}

impl BondTable {
    pub fn len(&self) -> usize {
        self.pairs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }

    pub fn contains(&self, a: u32, b: u32) -> bool {
        let key = if a < b { [a, b] } else { [b, a] };
        self.pairs.binary_search(&key).is_ok()
    }

    /// `bond`'s order, or `Single` if it isn't in the sparse `orders` list.
    pub fn order_of(&self, bond: u32) -> BondOrder {
        self.orders
            .binary_search_by_key(&bond, |&(i, _)| i)
            .map_or(BondOrder::Single, |k| self.orders[k].1)
    }

    /// Number of bonds per atom.
    pub fn degrees(&self, atom_count: usize) -> Vec<u32> {
        let mut degree = vec![0u32; atom_count];
        for [a, b] in &self.pairs {
            degree[*a as usize] += 1;
            degree[*b as usize] += 1;
        }
        degree
    }

    /// Connected component of the bond graph each atom belongs to (a
    /// "fragment"): dense 0-based ids, assigned in atom order so the same
    /// structure numbers its fragments the same way every time (union-find
    /// with path halving, effectively linear in bond count).
    pub fn fragments(&self, atom_count: usize) -> Vec<u32> {
        fn find(parent: &mut [u32], mut x: u32) -> u32 {
            while parent[x as usize] != x {
                parent[x as usize] = parent[parent[x as usize] as usize];
                x = parent[x as usize];
            }
            x
        }
        let mut parent: Vec<u32> = (0..atom_count as u32).collect();
        for &[a, b] in &self.pairs {
            let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
            if ra != rb {
                parent[ra as usize] = rb;
            }
        }
        let mut id_of_root = vec![u32::MAX; atom_count];
        let mut next_id = 0u32;
        (0..atom_count as u32)
            .map(|a| {
                let root = find(&mut parent, a) as usize;
                if id_of_root[root] == u32::MAX {
                    id_of_root[root] = next_id;
                    next_id += 1;
                }
                id_of_root[root]
            })
            .collect()
    }
}

/// Distance slack added to the sum of covalent radii, for the generic
/// (untemplated) distance-based path.
pub const COVALENT_TOLERANCE: f32 = 0.45;
/// A generic candidate pair closer than this fraction of the sum of
/// covalent radii is a clash, not a bond (task: "much shorter than
/// expected for that element pair").
const CLASH_FACTOR: f32 = 0.6;
/// Anything beyond this can never bond; also the spatial grid cell size.
const MAX_BOND: f32 = 3.0;

/// A template pair beyond this distance is not trusted even though the
/// names match: guards against a residue that is merely mislabeled (its
/// coordinates belong to something else). Generous relative to any real
/// covalent bond in `templates` (heaviest is ~1.9 A, C-S) so ordinary
/// structural noise or a strained model never loses a real bond to this
/// cap.
const TEMPLATE_HEAVY_MAX: f32 = 2.2;
const TEMPLATE_H_MAX: f32 = 1.3;
/// A template pair where one atom is a metal is coordination, not a
/// covalent bond, and runs longer: deoxyhemoglobin's high-spin Fe-N is
/// ~2.1-2.3 A. Matches the old flat metal cutoff this module used to
/// apply everywhere (now only here and nowhere else -- see module doc).
const TEMPLATE_METAL_MAX: f32 = 2.8;
/// A templated residue whose named atoms mostly fail to connect (wrong
/// naming convention, or really a different residue under a reused name)
/// falls back to distance-based perception for all its atoms instead of
/// silently leaving most of it unbonded.
const TEMPLATE_HEAVY_COVERAGE: f32 = 0.9;

/// Peptide C(i)-N(i+1): real length ~1.33 A (task: "~1.3-1.5 A + tolerance").
const PEPTIDE_BOND_MAX: f32 = 1.6;
/// Nucleic O3'(i)-P(i+1): real P-O ester length ~1.6 A.
const NUCLEIC_LINK_MAX: f32 = 1.75;
/// Disulfide SG-SG (task-specified).
const DISULFIDE_MAX: f32 = 2.3;
/// Glycosidic C1/C2-O/N: real length ~1.4-1.5 A (6X3Z's own LINK records
/// report 1.43-1.49 A).
const GLYCOSIDIC_MAX: f32 = 1.6;
/// Below this, two named-link atoms are the same site (guards against a
/// zero-length or duplicate-coordinate degenerate input), not a bond.
const LINK_MIN: f32 = 0.5;

/// Max simultaneous bonds per element (task: "C 4, N 4, O 2, S 6 ...
/// priority to shortest"), for the generic distance-based path only:
/// standard-residue templates and named links are chemistry, not
/// geometry, and are never capped.
fn valence_cap(e: Element) -> u32 {
    match e.atomic_number() {
        1 => 1,                // H
        5 | 14 => 4,           // B, Si
        6 => 4,                // C
        7 => 4,                // N
        8 => 2,                // O
        9 | 17 | 35 | 53 => 1, // F, Cl, Br, I
        15 | 16 => 6,          // P, S (hypervalent oxoanions)
        _ => 6,
    }
}

/// Geometric bond perception, combined with standard-residue chemistry;
/// see the module doc for the full pipeline. Runs in parallel over a
/// uniform grid for the distance-based fallback.
pub fn perceive(topology: &Topology, positions: &[Vec3]) -> BondTable {
    let n = positions.len();
    assert_eq!(n, topology.atom_count());
    if n == 0 {
        return BondTable::default();
    }
    if let Some(md_bonds) = &topology.md_bonds {
        return verbatim(md_bonds, n);
    }

    let element = &topology.element;
    let alt = &topology.alt_loc;

    let templated = match_templates(topology, positions);
    let mut pairs = link_glycosidic(&templated, positions, topology);
    pairs.extend(templated.pairs);
    pairs.extend(link_adjacent_residues(topology, positions));
    pairs.extend(link_disulfides(topology, positions));
    pairs.extend(generic_distance_bonds(
        positions,
        element,
        alt,
        &templated.is_templated,
        topology,
    ));

    // Priority low to high (see module doc #1-3): a legacy PDB CONECT
    // repeat (0), the built-in residue template (1), the file's own
    // `_chem_comp_bond` (2). `resolve_orders` keeps the highest priority
    // per pair.
    let mut order_overrides: Vec<([u32; 2], u8, BondOrder)> = templated
        .orders
        .iter()
        .map(|&(pair, ord)| (pair, 1, ord))
        .collect();

    for bond in &topology.explicit_bonds {
        let [a, b] = bond.atoms;
        if a != b && (a as usize) < n && (b as usize) < n {
            let pair = order(a, b);
            pairs.push(pair);
            if bond.order != BondOrder::Single {
                order_overrides.push((pair, 0, bond.order));
            }
        }
    }
    order_overrides.extend(
        topology
            .chem_comp_bond_order
            .iter()
            .map(|&(pair, ord)| (pair, 2, ord)),
    );

    pairs.par_sort_unstable();
    pairs.dedup();
    let orders = resolve_orders(order_overrides, &pairs);
    BondTable { pairs, orders }
}

/// Keeps the highest-priority order per pair (see `perceive`), then
/// resolves each into its index in the final, already sorted and
/// deduplicated `pairs` -- naturally still sorted by index, since the
/// overrides stay sorted by pair while `pairs.binary_search` walks the
/// same order. A pair an override names that perception never actually
/// bonded (a mismatched `_chem_comp_bond` reference) is dropped rather
/// than invented.
fn resolve_orders(
    mut overrides: Vec<([u32; 2], u8, BondOrder)>,
    pairs: &[[u32; 2]],
) -> Vec<(u32, BondOrder)> {
    if overrides.is_empty() {
        return Vec::new();
    }
    overrides.par_sort_unstable_by(|a, b| a.0.cmp(&b.0).then(b.1.cmp(&a.1)));
    overrides.dedup_by_key(|o| o.0);
    overrides
        .par_iter()
        .filter_map(|&(pair, _, ord)| pairs.binary_search(&pair).ok().map(|i| (i as u32, ord)))
        .collect()
}

fn verbatim(bonds: &[[u32; 2]], atom_count: usize) -> BondTable {
    let mut pairs: Vec<[u32; 2]> = bonds
        .iter()
        .filter(|&&[a, b]| a != b && (a as usize) < atom_count && (b as usize) < atom_count)
        .map(|&[a, b]| order(a, b))
        .collect();
    pairs.par_sort_unstable();
    pairs.dedup();
    BondTable {
        pairs,
        orders: Vec::new(),
    }
}

#[inline]
fn order(a: u32, b: u32) -> [u32; 2] {
    if a < b {
        [a, b]
    } else {
        [b, a]
    }
}

#[inline]
fn alt_compatible(alt: &[u8], a: u32, b: u32) -> bool {
    let aa = alt.get(a as usize).copied().unwrap_or(0);
    let ab = alt.get(b as usize).copied().unwrap_or(0);
    aa == 0 || ab == 0 || aa == ab
}

/// Below this many candidates, a parallel sort or scan's fork/join
/// overhead costs more than it saves -- real structures' named-link and
/// generic-path candidate lists are almost always this small.
const PARALLEL_THRESHOLD: usize = 50_000;

/// Sorts candidates nearest-first and greedily accepts each one whose
/// atoms are both still under their cap, incrementing running degree
/// counts as it goes. Used both for the valence-capped generic path and
/// for the "at most one partner" named links (disulfide, glycosidic).
fn greedy_by_distance(
    candidates: Vec<(f32, u32, u32)>,
    cap: impl Fn(u32) -> u32 + Sync,
    degree: &mut [u32],
) -> Vec<[u32; 2]> {
    if candidates.is_empty() {
        return Vec::new();
    }
    // An atom whose candidate count never exceeds its cap can never lose
    // a "priority to shortest" contest, so it never needs sorting -- only
    // an over-subscribed atom (more candidates than its cap allows) does.
    // Real structures rarely have any such atom; a packed synthetic one
    // can have many, so this degrades gracefully to a full sort. Touch
    // counting and the accept/contend split are both a pure function of
    // the whole candidate set, not of processing order, so both
    // parallelize; only the shortest-first fill is inherently sequential.
    let parallel = candidates.len() > PARALLEL_THRESHOLD;
    let cap_of = precompute_cap(degree.len(), &cap, parallel);
    let touches = count_touches(&candidates, degree.len(), parallel);
    let under_cap = |a: u32| touches[a as usize] <= cap_of[a as usize];
    let (mut out, contended) = split_by_cap(candidates, under_cap, parallel);

    for &[a, b] in &out {
        degree[a as usize] += 1;
        degree[b as usize] += 1;
    }
    out.extend(accept_nearest_first(contended, &cap_of, degree));
    out
}

/// `cap(a)` for every atom, computed once: both `greedy_by_distance` and
/// `accept_nearest_first` otherwise re-call it (an indirect call plus,
/// for the valence cap, an `element` lookup) for the same atom on every
/// candidate that names it.
fn precompute_cap(
    atom_count: usize,
    cap: &(impl Fn(u32) -> u32 + Sync),
    parallel: bool,
) -> Vec<u32> {
    if parallel {
        (0..atom_count as u32).into_par_iter().map(cap).collect()
    } else {
        (0..atom_count as u32).map(cap).collect()
    }
}

/// One touch count per atom across all candidates, for the cap check
/// above: an atom whose touch count already exceeds its cap has more
/// candidates than it can keep, and needs the shortest-first sort.
fn count_touches(candidates: &[(f32, u32, u32)], atom_count: usize, parallel: bool) -> Vec<u32> {
    if !parallel {
        let mut touches = vec![0u32; atom_count];
        for &(_, a, b) in candidates {
            touches[a as usize] += 1;
            touches[b as usize] += 1;
        }
        return touches;
    }
    let touches: Vec<AtomicU32> = (0..atom_count).map(|_| AtomicU32::new(0)).collect();
    candidates.par_iter().for_each(|&(_, a, b)| {
        touches[a as usize].fetch_add(1, Ordering::Relaxed);
        touches[b as usize].fetch_add(1, Ordering::Relaxed);
    });
    touches.into_iter().map(AtomicU32::into_inner).collect()
}

/// Splits into (both atoms under cap -> accepted bond, otherwise ->
/// still contended), without touching `degree`: an accepted candidate's
/// degree contribution is order-independent (each adds exactly 1 to two
/// atoms), so the caller applies it afterward in one pass over `out`.
fn split_by_cap(
    candidates: Vec<(f32, u32, u32)>,
    under_cap: impl Fn(u32) -> bool + Sync,
    parallel: bool,
) -> (Vec<[u32; 2]>, Vec<(f32, u32, u32)>) {
    let split_chunk = |chunk: &[(f32, u32, u32)]| {
        // Each side capped at `chunk.len()`, its own worst case: never
        // reallocates, unlike growing from empty on every rayon fold leaf.
        let mut out = Vec::with_capacity(chunk.len());
        let mut contended = Vec::with_capacity(chunk.len());
        for &c @ (_, a, b) in chunk {
            if under_cap(a) && under_cap(b) {
                out.push(order(a, b));
            } else {
                contended.push(c);
            }
        }
        (out, contended)
    };
    if !parallel {
        return split_chunk(&candidates);
    }
    let chunk_size = (candidates.len() / (rayon::current_num_threads() * 4)).max(1024);
    candidates.par_chunks(chunk_size).map(split_chunk).reduce(
        || (Vec::new(), Vec::new()),
        |mut a, mut b| {
            a.0.append(&mut b.0);
            a.1.append(&mut b.1);
            a
        },
    )
}

/// Every candidate's sort key is a non-negative, non-NaN distance (or
/// squared distance), so its bit pattern orders the same way as the
/// float itself (IEEE 754) -- sorting by that `u32` skips `total_cmp`'s
/// NaN/sign handling on every comparison. The parallel sort's fork/join
/// overhead only pays off once there is real work to split.
///
/// `degree` already reflects `greedy_by_distance`'s uncontended accepts;
/// nothing downstream reads it afterward, so this greedy fill tracks
/// slots left (`cap_of[a] - degree[a]`) in its own scratch array instead
/// of re-deriving "under cap" from two arrays on every candidate.
fn accept_nearest_first(
    mut candidates: Vec<(f32, u32, u32)>,
    cap_of: &[u32],
    degree: &[u32],
) -> Vec<[u32; 2]> {
    if candidates.len() > PARALLEL_THRESHOLD {
        candidates.par_sort_unstable_by_key(|c| c.0.to_bits());
    } else {
        candidates.sort_unstable_by_key(|c| c.0.to_bits());
    }
    let mut remaining: Vec<u32> = cap_of
        .iter()
        .zip(degree)
        .map(|(&cap, &used)| cap.saturating_sub(used))
        .collect();
    let mut out = Vec::with_capacity(candidates.len());
    for (_, a, b) in candidates {
        let (ai, bi) = (a as usize, b as usize);
        if remaining[ai] > 0 && remaining[bi] > 0 {
            remaining[ai] -= 1;
            remaining[bi] -= 1;
            out.push(order(a, b));
        }
    }
    out
}

struct TemplateResult {
    pairs: Vec<[u32; 2]>,
    /// Non-`Single` order per pair the template named (see
    /// `BondOrder::Aromatic` doc: a template stores the CCD's own Kekulé
    /// `Single`/`Double`, so this is the built-in-table priority source
    /// in `perceive`'s override merge).
    orders: Vec<([u32; 2], BondOrder)>,
    /// Atom is part of a residue whose template coverage passed: excluded
    /// from the generic distance-based path.
    is_templated: Vec<bool>,
    /// Glycan anomeric carbons (C1, or C2 for the sialic acids).
    glycan_anomeric: Vec<u32>,
    /// Candidate glycosidic acceptors: a glycan's other oxygens, and an
    /// ASN's ND2 for N-glycosylation.
    glycosidic_acceptors: Vec<u32>,
}

/// Phase 1: for each residue with a [`templates`] entry, bond its atoms
/// exactly per that template's atom-name pairs (trusting chemistry over
/// geometry, modulo the [`TEMPLATE_HEAVY_MAX`] sanity cap), alt-loc aware.
/// A hydrogen the template doesn't name (nonstandard naming) bonds to its
/// nearest same-residue heavy atom within normal covalent range instead.
/// Per-residue results, accumulated per rayon fold chunk and merged.
#[derive(Default)]
struct TemplateAcc {
    pairs: Vec<[u32; 2]>,
    orders: Vec<([u32; 2], BondOrder)>,
    /// Atom ranges whose residue's template coverage passed (see
    /// [`TemplateResult::is_templated`]); merged into that bitmap
    /// sequentially at the end, since setting bits is cheap and a
    /// `Vec<bool>` can't be split and merged across fold chunks.
    templated: Vec<Range<u32>>,
    glycan_anomeric: Vec<u32>,
    glycosidic_acceptors: Vec<u32>,
}

/// Scratch buffers reused across residues within one fold chunk so
/// matching a residue never allocates: a `HashMap` (or a fresh `Vec`)
/// per residue dominated this function's cost at scale (~350ms of
/// 577ms on a 1M-atom system, before this and the coverage pre-check
/// below).
#[derive(Default)]
struct TemplateScratch<'a> {
    names: Vec<(&'a str, u32)>,
    candidates: Vec<[u32; 2]>,
    order_candidates: Vec<([u32; 2], BondOrder)>,
    heavy_matched: Vec<bool>,
}

fn match_templates(topology: &Topology, positions: &[Vec3]) -> TemplateResult {
    let element = &topology.element;
    let alt = &topology.alt_loc;
    let n = topology.atom_count();

    let acc = topology
        .residues
        .par_iter()
        .fold(
            || (TemplateAcc::default(), TemplateScratch::default()),
            |(mut acc, mut scratch), res| {
                match_one_residue(
                    topology,
                    positions,
                    element,
                    alt,
                    res,
                    &mut acc,
                    &mut scratch,
                );
                (acc, scratch)
            },
        )
        .map(|(acc, _)| acc)
        .reduce(TemplateAcc::default, |mut a, mut b| {
            a.pairs.append(&mut b.pairs);
            a.orders.append(&mut b.orders);
            a.templated.append(&mut b.templated);
            a.glycan_anomeric.append(&mut b.glycan_anomeric);
            a.glycosidic_acceptors.append(&mut b.glycosidic_acceptors);
            a
        });

    let mut is_templated = vec![false; n];
    for range in acc.templated {
        is_templated[range.start as usize..range.end as usize].fill(true);
    }

    TemplateResult {
        pairs: acc.pairs,
        orders: acc.orders,
        is_templated,
        glycan_anomeric: acc.glycan_anomeric,
        glycosidic_acceptors: acc.glycosidic_acceptors,
    }
}

fn match_one_residue<'a>(
    topology: &'a Topology,
    positions: &[Vec3],
    element: &[Element],
    alt: &[u8],
    res: &ResidueRec,
    acc: &mut TemplateAcc,
    scratch: &mut TemplateScratch<'a>,
) {
    let TemplateScratch {
        names,
        candidates,
        order_candidates,
        heavy_matched,
    } = scratch;
    let resname = topology.names.get(res.comp);
    let Some(template) = templates::get(resname) else {
        return;
    };
    let base = res.atoms.start;
    names.clear();
    names.extend(
        res.atoms
            .clone()
            .map(|a| (topology.atom_name(a as usize), a)),
    );
    candidates.clear();
    order_candidates.clear();
    heavy_matched.clear();
    heavy_matched.resize(names.len(), false);

    // Cheap upper-bound reject before the O(pairs x atoms) pairing pass
    // below: a name absent from the template can never end up bonded, so
    // if even *being named right* already falls short of coverage, there
    // is no point computing distances at all. This is what keeps a
    // 1M-atom system where most residues carry generic, non-matching
    // atom names fast.
    let total_heavy = res
        .atoms
        .clone()
        .filter(|&a| !element[a as usize].is_hydrogen())
        .count();
    let present_heavy = names
        .iter()
        .filter(|(nm, a)| {
            !element[*a as usize].is_hydrogen() && template.heavy_names.binary_search(nm).is_ok()
        })
        .count();
    if total_heavy > 1 && (present_heavy as f32) < TEMPLATE_HEAVY_COVERAGE * total_heavy as f32 {
        return;
    }

    for &(n1, n2, bond_order) in template.bonds {
        for &(name_a, a) in names.iter() {
            if name_a != n1 {
                continue;
            }
            for &(name_b, b) in names.iter() {
                if name_b != n2 || a == b || !alt_compatible(alt, a, b) {
                    continue;
                }
                let d = positions[a as usize].distance(positions[b as usize]);
                let has_h = element[a as usize].is_hydrogen() || element[b as usize].is_hydrogen();
                let has_metal = element[a as usize].is_metal() || element[b as usize].is_metal();
                let cap = if has_metal {
                    TEMPLATE_METAL_MAX
                } else if has_h {
                    TEMPLATE_H_MAX
                } else {
                    TEMPLATE_HEAVY_MAX
                };
                if d > cap {
                    continue;
                }
                candidates.push(order(a, b));
                if bond_order != BondOrder::Single {
                    order_candidates.push((order(a, b), bond_order));
                }
                for x in [a, b] {
                    if !element[x as usize].is_hydrogen() {
                        heavy_matched[(x - base) as usize] = true;
                    }
                }
            }
        }
    }

    let matched_heavy = heavy_matched.iter().filter(|&&m| m).count();
    let coverage_ok =
        total_heavy <= 1 || matched_heavy as f32 >= TEMPLATE_HEAVY_COVERAGE * total_heavy as f32;
    if !coverage_ok {
        return; // whole residue falls through to distance-based perception
    }

    acc.templated.push(res.atoms.clone());
    for a in res.atoms.clone() {
        if element[a as usize].is_hydrogen() && !candidates.iter().any(|&[x, y]| x == a || y == a) {
            if let Some(h) = nearest_heavy_in_range(a, res.atoms.clone(), element, positions) {
                candidates.push(order(a, h));
            }
        }
    }
    acc.pairs.extend_from_slice(candidates);
    acc.orders.extend_from_slice(order_candidates);

    if templates::GLYCANS.contains(&resname) {
        let anomeric_name = if templates::SIALIC_ACIDS.contains(&resname) {
            "C2"
        } else {
            "C1"
        };
        acc.glycan_anomeric.extend(
            names
                .iter()
                .filter(|(nm, _)| *nm == anomeric_name)
                .map(|&(_, a)| a),
        );
        acc.glycosidic_acceptors.extend(
            res.atoms
                .clone()
                .filter(|&a| element[a as usize] == Element::OXYGEN),
        );
    } else if resname == "ASN" {
        acc.glycosidic_acceptors
            .extend(names.iter().filter(|(nm, _)| *nm == "ND2").map(|&(_, a)| a));
    }
}

fn nearest_heavy_in_range(
    h: u32,
    range: Range<u32>,
    element: &[Element],
    positions: &[Vec3],
) -> Option<u32> {
    let hp = positions[h as usize];
    let he = element[h as usize];
    range
        .filter(|&a| a != h && !element[a as usize].is_hydrogen())
        .map(|a| (a, hp.distance(positions[a as usize])))
        .filter(|&(a, d)| {
            d <= he.covalent_radius() + element[a as usize].covalent_radius() + COVALENT_TOLERANCE
        })
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .map(|(a, _)| a)
}

/// Phase 2a: peptide C(i)-N(i+1) and nucleic O3'(i)-P(i+1), for every
/// pair of residues adjacent in `topology.residues` within the same
/// chain (works for a templated or an untemplated/modified residue
/// alike -- MSE, PTR, SEP and friends still need their peptide bond). A
/// chain break simply fails the distance window.
fn link_adjacent_residues(topology: &Topology, positions: &[Vec3]) -> Vec<[u32; 2]> {
    // A structure with no phosphorus at all has no nucleic backbone, so
    // skip that scan entirely: on an all-protein system this halves the
    // per-residue-pair work (each `named_link` miss otherwise scans the
    // whole residue for a name that is never there).
    let has_nucleic = topology.element.iter().any(|e| e.atomic_number() == 15);
    // Chains are independent (a link never crosses one), and each pair's
    // name lookup is a few atoms' worth of work -- fine-grained enough
    // that a big system (many chains) benefits from spreading it out.
    topology
        .chains
        .par_iter()
        .fold(Vec::new, |mut out, chain| {
            for r in chain.residues.start..chain.residues.end.saturating_sub(1) {
                let (res, next) = (
                    &topology.residues[r as usize],
                    &topology.residues[r as usize + 1],
                );
                if let Some(bond) = named_link(
                    topology,
                    positions,
                    res.atoms.clone(),
                    next.atoms.clone(),
                    "C",
                    "N",
                    PEPTIDE_BOND_MAX,
                ) {
                    out.push(bond);
                }
                if has_nucleic {
                    if let Some(bond) = named_link(
                        topology,
                        positions,
                        res.atoms.clone(),
                        next.atoms.clone(),
                        "O3'",
                        "P",
                        NUCLEIC_LINK_MAX,
                    ) {
                        out.push(bond);
                    }
                }
            }
            out
        })
        .reduce(Vec::new, |mut a, mut b| {
            a.append(&mut b);
            a
        })
}

fn named_link(
    topology: &Topology,
    positions: &[Vec3],
    mut from: Range<u32>,
    mut to: Range<u32>,
    from_name: &str,
    to_name: &str,
    max_dist: f32,
) -> Option<[u32; 2]> {
    let a = from.find(|&a| topology.atom_name(a as usize) == from_name)?;
    let b = to.find(|&b| topology.atom_name(b as usize) == to_name)?;
    let d = positions[a as usize].distance(positions[b as usize]);
    (LINK_MIN..=max_dist).contains(&d).then(|| order(a, b))
}

/// Phase 2b: disulfide SG-SG, at most one partner each, nearest first.
/// Not restricted to adjacent residues or a particular residue name: any
/// two sulfurs named `SG` within [`DISULFIDE_MAX`] are a disulfide.
fn link_disulfides(topology: &Topology, positions: &[Vec3]) -> Vec<[u32; 2]> {
    let sg: Vec<u32> = (0..topology.atom_count() as u32)
        .filter(|&a| {
            topology.element[a as usize] == Element::SULFUR
                && topology.atom_name(a as usize) == "SG"
        })
        .collect();
    let mut candidates = Vec::new();
    for (i, &a) in sg.iter().enumerate() {
        for &b in &sg[i + 1..] {
            if topology.residue_index[a as usize] == topology.residue_index[b as usize] {
                continue;
            }
            let d = positions[a as usize].distance(positions[b as usize]);
            if (LINK_MIN..=DISULFIDE_MAX).contains(&d) {
                candidates.push((d, a, b));
            }
        }
    }
    let mut degree = vec![0u32; topology.atom_count()];
    greedy_by_distance(candidates, |_| 1, &mut degree)
}

/// Phase 2c: glycosidic links (glycan anomeric carbon to an acceptor
/// oxygen or an N-glycosylation ASN's ND2), at most one partner each,
/// nearest first.
fn link_glycosidic(
    templated: &TemplateResult,
    positions: &[Vec3],
    topology: &Topology,
) -> Vec<[u32; 2]> {
    let mut candidates = Vec::new();
    for &a in &templated.glycan_anomeric {
        for &b in &templated.glycosidic_acceptors {
            if topology.residue_index[a as usize] == topology.residue_index[b as usize] {
                continue;
            }
            let d = positions[a as usize].distance(positions[b as usize]);
            if (LINK_MIN..=GLYCOSIDIC_MAX).contains(&d) {
                candidates.push((d, a, b));
            }
        }
    }
    let mut degree = vec![0u32; topology.atom_count()];
    greedy_by_distance(candidates, |_| 1, &mut degree)
}

/// Phase 3: distance-based perception (the original algorithm) restricted
/// to atoms no template covers, excluding metals and single-atom (ion)
/// residues entirely -- those bond only via a cofactor template or an
/// explicit file link, never by proximity (task #4). A clash guard drops
/// pairs much closer than any real bond of that element pair, and a
/// per-element valence cap keeps the nearest candidates and drops the
/// rest once an atom is full.
fn generic_distance_bonds(
    positions: &[Vec3],
    element: &[Element],
    alt: &[u8],
    is_templated: &[bool],
    topology: &Topology,
) -> Vec<[u32; 2]> {
    let n = positions.len();
    // Atom indices only, not a `Vec<bool>` the size of the whole system:
    // the grid below is then built from just these, so a near-fully
    // templated real structure (almost nothing eligible) gets a small,
    // cache-tight grid instead of one sized to every atom, and every
    // neighbor `j` the inner loop visits is eligible by construction --
    // no per-candidate lookup to re-check it.
    let eligible: Vec<u32> = (0..n as u32)
        .filter(|&a| {
            let a = a as usize;
            !is_templated[a]
                && !element[a].is_metal()
                && topology.residues[topology.residue_index[a] as usize]
                    .atoms
                    .len()
                    > 1
        })
        .collect();
    if eligible.is_empty() {
        return Vec::new();
    }

    let grid = Grid::new(&eligible, positions);
    let candidates: Vec<(f32, u32, u32)> = eligible
        .par_iter()
        .fold(Vec::new, |mut out, &i| {
            let i = i as usize;
            let p = positions[i];
            let e = element[i];
            let r = e.covalent_radius();
            let (cx, cy, cz) = grid.cell_of(p);
            for dz in -1..=1 {
                for dy in -1..=1 {
                    for dx in -1..=1 {
                        let Some(cell) = grid.cell_index(cx + dx, cy + dy, cz + dz) else {
                            continue;
                        };
                        for &j in grid.atoms_in(cell) {
                            let j = j as usize;
                            if j <= i || !alt_compatible(alt, i as u32, j as u32) {
                                continue;
                            }
                            let rj = element[j].covalent_radius();
                            let ideal = r + rj;
                            // Squared distance: avoids a sqrt on every
                            // checked neighbor, not just accepted ones
                            // (this loop runs tens of millions of times
                            // on a 1M-atom system). Ordering by d2 sorts
                            // identically to ordering by d.
                            let d2 = p.distance_squared(positions[j]);
                            let lo = CLASH_FACTOR * ideal;
                            let hi = ideal + COVALENT_TOLERANCE;
                            if d2 >= lo * lo && d2 <= hi * hi {
                                out.push((d2, i as u32, j as u32));
                            }
                        }
                    }
                }
            }
            out
        })
        .reduce(Vec::new, |mut a, mut b| {
            a.append(&mut b);
            a
        });

    let mut degree = vec![0u32; n];
    greedy_by_distance(
        candidates,
        |a| valence_cap(element[a as usize]),
        &mut degree,
    )
}

/// Uniform grid over the bounding box with `MAX_BOND` cells, built with a
/// counting sort so every cell is a contiguous slice of atom indices.
/// Covers only the given atoms (their original indices, kept as-is in
/// `atoms` -- not their position within the input slice).
struct Grid {
    origin: Vec3,
    dims: [i32; 3],
    starts: Vec<u32>,
    atoms: Vec<u32>,
}

impl Grid {
    /// `atoms` must be non-empty; the caller (`generic_distance_bonds`)
    /// returns early instead of building an empty grid.
    fn new(atoms: &[u32], positions: &[Vec3]) -> Self {
        let first = positions[atoms[0] as usize];
        let (min, max) = atoms.iter().fold((first, first), |(lo, hi), &a| {
            let p = positions[a as usize];
            (lo.min(p), hi.max(p))
        });
        let extent = max - min;
        let dims = [
            (extent.x / MAX_BOND).floor() as i32 + 1,
            (extent.y / MAX_BOND).floor() as i32 + 1,
            (extent.z / MAX_BOND).floor() as i32 + 1,
        ];
        let cell_count = dims[0] as usize * dims[1] as usize * dims[2] as usize;
        let mut grid = Grid {
            origin: min,
            dims,
            starts: vec![0; cell_count + 1],
            atoms: Vec::new(),
        };

        let cells: Vec<u32> = atoms
            .par_iter()
            .map(|&a| {
                let (x, y, z) = grid.cell_of(positions[a as usize]);
                grid.cell_index(x, y, z).expect("in bounds") as u32
            })
            .collect();
        for &c in &cells {
            grid.starts[c as usize + 1] += 1;
        }
        for c in 0..cell_count {
            grid.starts[c + 1] += grid.starts[c];
        }
        let mut cursor = grid.starts.clone();
        grid.atoms = vec![0; atoms.len()];
        for (&a, &c) in atoms.iter().zip(cells.iter()) {
            let slot = &mut cursor[c as usize];
            grid.atoms[*slot as usize] = a;
            *slot += 1;
        }
        grid
    }

    fn cell_of(&self, p: Vec3) -> (i32, i32, i32) {
        let c = (p - self.origin) / MAX_BOND;
        (c.x.floor() as i32, c.y.floor() as i32, c.z.floor() as i32)
    }

    fn cell_index(&self, x: i32, y: i32, z: i32) -> Option<usize> {
        if x < 0 || y < 0 || z < 0 || x >= self.dims[0] || y >= self.dims[1] || z >= self.dims[2] {
            return None;
        }
        Some((z as usize * self.dims[1] as usize + y as usize) * self.dims[0] as usize + x as usize)
    }

    fn atoms_in(&self, cell: usize) -> &[u32] {
        &self.atoms[self.starts[cell] as usize..self.starts[cell + 1] as usize]
    }
}

impl ExplicitBond {
    pub fn new(a: u32, b: u32, kind: ExplicitBondKind) -> Self {
        Self {
            atoms: [a, b],
            kind,
            order: BondOrder::Single,
        }
    }

    pub fn with_order(a: u32, b: u32, kind: ExplicitBondKind, order: BondOrder) -> Self {
        Self {
            atoms: [a, b],
            kind,
            order,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ChainRec, Interner, ResidueRec, SecondaryStructure};

    fn topology(elements: &[Element]) -> Topology {
        Topology {
            element: elements.to_vec(),
            alt_loc: vec![0; elements.len()],
            names: Interner::new(),
            ..Default::default()
        }
    }

    #[test]
    fn bonds_by_covalent_distance() {
        // Water-like: O with two H at 0.96 A, plus a far-away C. None of
        // these atoms have a residue name, so this is the untemplated path.
        let mut t = topology(&[
            Element::OXYGEN,
            Element::HYDROGEN,
            Element::HYDROGEN,
            Element::CARBON,
        ]);
        one_residue(&mut t, "LIG");
        let p = [
            Vec3::ZERO,
            Vec3::new(0.96, 0.0, 0.0),
            Vec3::new(-0.24, 0.93, 0.0),
            Vec3::new(5.0, 0.0, 0.0),
        ];
        let bonds = perceive(&t, &p);
        assert_eq!(bonds.pairs, vec![[0, 1], [0, 2]]);
        assert!(bonds.contains(1, 0));
        assert!(!bonds.contains(0, 3));
        assert_eq!(bonds.degrees(4), vec![2, 1, 1, 0]);
        // The water (atoms 0,1,2) is one fragment; the far-away, unbonded
        // carbon (atom 3) is its own.
        let fragments = bonds.fragments(4);
        assert_eq!(fragments[0], fragments[1]);
        assert_eq!(fragments[0], fragments[2]);
        assert_ne!(fragments[0], fragments[3]);
    }

    #[test]
    fn fragments_are_numbered_in_atom_order_regardless_of_bond_order() {
        // Atom 3 bonds to atom 0 first, so its fragment must still get id 0
        // (the lowest atom index it touches), not id 1 from being visited
        // via a later pair.
        let bonds = BondTable {
            pairs: vec![[0, 3], [1, 2]],
            orders: Vec::new(),
        };
        assert_eq!(bonds.fragments(4), vec![0, 1, 1, 0]);
    }

    #[test]
    fn hydrogen_keeps_only_its_shortest_bond() {
        let mut t = topology(&[Element::CARBON, Element::HYDROGEN, Element::OXYGEN]);
        one_residue(&mut t, "LIG");
        // H is 1.09 from C and 1.3 from O: both within cutoff, keep C-H only.
        let p = [
            Vec3::ZERO,
            Vec3::new(1.09, 0.0, 0.0),
            Vec3::new(2.39, 0.0, 0.0),
        ];
        let bonds = perceive(&t, &p);
        assert_eq!(bonds.pairs, vec![[0, 1]]);
    }

    #[test]
    fn alternate_locations_do_not_bond_across_each_other() {
        let mut t = topology(&[Element::CARBON, Element::CARBON, Element::CARBON]);
        one_residue(&mut t, "LIG");
        t.alt_loc = vec![b'A', b'B', 0];
        // Atom 2 is 1.5 A from both; atoms 0 and 1 are 1.5 A apart but in
        // different alternate locations.
        let p = [
            Vec3::ZERO,
            Vec3::new(1.5, 0.0, 0.0),
            Vec3::new(0.75, 1.3, 0.0),
        ];
        let bonds = perceive(&t, &p);
        assert_eq!(bonds.pairs, vec![[0, 2], [1, 2]]);
    }

    #[test]
    fn explicit_bonds_are_merged_and_deduplicated() {
        let mut t = topology(&[Element::SULFUR, Element::SULFUR, Element::CARBON]);
        one_residue(&mut t, "LIG");
        t.explicit_bonds
            .push(ExplicitBond::new(2, 0, ExplicitBondKind::Covalent));
        // S-S at 2.03 A is perceived; C at 4 A is only bonded explicitly.
        let p = [
            Vec3::ZERO,
            Vec3::new(2.03, 0.0, 0.0),
            Vec3::new(0.0, 4.0, 0.0),
        ];
        let bonds = perceive(&t, &p);
        assert_eq!(bonds.pairs, vec![[0, 1], [0, 2]]);
    }

    #[test]
    fn md_topology_bonds_are_used_verbatim_and_skip_perception() {
        // Atoms placed 10 A apart: geometry alone would perceive nothing.
        let mut t = topology(&[Element::CARBON, Element::CARBON, Element::CARBON]);
        one_residue(&mut t, "LIG");
        t.md_bonds = Some(vec![[0, 1], [1, 2]]);
        let p = [
            Vec3::ZERO,
            Vec3::new(10.0, 0.0, 0.0),
            Vec3::new(20.0, 0.0, 0.0),
        ];
        assert_eq!(perceive(&t, &p).pairs, vec![[0, 1], [1, 2]]);
    }

    #[test]
    fn md_bonds_are_normalized_and_deduplicated() {
        let mut t = topology(&[Element::CARBON, Element::CARBON, Element::CARBON]);
        one_residue(&mut t, "LIG");
        // Reversed, duplicated, self and out-of-range pairs.
        t.md_bonds = Some(vec![[1, 0], [0, 1], [2, 2], [0, 9]]);
        let p = [
            Vec3::ZERO,
            Vec3::new(20.0, 0.0, 0.0),
            Vec3::new(40.0, 0.0, 0.0),
        ];
        assert_eq!(perceive(&t, &p).pairs, vec![[0, 1]]);
    }

    #[test]
    fn metals_get_no_geometric_bonds() {
        let mut t = topology(&[Element::from_atomic_number(26).unwrap(), Element::NITROGEN]);
        one_residue_each(&mut t, &["FE", "LIG"]);
        let p = [Vec3::ZERO, Vec3::new(2.1, 0.0, 0.0)];
        assert!(
            perceive(&t, &p).is_empty(),
            "a bare metal ion gets no default geometric bond, even close"
        );
    }

    #[test]
    fn metal_near_water_does_not_bond() {
        let mut t = topology(&[Element::from_atomic_number(30).unwrap(), Element::OXYGEN]);
        one_residue_each(&mut t, &["ZN", "HOH"]);
        let p = [Vec3::ZERO, Vec3::new(2.2, 0.0, 0.0)];
        assert!(perceive(&t, &p).is_empty());
    }

    #[test]
    fn clash_guard_rejects_atoms_too_close_for_a_real_bond() {
        let mut t = topology(&[Element::CARBON, Element::CARBON]);
        one_residue(&mut t, "LIG");
        // Ideal C-C is 1.52 A; 0.5 A is a clash (two different residues'
        // atoms placed too close), not a bond.
        let p = [Vec3::ZERO, Vec3::new(0.5, 0.0, 0.0)];
        assert!(perceive(&t, &p).is_empty());
    }

    #[test]
    fn valence_cap_keeps_only_the_nearest_bonds() {
        // A lone carbon surrounded by 5 others (triangular-bipyramid
        // directions, so the 5 satellites are always >2.6 A from each
        // other -- never candidates for a bond among themselves) each
        // within carbon's bonding distance but at increasing radius: only
        // 4 (the cap) survive, and the farthest (largest radius) is the
        // one dropped.
        let elements = vec![Element::CARBON; 6];
        let mut t = topology(&elements);
        one_residue(&mut t, "LIG");
        let dirs = [
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(-0.5, 0.866, 0.0),
            Vec3::new(-0.5, -0.866, 0.0),
        ];
        let mut p = vec![Vec3::ZERO];
        p.extend((0..5).map(|i| dirs[i] * (1.85 + i as f32 * 0.02)));
        let bonds = perceive(&t, &p);
        assert_eq!(bonds.degrees(6)[0], 4, "carbon caps at 4 bonds");
        assert!(
            !bonds.contains(0, 5),
            "the farthest (5th) satellite is dropped"
        );
        for satellite in 1..5 {
            assert!(
                bonds.contains(0, satellite),
                "nearer satellite {satellite} kept"
            );
        }
    }

    #[test]
    fn standard_residue_template_ignores_a_same_residue_clash() {
        // ALA: N, CA, C, O, CB at real geometry, but CB placed 1.7 A from
        // O -- an unrelated same-residue clash a purely geometric pass
        // would wrongly bond (O's covalent radius + C's + tolerance is
        // well past 1.7 A). The template must not draw that extra bond.
        let mut t = topology(&[
            Element::NITROGEN,
            Element::CARBON,
            Element::CARBON,
            Element::OXYGEN,
            Element::CARBON,
        ]);
        one_residue(&mut t, "ALA");
        for (a, name) in ["N", "CA", "C", "O", "CB"].iter().enumerate() {
            set_name(&mut t, a, name);
        }
        let p = [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.46, 0.0, 0.0),
            Vec3::new(1.46 + 1.52, 0.0, 0.0),
            Vec3::new(1.46 + 1.52, 1.23, 0.0),
            Vec3::new(1.46, 1.53, 0.0),
        ];
        let bonds = perceive(&t, &p);
        let ca = 1;
        let cb = 4;
        let o = 3;
        assert!(bonds.contains(ca, cb), "template CA-CB");
        assert!(!bonds.contains(o, cb), "no clash bond outside the template");
    }

    /// The headline complaint this module exists to fix: two carbons of
    /// *different* residues 1.7 A apart (a real crystal-packing or
    /// compressed-MD-frame clash) must not become a bond. Two full ALA
    /// residues, each internally correct (so both pass template
    /// coverage), in separate chains (so the peptide-link rule never
    /// applies), with only their CB atoms close: 1.7 A clears both the
    /// clash floor (0.91 A) and the old flat cutoff (1.97 A), so only the
    /// "a templated residue bonds another residue solely through a named
    /// link" rule -- not distance -- keeps them apart. Renaming both
    /// residues to an unrecognized name (no template) proves it really is
    /// that rule: the identical geometry then does bond.
    #[test]
    fn two_different_residues_1_7_angstroms_apart_do_not_clash_bond() {
        let (t, p) = two_ala_residues_with_close_cb();
        let bonds = perceive(&t, &p);
        let (cb0, cb1) = (4u32, 9u32);
        assert!((p[cb0 as usize].distance(p[cb1 as usize]) - 1.7).abs() < 1e-4);
        assert!(
            !bonds.contains(cb0, cb1),
            "standard residues never clash-bond by distance"
        );

        let mut lig = t;
        let lig_name = lig.names.intern("LIG");
        for res in &mut lig.residues {
            res.comp = lig_name;
        }
        let bonds = perceive(&lig, &p);
        assert!(
            bonds.contains(cb0, cb1),
            "the same geometry bonds once neither residue has a template"
        );
    }

    fn two_ala_residues_with_close_cb() -> (Topology, Vec<Vec3>) {
        let one = [
            Element::NITROGEN,
            Element::CARBON,
            Element::CARBON,
            Element::OXYGEN,
            Element::CARBON,
        ];
        let mut t = topology(&one.repeat(2));
        let a = t.names.intern("A");
        let b = t.names.intern("B");
        let ala = t.names.intern("ALA");
        t.name = vec![*b"X   "; 10];
        for atoms in [0..5u32, 5..10u32] {
            for (i, atom) in atoms.enumerate() {
                set_name(&mut t, atom as usize, ["N", "CA", "C", "O", "CB"][i]);
            }
        }
        t.residue_index = vec![0, 0, 0, 0, 0, 1, 1, 1, 1, 1];
        t.residues = vec![
            ResidueRec {
                atoms: 0..5,
                chain: 0,
                comp: ala,
                seq_id: 1,
                auth_seq_id: 1,
                ins_code: 0,
                ss: SecondaryStructure::Unknown,
            },
            ResidueRec {
                atoms: 5..10,
                chain: 1,
                comp: ala,
                seq_id: 1,
                auth_seq_id: 1,
                ins_code: 0,
                ss: SecondaryStructure::Unknown,
            },
        ];
        t.chains = vec![
            ChainRec {
                residues: 0..1,
                label_asym: a,
                auth_asym: a,
                entity: 1,
            },
            ChainRec {
                residues: 1..2,
                label_asym: b,
                auth_asym: b,
                entity: 2,
            },
        ];
        // Residue 0 at real ALA geometry; residue 1 is the same rigid
        // shape translated so its CB lands exactly 1.7 A from residue 0's
        // CB (see the test doc), with nothing else of residue 1 anywhere
        // near residue 0 except (deliberately) N0-N1 also under 2 A --
        // a second clash the template must also ignore.
        let shape = [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.46, 0.0, 0.0),
            Vec3::new(1.46 + 1.52, 0.0, 0.0),
            Vec3::new(1.46 + 1.52, 1.23, 0.0),
            Vec3::new(1.46, 1.53, 0.0),
        ];
        let translation = Vec3::new(1.7, 0.0, 0.0);
        let mut p = shape.to_vec();
        p.extend(shape.iter().map(|&v| v + translation));
        (t, p)
    }

    #[test]
    fn peptide_bond_links_adjacent_residues() {
        let mut t = topology(&[Element::CARBON, Element::NITROGEN]);
        two_residues(&mut t, "ALA", "GLY");
        set_name(&mut t, 0, "C");
        set_name(&mut t, 1, "N");
        let p = [Vec3::ZERO, Vec3::new(1.33, 0.0, 0.0)];
        assert_eq!(perceive(&t, &p).pairs, vec![[0, 1]]);
    }

    #[test]
    fn disulfide_bond_between_far_apart_cysteines() {
        let mut t = topology(&[Element::SULFUR, Element::SULFUR]);
        two_residues(&mut t, "CYS", "CYS");
        set_name(&mut t, 0, "SG");
        set_name(&mut t, 1, "SG");
        let p = [Vec3::ZERO, Vec3::new(2.05, 0.0, 0.0)];
        assert_eq!(perceive(&t, &p).pairs, vec![[0, 1]]);
    }

    #[test]
    fn grid_finds_neighbors_across_cell_borders() {
        // A chain of carbons 1.5 A apart spanning many 3 A cells.
        let n = 200;
        let t = {
            let mut t = topology(&vec![Element::CARBON; n]);
            one_residue(&mut t, "LIG");
            t
        };
        let p: Vec<Vec3> = (0..n)
            .map(|i| Vec3::new(i as f32 * 1.5, 0.0, 0.0))
            .collect();
        let bonds = perceive(&t, &p);
        assert_eq!(bonds.len(), n - 1);
        assert!((0..n - 1).all(|i| bonds.contains(i as u32, i as u32 + 1)));
    }

    // -- test topology helpers -------------------------------------------

    /// Puts every atom into one multi-atom residue named `resname` (a
    /// ligand-like untemplated residue, unless `resname` has a template).
    fn one_residue(t: &mut Topology, resname: &str) {
        let n = t.element.len();
        let a = t.names.intern("A");
        let comp = t.names.intern(resname);
        t.name = vec![*b"X   "; n];
        t.residue_index = vec![0; n];
        t.residues = vec![ResidueRec {
            atoms: 0..n as u32,
            chain: 0,
            comp,
            seq_id: 1,
            auth_seq_id: 1,
            ins_code: 0,
            ss: SecondaryStructure::Unknown,
        }];
        t.chains = vec![ChainRec {
            residues: 0..1,
            label_asym: a,
            auth_asym: a,
            entity: 1,
        }];
    }

    /// Puts every atom in its own single-atom residue named `names[i]`,
    /// so a metal ion / ligand / water each get their own residue.
    fn one_residue_each(t: &mut Topology, names: &[&str]) {
        let n = t.element.len();
        let a = t.names.intern("A");
        t.name = vec![*b"X   "; n];
        t.residue_index = (0..n as u32).collect();
        t.residues = names
            .iter()
            .enumerate()
            .map(|(i, name)| {
                let comp = t.names.intern(name);
                ResidueRec {
                    atoms: i as u32..i as u32 + 1,
                    chain: 0,
                    comp,
                    seq_id: i as i32 + 1,
                    auth_seq_id: i as i32 + 1,
                    ins_code: 0,
                    ss: SecondaryStructure::Unknown,
                }
            })
            .collect();
        t.chains = vec![ChainRec {
            residues: 0..names.len() as u32,
            label_asym: a,
            auth_asym: a,
            entity: 1,
        }];
    }

    fn two_residues(t: &mut Topology, first: &str, second: &str) {
        let a = t.names.intern("A");
        let n1 = t.names.intern(first);
        let n2 = t.names.intern(second);
        t.name = vec![*b"X   "; 2];
        t.residue_index = vec![0, 1];
        t.residues = vec![
            ResidueRec {
                atoms: 0..1,
                chain: 0,
                comp: n1,
                seq_id: 1,
                auth_seq_id: 1,
                ins_code: 0,
                ss: SecondaryStructure::Unknown,
            },
            ResidueRec {
                atoms: 1..2,
                chain: 0,
                comp: n2,
                seq_id: 2,
                auth_seq_id: 2,
                ins_code: 0,
                ss: SecondaryStructure::Unknown,
            },
        ];
        t.chains = vec![ChainRec {
            residues: 0..2,
            label_asym: a,
            auth_asym: a,
            entity: 1,
        }];
    }

    fn set_name(t: &mut Topology, atom: usize, name: &str) {
        let mut buf = [b' '; 4];
        buf[..name.len()].copy_from_slice(name.as_bytes());
        t.name[atom] = buf;
    }
}
