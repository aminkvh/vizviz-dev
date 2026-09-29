//! Deterministic synthetic "protein-like" structures for benchmarks.
//!
//! Not chemistry: chains are random walks of pseudo-residues (8 heavy atoms
//! each) packed at real protein atom density inside a sphere. Same
//! `SynthParams` always produces the same structure, on any platform.

use rayon::prelude::*;
use vv_core::glam::Vec3;
use vv_core::{
    ChainRec, CoordSet, Element, Interner, ResidueRec, SecondaryStructure, Structure, Topology,
};

pub const ATOMS_PER_RESIDUE: usize = 8;
/// Heavy-atom density of folded protein, cubic angstroms per atom.
const VOLUME_PER_ATOM: f32 = 20.0;
const CA_SPACING: f32 = 3.8;

#[derive(Clone, Debug)]
pub struct SynthParams {
    pub atom_count: usize,
    pub seed: u64,
    pub residues_per_chain: usize,
}

impl SynthParams {
    pub fn new(atom_count: usize) -> Self {
        Self {
            atom_count,
            seed: 0x5EED,
            residues_per_chain: 300,
        }
    }
}

const RESIDUE_NAMES: [&str; 20] = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
];

const NAMES: [[u8; 4]; ATOMS_PER_RESIDUE] = [
    *b"N   ", *b"CA  ", *b"C   ", *b"O   ", *b"CB  ", *b"CG  ", *b"CD  ", *b"CE  ",
];

/// Element pattern per residue: N, CA, C, O, CB, CG, CD*, CE. CD* varies
/// per residue (O/N/S) so element-based coloring has some variety.
const ELEMENTS: [Element; ATOMS_PER_RESIDUE] = [
    Element::NITROGEN,
    Element::CARBON,
    Element::CARBON,
    Element::OXYGEN,
    Element::CARBON,
    Element::CARBON,
    Element::CARBON,
    Element::CARBON,
];

/// Sphere radius that holds `atom_count` atoms at protein density.
pub fn packing_radius(atom_count: usize) -> f32 {
    let volume = atom_count as f32 * VOLUME_PER_ATOM;
    (volume * 3.0 / (4.0 * std::f32::consts::PI)).cbrt()
}

pub fn protein_like(params: &SynthParams) -> Structure {
    let atom_count = params.atom_count;
    let residues_per_chain = params.residues_per_chain.max(1);
    let residue_count = atom_count.div_ceil(ATOMS_PER_RESIDUE);
    let chain_count = residue_count.div_ceil(residues_per_chain);
    let atoms_per_chain = residues_per_chain * ATOMS_PER_RESIDUE;
    let radius = packing_radius(atom_count);

    let mut positions = vec![Vec3::ZERO; atom_count];
    let mut element = vec![Element::UNKNOWN; atom_count];

    positions
        .par_chunks_mut(atoms_per_chain)
        .zip(element.par_chunks_mut(atoms_per_chain))
        .enumerate()
        .for_each(|(chain, (pos, elem))| {
            let mut rng =
                SplitMix64::new(params.seed ^ (chain as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15));
            fill_chain(&mut rng, radius, pos, elem);
        });

    let mut names = Interner::new();
    let comps: Vec<_> = RESIDUE_NAMES.iter().map(|n| names.intern(n)).collect();

    let mut residue_index = vec![0u32; atom_count];
    let mut residues = Vec::with_capacity(residue_count);
    let mut chains = Vec::with_capacity(chain_count);
    let mut name_buf = String::new();
    for chain in 0..chain_count {
        let first_res = chain * residues_per_chain;
        let last_res = ((chain + 1) * residues_per_chain).min(residue_count);
        for r in first_res..last_res {
            let start = r * ATOMS_PER_RESIDUE;
            let end = ((r + 1) * ATOMS_PER_RESIDUE).min(atom_count);
            residue_index[start..end].fill(r as u32);
            let mut h = SplitMix64::new(params.seed ^ (r as u64 + 1));
            let seq_id = (r - first_res) as i32 + 1;
            residues.push(ResidueRec {
                atoms: start as u32..end as u32,
                chain: chain as u32,
                comp: comps[(h.next_u64() % RESIDUE_NAMES.len() as u64) as usize],
                seq_id,
                auth_seq_id: seq_id,
                ins_code: 0,
                ss: SecondaryStructure::Unknown,
            });
        }
        chain_name(chain, &mut name_buf);
        let id = names.intern(&name_buf);
        chains.push(ChainRec {
            residues: first_res as u32..last_res as u32,
            label_asym: id,
            auth_asym: id,
            entity: 1,
        });
    }

    let topology = Topology {
        element,
        name: (0..atom_count)
            .map(|a| NAMES[a % ATOMS_PER_RESIDUE])
            .collect(),
        serial: (1..=atom_count as u32).collect(),
        b_factor: vec![0.0; atom_count],
        occupancy: vec![1.0; atom_count],
        alt_loc: vec![0; atom_count],
        charge: vec![0; atom_count],
        flags: vec![0; atom_count],
        residue_index,
        residues,
        chains,
        names,
        id: "SYNTH".to_string(),
        ..Default::default()
    };
    debug_assert_eq!(topology.validate(), Ok(()));
    Structure::new(topology, CoordSet::new(positions)).expect("generator produced consistent sizes")
}

/// Random walk of CA positions confined to the packing sphere; each residue
/// hangs its other atoms off a local frame at the CA.
fn fill_chain(rng: &mut SplitMix64, radius: f32, pos: &mut [Vec3], elem: &mut [Element]) {
    let mut ca = rng.in_ball() * radius * 0.9;
    let mut dir = rng.unit_vector();
    for (r, (res_pos, res_elem)) in pos
        .chunks_mut(ATOMS_PER_RESIDUE)
        .zip(elem.chunks_mut(ATOMS_PER_RESIDUE))
        .enumerate()
    {
        if r > 0 {
            dir = (dir + rng.unit_vector() * 0.6).normalize_or(dir);
            let next = ca + dir * CA_SPACING;
            if next.length() > radius {
                dir = (-ca.normalize_or(Vec3::X) + rng.unit_vector() * 0.3).normalize_or(-dir);
            }
            ca += dir * CA_SPACING;
        }
        let t = dir;
        let up = if t.y.abs() < 0.9 { Vec3::Y } else { Vec3::X };
        let n = t.cross(up).normalize();
        let b = t.cross(n);
        let atoms = [
            ca - t * 1.46 + n * 0.4,
            ca,
            ca + t * 1.52,
            ca + t * 1.52 + b * 1.23,
            ca + n * 1.53,
            ca + n * 1.53 + (n * 0.5 + b * 0.87) * 1.52,
            ca + n * 1.53 + (n * 0.5 + b * 0.87) * 1.52 + t * 1.5,
            ca + n * 1.53 + (n * 0.5 + b * 0.87) * 1.52 + t * 1.5 + n * 1.5,
        ];
        let delta = match rng.next_u64() % 4 {
            0 => Element::OXYGEN,
            1 => Element::NITROGEN,
            2 => Element::SULFUR,
            _ => Element::CARBON,
        };
        for (i, (p, e)) in res_pos.iter_mut().zip(res_elem.iter_mut()).enumerate() {
            *p = atoms[i];
            *e = if i == 6 { delta } else { ELEMENTS[i] };
        }
    }
}

/// mmCIF-style chain ids: A..Z, AA..AZ, BA.., then three letters.
fn chain_name(index: usize, out: &mut String) {
    out.clear();
    let mut n = index;
    loop {
        out.insert(0, (b'A' + (n % 26) as u8) as char);
        n /= 26;
        if n == 0 {
            break;
        }
        n -= 1;
    }
}

/// SplitMix64: tiny, fast, good enough for layout noise; deterministic
/// across platforms (pure integer arithmetic).
struct SplitMix64(u64);

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in [0, 1).
    fn next_f32(&mut self) -> f32 {
        (self.next_u64() >> 40) as f32 / (1u64 << 24) as f32
    }

    /// Uniform in [-1, 1).
    fn symmetric(&mut self) -> f32 {
        self.next_f32() * 2.0 - 1.0
    }

    fn unit_vector(&mut self) -> Vec3 {
        loop {
            let v = Vec3::new(self.symmetric(), self.symmetric(), self.symmetric());
            let len_sq = v.length_squared();
            if len_sq > 1e-4 && len_sq <= 1.0 {
                return v / len_sq.sqrt();
            }
        }
    }

    /// Uniform inside the unit ball.
    fn in_ball(&mut self) -> Vec3 {
        loop {
            let v = Vec3::new(self.symmetric(), self.symmetric(), self.symmetric());
            if v.length_squared() <= 1.0 {
                return v;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_atom_count_and_valid_hierarchy() {
        for n in [1usize, 7, 8, 9, 1234, 24_000] {
            let s = protein_like(&SynthParams {
                residues_per_chain: 50,
                ..SynthParams::new(n)
            });
            assert_eq!(s.atom_count(), n);
            assert_eq!(s.topology.validate(), Ok(()));
            assert_eq!(s.topology.residue_count(), n.div_ceil(ATOMS_PER_RESIDUE));
            assert!(s.topology.element.iter().all(|e| !e.is_unknown()));
        }
    }

    #[test]
    fn deterministic_for_same_seed() {
        let a = protein_like(&SynthParams::new(5000));
        let b = protein_like(&SynthParams::new(5000));
        assert_eq!(a.frame(0), b.frame(0));
        let c = protein_like(&SynthParams {
            seed: 7,
            ..SynthParams::new(5000)
        });
        assert_ne!(a.frame(0), c.frame(0));
    }

    #[test]
    fn packed_at_roughly_protein_density() {
        let n = 100_000;
        let s = protein_like(&SynthParams::new(n));
        let (_, r) = s.frame(0).bounding_sphere().unwrap();
        let expected = packing_radius(n);
        assert!(
            r > expected * 0.7 && r < expected * 1.3,
            "radius {r} vs expected {expected}"
        );
    }

    #[test]
    fn chain_names_follow_mmcif_convention() {
        let mut s = String::new();
        let names: Vec<String> = [0, 1, 25, 26, 27, 51, 52, 701, 702]
            .iter()
            .map(|&i| {
                chain_name(i, &mut s);
                s.clone()
            })
            .collect();
        assert_eq!(names, ["A", "B", "Z", "AA", "AB", "AZ", "BA", "ZZ", "AAA"]);
    }
}
