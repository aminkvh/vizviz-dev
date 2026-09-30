//! A cartoon's nucleotide bases (`vv_core::bases`): slabs on the ring
//! atoms with a stick to the sugar, one rung per base pair, or the plain
//! stick to the pairing atom.

use super::*;
use vv_core::bases::{self, Base};
use vv_core::PolytopeMesh;
use vv_render::renderer::CartoonBindings;
use vv_render::GlycanGpu;

/// Radius of the connecting sticks and ladder rungs.
const STICK: f32 = 0.2;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum BaseStyle {
    Stick,
    Plate,
    Ladder,
}

impl BaseStyle {
    pub(crate) fn of(rep: &Rep) -> Self {
        match rep.option("bases").unwrap_or(1.0) as usize {
            0 => BaseStyle::Stick,
            2 => BaseStyle::Ladder,
            _ => BaseStyle::Plate,
        }
    }
}

pub(crate) fn stick_sizes() -> AtomSizes {
    AtomSizes {
        radius_scale: 0.0,
        radius_offset: STICK,
        bond_radius: STICK,
    }
}

/// What the bases of `keep`'s nucleotides are made of at `frame`: slab
/// triangles, and the sticks (as a subset of the structure's atoms).
pub(crate) struct BaseParts {
    pub plates: PolytopeMesh,
    pub sticks: LigandSubset,
}

fn unpack(color: u32) -> [u8; 3] {
    [color as u8, (color >> 8) as u8, (color >> 16) as u8]
}

/// The bases of the shown conformer whose ring atoms all pass `keep`.
fn kept_bases(loaded: &LoadedStructure, keep: &Option<Vec<bool>>, positions: &[Vec3]) -> Vec<Base> {
    let shown = loaded.shown_atoms();
    let visible = |a: u32| shown.as_ref().is_none_or(|s| s.contains(a as usize));
    bases::bases(&loaded.structure.topology, positions, &visible)
        .into_iter()
        .filter(|b| {
            keep.as_ref()
                .is_none_or(|k| b.ring.iter().all(|&a| k[a as usize]))
        })
        .collect()
}

fn sticks_between(structure: &Structure, rungs: &[[u32; 2]]) -> LigandSubset {
    let mut ends = vec![false; structure.atom_count()];
    for &[a, b] in rungs {
        ends[a as usize] = true;
        ends[b as usize] = true;
    }
    subset_of(&ends, rungs)
}

/// The legacy stick: sugar trace atom to the pairing atom, on the shown
/// conformers.
pub(super) fn stick_rungs(loaded: &LoadedStructure, keep: &Option<Vec<bool>>) -> Vec<[u32; 2]> {
    loaded
        .as_shown(&vv_core::nucleic_ladder(&loaded.structure.topology))
        .into_iter()
        .filter(|[a, b]| {
            keep.as_ref()
                .is_none_or(|k| k[*a as usize] && k[*b as usize])
        })
        .collect()
}

/// Slab colours are `colors` at each base's glycosidic atom.
pub(crate) fn parts(
    loaded: &LoadedStructure,
    keep: &Option<Vec<bool>>,
    frame: usize,
    colors: &[u32],
    style: BaseStyle,
) -> BaseParts {
    let coords = loaded.structure.frame(frame);
    let positions = coords.positions();
    let mut plates = PolytopeMesh::default();
    let rungs = match style {
        BaseStyle::Stick => stick_rungs(loaded, keep),
        BaseStyle::Plate => {
            let all = kept_bases(loaded, keep, positions);
            for base in &all {
                let color = unpack(colors[base.glycosidic as usize]);
                bases::push_plate(base, positions, color, &mut plates);
            }
            bases::connectors(&all)
        }
        BaseStyle::Ladder => {
            let all = kept_bases(loaded, keep, positions);
            bases::ladder(&all, &bases::pairs(&all, positions))
        }
    };
    BaseParts {
        plates,
        sticks: sticks_between(&loaded.structure, &rungs),
    }
}

struct Plates {
    gpu: GlycanGpu,
    bindings: CartoonBindings,
}

/// A cartoon's bases on the GPU.
#[derive(Default)]
pub(super) struct Bases {
    plates: Option<Plates>,
    sticks: Option<Derived>,
}

impl Bases {
    pub(super) fn build(
        ctx: &GpuContext,
        renderer: &Renderer,
        loaded: &LoadedStructure,
        keep: &Option<Vec<bool>>,
        frame: usize,
        colors: &[u32],
        style: BaseStyle,
    ) -> Result<Self, OutOfGpuMemory> {
        let BaseParts { plates, sticks } = parts(loaded, keep, frame, colors, style);
        let plates = if plates.positions.is_empty() {
            None
        } else {
            let gpu = GlycanGpu::upload(ctx, &plates)?;
            let bindings = renderer.bind_glycan(&gpu);
            Some(Plates { gpu, bindings })
        };
        // Sugar-to-base links are synthetic, never in the `BondTable`.
        let sticks = (!sticks.atoms.is_empty()).then(|| {
            Derived::atoms(
                ctx,
                renderer,
                &loaded.structure,
                frame,
                sticks,
                colors,
                None,
            )
        });
        Ok(Bases { plates, sticks })
    }

    pub(super) fn sticks(&self) -> Option<&Derived> {
        self.sticks.as_ref()
    }

    /// The slab mesh and its bindings, when there are slabs to draw.
    pub(super) fn plates(&self) -> Option<(&GlycanGpu, &CartoonBindings)> {
        self.plates.as_ref().map(|p| (&p.gpu, &p.bindings))
    }
}

/// Adds the bases of `keep`'s nucleotides to the path tracer's scene.
pub(crate) fn push_traced(
    traced: &mut vv_render::path_trace::TraceScene,
    loaded: &LoadedStructure,
    keep: &Option<Vec<bool>>,
    frame: usize,
    colors: &[u32],
    style: BaseStyle,
    material: vv_render::Material,
) {
    let BaseParts { plates, sticks } = parts(loaded, keep, frame, colors, style);
    traced.push_glycan_mesh(&plates, material);
    let sizes = stick_sizes();
    let radii: Vec<f32> = loaded
        .structure
        .topology
        .element
        .iter()
        .map(|e| sizes.atom_radius(e.vdw_radius()))
        .collect();
    traced.push_atoms(
        loaded.structure.frame(frame).positions(),
        &radii,
        colors,
        sticks.atoms.iter().map(|&a| a as usize),
        sticks.real_pairs.iter().copied(),
        sizes.bond_radius,
        material,
    );
}
