//! The non-polymer content a cartoon or tube draws beside the polymer
//! (`vv_core::Companion`): which kinds a rep shows, how each is styled,
//! and the pieces of geometry that carry them.

use super::*;
use vv_core::Companion;

const ALL: [Companion; 7] = [
    Companion::Ligand,
    Companion::Additive,
    Companion::Ion,
    Companion::Glycan,
    Companion::BoundLipid,
    Companion::MembraneLipid,
    Companion::Water,
];

/// Radius of glycan and bound-lipid sticks, thinner than a licorice rep so
/// a large glycan stays legible beside the ribbon.
const THIN: f32 = 0.15;
/// Ions draw at this fraction of their van der Waals radius: readable as
/// spheres without burying the site they sit in.
const ION_SCALE: f32 = 0.5;

/// The `repopt` name that switches `kind` on or off.
fn option_name(kind: Companion) -> &'static str {
    match kind {
        Companion::Ligand => "ligands",
        Companion::Additive => "additives",
        Companion::Ion => "ions",
        Companion::Glycan => "glycans",
        Companion::BoundLipid | Companion::MembraneLipid => "lipids",
        Companion::Water => "water",
    }
}

fn shown(rep: &Rep, kind: Companion) -> bool {
    rep.option(option_name(kind)).unwrap_or(0.0) >= 0.5
}

/// How `kind` is drawn: ligands, cofactors and additives as ball-and-stick,
/// ions as spheres, glycans and bound lipids as thin sticks, membrane
/// lipids as lines, water as ball-and-stick (an O with its Hs, when it has
/// any).
pub(crate) fn style(kind: Companion) -> (GpuRepresentation, AtomSizes) {
    let sticks = AtomSizes {
        radius_scale: 0.0,
        radius_offset: THIN,
        bond_radius: THIN,
    };
    match kind {
        Companion::Ligand | Companion::Additive | Companion::Water => (
            GpuRepresentation::BallAndStick,
            AtomSizes::of(GpuRepresentation::BallAndStick),
        ),
        Companion::Ion => (
            GpuRepresentation::Spacefill,
            AtomSizes {
                radius_scale: ION_SCALE,
                ..AtomSizes::of(GpuRepresentation::Spacefill)
            },
        ),
        Companion::Glycan | Companion::BoundLipid => (GpuRepresentation::Sticks, sticks),
        Companion::MembraneLipid => (
            GpuRepresentation::Lines,
            AtomSizes::of(GpuRepresentation::Lines),
        ),
    }
}

/// The atoms of every kind `rep` shows, kept to `keep`, with their bonds.
/// Shared by the viewport and the path tracer.
pub(crate) fn subsets(
    structure: &Structure,
    bonds: &[[u32; 2]],
    rep: &Rep,
    keep: &Option<Vec<bool>>,
) -> Vec<(Companion, LigandSubset)> {
    let topology = &structure.topology;
    let kinds = vv_core::companion::companions(topology);
    let shown_kinds: Vec<Companion> = ALL.into_iter().filter(|&k| shown(rep, k)).collect();
    let mut masks: Vec<Vec<bool>> = vec![vec![false; topology.atom_count()]; shown_kinds.len()];
    for (r, kind) in kinds.iter().enumerate() {
        let Some(slot) = shown_kinds.iter().position(|k| Some(*k) == *kind) else {
            continue;
        };
        for a in topology.residues[r].atoms.clone() {
            masks[slot][a as usize] = keep.as_ref().is_none_or(|k| k[a as usize]);
        }
    }
    shown_kinds
        .into_iter()
        .zip(masks)
        .filter(|(_, mask)| mask.contains(&true))
        .map(|(kind, mask)| (kind, subset_of(&mask, bonds)))
        .collect()
}

struct Piece {
    kind: Companion,
    derived: Derived,
}

/// A rep's companions: one [`Derived`] per kind shown.
#[derive(Default)]
pub(super) struct Companions {
    pieces: Vec<Piece>,
}

impl Companions {
    pub(super) fn build(
        ctx: &GpuContext,
        renderer: &Renderer,
        loaded: &LoadedStructure,
        rep: &Rep,
        keep: &Option<Vec<bool>>,
        frame: usize,
        all_colors: &[u32],
        adjacency: Option<&vv_core::Adjacency>,
    ) -> Self {
        let pieces = subsets(&loaded.structure, &loaded.bonds.pairs, rep, keep)
            .into_iter()
            .map(|(kind, subset)| {
                let (representation, sizes) = style(kind);
                let strand_input = adjacency
                    .filter(|_| representation != GpuRepresentation::Lines)
                    .filter(|_| sizes.bond_radius > 0.0)
                    .map(|adjacency| StrandInput {
                        bonds: &loaded.bonds,
                        adjacency,
                        separation: sizes.bond_radius * STRAND_SEPARATION_SCALE,
                    });
                let derived = Derived::atoms(
                    ctx,
                    renderer,
                    &loaded.structure,
                    frame,
                    subset,
                    all_colors,
                    strand_input,
                );
                Piece { kind, derived }
            })
            .collect();
        Companions { pieces }
    }

    pub(super) fn set_frame(&self, ctx: &GpuContext, structure: &Structure, frame: usize) {
        for piece in &self.pieces {
            piece.derived.set_frame(ctx, structure, frame);
        }
    }

    /// Each piece with the representation and sizes it draws with.
    pub(super) fn draw(&self) -> impl Iterator<Item = (&Derived, GpuRepresentation, AtomSizes)> {
        self.pieces.iter().map(|p| {
            let (representation, sizes) = style(p.kind);
            (&p.derived, representation, sizes)
        })
    }

    /// Occlusion spheres at the pieces' drawn radii.
    pub(super) fn proxies(
        &self,
        positions: &[Vec3],
        elements: &[vv_core::Element],
        out: &mut Vec<(Vec3, f32)>,
    ) {
        for (derived, _, sizes) in self.draw() {
            out.extend(derived.atom_map.iter().filter_map(|&a| {
                let radius = sizes.atom_radius(elements[a as usize].vdw_radius());
                (radius > 0.0).then_some((positions[a as usize], radius))
            }));
        }
    }
}

/// Adds every companion `rep` shows to the path tracer's scene, drawn as
/// the viewport draws them.
pub(crate) fn push_traced(
    traced: &mut vv_render::path_trace::TraceScene,
    loaded: &LoadedStructure,
    rep: &Rep,
    keep: &Option<Vec<bool>>,
    positions: &[Vec3],
    colors: &[u32],
    material: vv_render::Material,
) {
    const LINE_RADIUS: f32 = 0.05;
    for (kind, subset) in subsets(&loaded.structure, &loaded.bonds.pairs, rep, keep) {
        let (representation, sizes) = style(kind);
        let radii: Vec<f32> = loaded
            .structure
            .topology
            .element
            .iter()
            .map(|e| sizes.atom_radius(e.vdw_radius()))
            .collect();
        let bond_radius = match representation {
            GpuRepresentation::Lines => LINE_RADIUS,
            _ => sizes.bond_radius,
        };
        traced.push_atoms(
            positions,
            &radii,
            colors,
            subset.atoms.iter().map(|&a| a as usize),
            subset.real_pairs.iter().copied(),
            bond_radius,
            material,
        );
    }
}
