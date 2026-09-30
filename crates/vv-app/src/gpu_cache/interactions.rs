//! The contact overlay (`vv_core::interactions`): each contact drawn as a
//! run of short cylinders, so the dashes are depth-tested with the rest of
//! the scene and appear in path-traced renders.

use std::collections::BTreeSet;

use super::*;
use vv_core::interactions::{self, InteractionKind};

/// Rep id of the overlay's pick sources: matches no rep, so a structure
/// row's highlight never lights it.
pub(crate) const OVERLAY_REP: RepId = RepId(u32::MAX);

const DASH_LENGTH: f32 = 0.3;
const DASH_GAP: f32 = 0.25;
pub(crate) const DASH_RADIUS: f32 = 0.07;

pub(crate) fn color(kind: InteractionKind) -> u32 {
    let [r, g, b] = match kind {
        InteractionKind::Hbond => [0x3C, 0x8C, 0xE6],
        InteractionKind::Metal => [0xC8, 0x78, 0x14],
        InteractionKind::SaltBridge => [0xC8, 0x32, 0x96],
    };
    vv_render::color::rgba(r, g, b)
}

/// One dash of a contact: its two ends, and the contact's real atoms.
pub(crate) struct Dash {
    pub from: Vec3,
    pub to: Vec3,
    pub atoms: [u32; 2],
    pub kind: InteractionKind,
}

fn contacts(loaded: &LoadedStructure, frame: usize, kind: InteractionKind) -> Vec<[u32; 2]> {
    let structure = &loaded.structure;
    let coords = structure.frame(frame);
    let positions = coords.positions();
    let shown = loaded.shown_atoms();
    let visible = |a: u32| shown.as_ref().is_none_or(|s| s.contains(a as usize));
    let topology = &structure.topology;
    match kind {
        InteractionKind::Hbond => {
            interactions::hydrogen_bonds(topology, &loaded.bonds, positions, &visible)
        }
        InteractionKind::Metal => interactions::metal_coordination(topology, positions, &visible),
        InteractionKind::SaltBridge => interactions::salt_bridges(topology, positions, &visible),
    }
}

/// The dashes of `pair`, centred along the line between its atoms.
fn dashes_of(a: Vec3, b: Vec3, atoms: [u32; 2], kind: InteractionKind) -> Vec<Dash> {
    let length = a.distance(b);
    let step = DASH_LENGTH + DASH_GAP;
    let count = ((length + DASH_GAP) / step).floor().max(1.0) as usize;
    let used = count as f32 * step - DASH_GAP;
    let start = ((length - used) * 0.5).max(0.0);
    let along = (b - a) / length.max(f32::EPSILON);
    (0..count)
        .map(|k| {
            let s = start + k as f32 * step;
            Dash {
                from: a + along * s,
                to: a + along * (s + DASH_LENGTH).min(length),
                atoms,
                kind,
            }
        })
        .collect()
}

/// Every dash of every overlay `loaded` has on, at `frame`.
pub(crate) fn dashes(loaded: &LoadedStructure, frame: usize) -> Vec<Dash> {
    let coords = loaded.structure.frame(frame);
    let positions = coords.positions();
    loaded
        .interactions
        .iter()
        .flat_map(|&kind| {
            contacts(loaded, frame, kind)
                .into_iter()
                .flat_map(move |[a, b]| {
                    dashes_of(positions[a as usize], positions[b as usize], [a, b], kind)
                })
        })
        .collect()
}

/// The overlay on the GPU, for the kinds and frame it was built at.
pub(super) struct Overlay {
    key: (BTreeSet<InteractionKind>, usize),
    gpu: GpuStructure,
    bindings: PageBindings,
    atom_map: Arc<Vec<u32>>,
    bond_atoms: Arc<Vec<[u32; 2]>>,
}

impl Overlay {
    /// `slot` brought up to date with `loaded`'s overlays at `frame`.
    pub(super) fn sync(
        slot: &mut Option<Overlay>,
        ctx: &GpuContext,
        renderer: &Renderer,
        loaded: &LoadedStructure,
        frame: usize,
    ) {
        let key = (loaded.interactions.clone(), frame);
        if loaded.interactions.is_empty() {
            *slot = None;
        } else if slot.as_ref().is_none_or(|o| o.key != key) {
            *slot = Some(Overlay::build(ctx, renderer, loaded, frame, key));
        }
    }

    fn build(
        ctx: &GpuContext,
        renderer: &Renderer,
        loaded: &LoadedStructure,
        frame: usize,
        key: (BTreeSet<InteractionKind>, usize),
    ) -> Self {
        let dashes = dashes(loaded, frame);
        let positions: Vec<Vec3> = dashes.iter().flat_map(|d| [d.from, d.to]).collect();
        let colors: Vec<u32> = dashes.iter().flat_map(|d| [color(d.kind); 2]).collect();
        let bonds: Vec<[u32; 2]> = (0..dashes.len() as u32)
            .map(|k| [2 * k, 2 * k + 1])
            .collect();
        let radii = vec![0.0; positions.len()];
        let gpu = GpuStructure::from_parts(ctx, &positions, &radii, &colors, &bonds);
        Overlay {
            key,
            bindings: renderer.bind(&gpu),
            gpu,
            atom_map: Arc::new(dashes.iter().flat_map(|d| d.atoms).collect()),
            bond_atoms: Arc::new(dashes.iter().map(|d| d.atoms).collect()),
        }
    }

    pub(super) fn draw(&self, id: StructureId) -> (DrawSource, DrawItem<'_>) {
        let source = DrawSource {
            id,
            rep: OVERLAY_REP,
            atom_map: Some(self.atom_map.clone()),
            bond_atoms: Some(self.bond_atoms.clone()),
        };
        let item = DrawItem {
            structure: &self.gpu,
            bindings: &self.bindings,
            representation: GpuRepresentation::BallAndStick,
            sizes: AtomSizes {
                radius_scale: 0.0,
                radius_offset: 0.0,
                bond_radius: DASH_RADIUS,
            },
            material: vv_render::Material::default(),
        };
        (source, item)
    }
}

/// Adds every dash to the path tracer's scene.
pub(crate) fn push_traced(
    traced: &mut vv_render::path_trace::TraceScene,
    loaded: &LoadedStructure,
    frame: usize,
) {
    let material = traced.materials.len() as u32;
    traced.materials.push(vv_render::Material::default());
    traced
        .cylinders
        .extend(
            dashes(loaded, frame)
                .into_iter()
                .map(|d| vv_render::path_trace::TraceCylinder {
                    a: d.from,
                    b: d.to,
                    radius: DASH_RADIUS,
                    color_a: color(d.kind),
                    color_b: color(d.kind),
                    material,
                }),
        );
}
