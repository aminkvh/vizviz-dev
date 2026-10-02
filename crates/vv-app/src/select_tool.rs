//! The Home ▸ Select tool: drag a box, circle or lasso in the viewport and
//! select the atoms that project inside it, as atoms, whole residues or
//! whole chains. The state is view state (like the mouse mode); the result
//! goes through `Command::Select`, so it is undoable and shows in the
//! Selections panel like any other selection.
//!
//! The test is a projection test, not a visibility test: atoms hidden
//! behind others are selected too (there is no cheap occlusion query), but
//! atoms in no visible rep are not.

use egui::{Color32, Pos2, Stroke};
use rayon::prelude::*;
use vv_core::fixedbitset::{Block, FixedBitSet};
use vv_core::glam::{Mat4, Vec3};
use vv_core::{ResidueRec, Topology};
use vv_scene::{Command, Mask, StructureId};

use crate::gpu_cache::select_atoms;
use crate::ui::{project_to_pixel, AppUi};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SelectShape {
    /// Clicking picks; dragging orbits (the ordinary mouse mode).
    #[default]
    Click,
    Box,
    Circle,
    Lasso,
}

impl SelectShape {
    pub const ALL: [Self; 4] = [Self::Click, Self::Box, Self::Circle, Self::Lasso];

    pub fn word(self) -> &'static str {
        match self {
            Self::Click => "click",
            Self::Box => "box",
            Self::Circle => "circle",
            Self::Lasso => "lasso",
        }
    }

    pub fn from_word(word: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|s| s.word() == word)
    }

    pub fn next(self) -> Self {
        cycle(&Self::ALL, self)
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SelectLevel {
    #[default]
    Atom,
    Residue,
    Chain,
    /// A connected component of the bond graph.
    Molecule,
}

impl SelectLevel {
    pub const ALL: [Self; 4] = [Self::Atom, Self::Residue, Self::Chain, Self::Molecule];

    pub fn word(self) -> &'static str {
        match self {
            Self::Atom => "atom",
            Self::Residue => "residue",
            Self::Chain => "chain",
            Self::Molecule => "molecule",
        }
    }

    pub fn from_word(word: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|l| l.word() == word)
    }

    pub fn next(self) -> Self {
        cycle(&Self::ALL, self)
    }
}

fn cycle<T: Copy + PartialEq>(all: &[T], now: T) -> T {
    let at = all.iter().position(|&x| x == now).unwrap_or(0);
    all[(at + 1) % all.len()]
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SelectTool {
    pub shape: SelectShape,
    pub level: SelectLevel,
    /// The drawing shape the tool button switches back to.
    remembered: SelectShape,
}

impl Default for SelectTool {
    fn default() -> Self {
        Self {
            shape: SelectShape::Click,
            level: SelectLevel::Atom,
            remembered: SelectShape::Box,
        }
    }
}

/// How a drag combines with the current selection.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Combine {
    Replace,
    Add,
    Subtract,
}

/// Seconds a press on the tool button is held before its flyout opens.
pub const HOLD_SECS: f64 = 0.4;

/// What a press on the grouped tool button does.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ToolPress {
    Activate,
    OpenFlyout,
}

/// A press on the button's corner marker, or one held past `HOLD_SECS`,
/// opens the flyout; a shorter press elsewhere is a click.
pub fn classify_press(held_secs: f64, on_corner: bool) -> ToolPress {
    if on_corner || held_secs >= HOLD_SECS {
        ToolPress::OpenFlyout
    } else {
        ToolPress::Activate
    }
}

impl SelectTool {
    /// Whether a left drag draws a shape instead of orbiting.
    pub fn draws(self) -> bool {
        self.shape != SelectShape::Click
    }

    pub fn set_shape(&mut self, shape: SelectShape) {
        self.shape = shape;
        if shape != SelectShape::Click {
            self.remembered = shape;
        }
    }

    /// The shape a click on the tool button switches to: back to plain
    /// picking when a shape is on, else the last shape used.
    pub fn toggled_shape(self) -> SelectShape {
        if self.draws() {
            SelectShape::Click
        } else {
            self.remembered
        }
    }

    /// Applies the words after `selectmode`: `shape NAME` or `level NAME`,
    /// where NAME may also be `next` (cycle) and, for a shape, `toggle`.
    pub fn set(&mut self, words: &str) -> Result<(), ()> {
        let mut it = words.split_whitespace();
        match (it.next(), it.next(), it.next()) {
            (Some("shape"), Some(name), None) => {
                let shape = match name {
                    "next" => self.shape.next(),
                    "toggle" => self.toggled_shape(),
                    _ => SelectShape::from_word(name).ok_or(())?,
                };
                self.set_shape(shape);
            }
            (Some("level"), Some(name), None) => {
                self.level = match name {
                    "next" => self.level.next(),
                    _ => SelectLevel::from_word(name).ok_or(())?,
                };
            }
            _ => return Err(()),
        }
        Ok(())
    }

    pub fn shape_word(self) -> &'static str {
        self.shape.word()
    }

    pub fn level_word(self) -> &'static str {
        self.level.word()
    }

    pub fn drag_label(self) -> String {
        format!("Select {} by {}", self.shape_word(), self.level_word())
    }
}

/// A screen region in viewport-local pixels.
#[derive(Clone, Debug, PartialEq)]
pub enum Region {
    Box { min: Pos2, max: Pos2 },
    Circle { center: Pos2, radius: f32 },
    Lasso(Vec<Pos2>),
}

impl Region {
    /// The region a drag traced: a box between its ends, a circle centered
    /// where it began, or the freehand polygon itself. `None` when the
    /// drag is too short to enclose anything.
    pub fn from_path(shape: SelectShape, path: &[Pos2]) -> Option<Self> {
        let (&first, &last) = (path.first()?, path.last()?);
        match shape {
            SelectShape::Click => None,
            SelectShape::Box => Some(Region::Box {
                min: first.min(last),
                max: first.max(last),
            }),
            SelectShape::Circle => Some(Region::Circle {
                center: first,
                radius: first.distance(last),
            }),
            SelectShape::Lasso => (path.len() >= 3).then(|| Region::Lasso(path.to_vec())),
        }
    }

    pub fn contains(&self, x: f32, y: f32) -> bool {
        match self {
            Region::Box { min, max } => x >= min.x && x <= max.x && y >= min.y && y <= max.y,
            Region::Circle { center, radius } => {
                let (dx, dy) = (x - center.x, y - center.y);
                dx * dx + dy * dy <= radius * radius
            }
            Region::Lasso(points) => polygon_contains(points, x, y),
        }
    }

    /// Corners of the axis-aligned box around the region, for a cheap
    /// reject before the exact test.
    fn bounds(&self) -> (Pos2, Pos2) {
        match self {
            Region::Box { min, max } => (*min, *max),
            Region::Circle { center, radius } => (
                *center - egui::vec2(*radius, *radius),
                *center + egui::vec2(*radius, *radius),
            ),
            Region::Lasso(points) => points.iter().fold(
                (
                    Pos2::new(f32::INFINITY, f32::INFINITY),
                    Pos2::new(f32::NEG_INFINITY, f32::NEG_INFINITY),
                ),
                |(lo, hi), p| (lo.min(*p), hi.max(*p)),
            ),
        }
    }
}

/// Even-odd ray casting; the polygon closes from its last point.
fn polygon_contains(points: &[Pos2], x: f32, y: f32) -> bool {
    let Some(&last) = points.last() else {
        return false;
    };
    let mut prev = last;
    let mut inside = false;
    for &p in points {
        if (p.y > y) != (prev.y > y) && x < (prev.x - p.x) * (y - p.y) / (prev.y - p.y) + p.x {
            inside = !inside;
        }
        prev = p;
    }
    inside
}

const BLOCK_BITS: usize = Block::BITS as usize;
const BLOCKS_PER_TASK: usize = 1024;

/// Atoms whose projection lands in `region`: one parallel pass over
/// `positions` writing whole mask blocks, with no per-atom allocation.
pub fn atoms_inside(
    positions: &[Vec3],
    view_proj: Mat4,
    size: (f32, f32),
    region: &Region,
) -> FixedBitSet {
    let mut mask = FixedBitSet::with_capacity(positions.len());
    let (lo, hi) = region.bounds();
    let inside = |i: usize| {
        project_to_pixel(view_proj, size.0, size.1, positions[i]).is_some_and(|(x, y, _)| {
            x >= lo.x && x <= hi.x && y >= lo.y && y <= hi.y && region.contains(x, y)
        })
    };
    mask.as_mut_slice()
        .par_chunks_mut(BLOCKS_PER_TASK)
        .enumerate()
        .for_each(|(task, blocks)| {
            for (j, block) in blocks.iter_mut().enumerate() {
                let start = (task * BLOCKS_PER_TASK + j) * BLOCK_BITS;
                let end = (start + BLOCK_BITS).min(positions.len());
                *block = (start..end).fold(0, |b, i| b | (Block::from(inside(i)) << (i - start)));
            }
        });
    mask
}

/// Grows `mask` to every atom of each connected fragment it touches
/// (`fragments[atom]` is the atom's fragment id, `BondTable::fragments`).
pub fn expand_fragments(mask: &mut FixedBitSet, fragments: &[u32]) {
    let mut hit = vec![false; fragments.iter().map(|&f| f as usize + 1).max().unwrap_or(0)];
    for atom in mask.ones() {
        hit[fragments[atom] as usize] = true;
    }
    for (atom, &fragment) in fragments.iter().enumerate() {
        if hit[fragment as usize] {
            mask.insert(atom);
        }
    }
}

/// Grows `mask` to every atom of each residue or chain it touches; the
/// other levels leave it alone.
pub fn expand(mask: &mut FixedBitSet, topology: &Topology, level: SelectLevel) {
    let key: fn(&ResidueRec, usize) -> usize = match level {
        SelectLevel::Atom | SelectLevel::Molecule => return,
        SelectLevel::Residue => |_, r| r,
        SelectLevel::Chain => |res, _| res.chain as usize,
    };
    let mut hit = vec![false; topology.residues.len().max(topology.chains.len())];
    for atom in mask.ones() {
        let r = topology.residue_index[atom] as usize;
        hit[key(&topology.residues[r], r)] = true;
    }
    for (r, res) in topology.residues.iter().enumerate() {
        if hit[key(res, r)] {
            mask.insert_range(res.atoms.start as usize..res.atoms.end as usize);
        }
    }
}

fn combine(current: Option<&Mask>, hits: FixedBitSet, how: Combine) -> FixedBitSet {
    let (Some(current), Combine::Add | Combine::Subtract) = (current, how) else {
        return hits;
    };
    let mut out = (**current).clone();
    match how {
        Combine::Subtract => out.difference_with(&hits),
        _ => out.union_with(&hits),
    }
    out
}

/// Atoms some visible rep selects and the altloc policy displays.
fn drawn_atoms(loaded: &vv_scene::LoadedStructure, frame: usize) -> FixedBitSet {
    let n = loaded.structure.atom_count();
    let mut drawn = FixedBitSet::with_capacity(n);
    for rep in loaded.reps.iter().filter(|r| r.visible) {
        match select_atoms(loaded, rep, frame) {
            Ok(None) => drawn.insert_range(..),
            Ok(Some(list)) => list.iter().for_each(|&a| drawn.insert(a as usize)),
            Err(_) => {}
        }
    }
    drawn
}

/// The hits among a structure's atoms that some visible rep draws.
fn visible_hits(
    loaded: &vv_scene::LoadedStructure,
    view_proj: Mat4,
    size: (f32, f32),
    region: &Region,
    level: SelectLevel,
) -> FixedBitSet {
    let frame = loaded
        .frame
        .min(loaded.structure.frame_count().saturating_sub(1));
    let coords = loaded.structure.frame(frame);
    let positions = coords.positions();
    let mut hits = atoms_inside(positions, view_proj, size, region);
    hits.intersect_with(&drawn_atoms(loaded, frame));
    grow_to_level(loaded, &mut hits, level);
    hits
}

/// Grows `hits` to the residues, chains or molecules of `level`.
pub fn grow_to_level(
    loaded: &vv_scene::LoadedStructure,
    hits: &mut FixedBitSet,
    level: SelectLevel,
) {
    match level {
        SelectLevel::Molecule => {
            expand_fragments(hits, &loaded.bonds.fragments(loaded.structure.atom_count()))
        }
        _ => expand(hits, &loaded.structure.topology, level),
    }
}

/// What a click on `atom` picks at `level`, without the hidden conformers.
pub fn picked_atoms(
    loaded: &vv_scene::LoadedStructure,
    atom: u32,
    level: SelectLevel,
) -> FixedBitSet {
    let mut hits = FixedBitSet::with_capacity(loaded.structure.atom_count());
    hits.insert(atom as usize);
    grow_to_level(loaded, &mut hits, level);
    if let Some(shown) = loaded.shown_atoms() {
        hits.intersect_with(&shown);
    }
    hits
}

/// How a click or drag combines with the selection: Shift adds, Alt
/// subtracts (Ctrl subtracts drags too, but a Ctrl click is the measure
/// chain's).
pub fn combine_for(modifiers: egui::Modifiers, drag: bool) -> Combine {
    if modifiers.alt || (drag && modifiers.command) {
        Combine::Subtract
    } else if modifiers.shift {
        Combine::Add
    } else {
        Combine::Replace
    }
}

impl AppUi<'_> {
    /// Selects what `region` (viewport-local pixels; `size` is the
    /// viewport) encloses, at the tool's granularity, in the visible
    /// structure with the most atoms inside.
    pub(crate) fn select_in_region(&mut self, region: &Region, size: (f32, f32), how: Combine) {
        let aspect = size.0 / size.1.max(1.0);
        let view_proj = self.camera.proj(aspect) * self.camera.view();
        let level = self.view.select_tool.level;
        let best = self
            .scene
            .structures()
            .filter(|(_, loaded)| loaded.visible)
            .map(|(id, loaded)| (id, visible_hits(loaded, view_proj, size, region, level)))
            .max_by_key(|(_, hits)| hits.count_ones(..));
        if let Some((id, hits)) = best {
            self.apply_hits(id, hits, how);
        }
    }

    /// A click on `atom` at the tool's level.
    pub(crate) fn pick_at_level(&mut self, id: StructureId, atom: u32, how: Combine) {
        let level = self.view.select_tool.level;
        let Some(loaded) = self.scene.structure(id) else {
            return;
        };
        let hits = picked_atoms(loaded, atom, level);
        self.apply_hits(id, hits, how);
    }

    /// Grows the active selection to whole residues, chains or molecules
    /// when the pick level is raised, so the level and the selection agree.
    pub(crate) fn snap_selection_to_level(&mut self) {
        let level = self.view.select_tool.level;
        let Some(active) = self.scene.active_selection() else {
            return;
        };
        let Some(loaded) = self.scene.structure(active.structure) else {
            return;
        };
        let mut grown = (*active.mask).clone();
        grow_to_level(loaded, &mut grown, level);
        if grown != *active.mask {
            let id = active.structure;
            self.dispatch(Command::Select {
                id,
                mask: std::sync::Arc::new(grown),
            });
        }
    }

    fn apply_hits(&mut self, id: StructureId, hits: FixedBitSet, how: Combine) {
        let same = self.scene.active_selection().filter(|a| a.structure == id);
        let mask = combine(same.map(|a| &a.mask), hits, how);
        self.pick_order.clear();
        if mask.is_clear() {
            self.dispatch(Command::ClearSelection);
        } else {
            self.dispatch(Command::Select {
                id,
                mask: std::sync::Arc::new(mask),
            });
        }
    }
}

/// The drag's outline over the viewport at `origin`, in the primary color.
pub fn paint_region(painter: &egui::Painter, origin: Pos2, region: &Region, primary: Color32) {
    let stroke = Stroke::new(1.5, primary);
    let fill = primary.gamma_multiply(0.12);
    let at = |p: &Pos2| origin + p.to_vec2();
    match region {
        Region::Box { min, max } => {
            let rect = egui::Rect::from_min_max(at(min), at(max));
            painter.rect(rect, 0.0, fill, stroke, egui::StrokeKind::Middle);
        }
        Region::Circle { center, radius } => {
            painter.circle(at(center), *radius, fill, stroke);
        }
        Region::Lasso(points) => {
            painter.add(egui::Shape::closed_line(
                points.iter().map(at).collect(),
                stroke,
            ));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use vv_core::{ChainRec, Interner};

    fn p(x: f32, y: f32) -> Pos2 {
        Pos2::new(x, y)
    }

    #[test]
    fn a_box_takes_either_drag_direction() {
        let region = Region::from_path(SelectShape::Box, &[p(10.0, 20.0), p(0.0, 0.0)]).unwrap();
        assert!(region.contains(5.0, 10.0));
        assert!(!region.contains(11.0, 10.0));
    }

    #[test]
    fn a_circle_is_centered_where_the_drag_began() {
        let region = Region::from_path(SelectShape::Circle, &[p(0.0, 0.0), p(3.0, 4.0)]).unwrap();
        assert!(region.contains(0.0, 4.9));
        assert!(!region.contains(0.0, 5.1));
    }

    #[test]
    fn a_lasso_handles_a_concave_polygon() {
        let points = [
            p(0.0, 0.0),
            p(10.0, 0.0),
            p(10.0, 10.0),
            p(5.0, 4.0),
            p(0.0, 10.0),
        ];
        let region = Region::from_path(SelectShape::Lasso, &points).unwrap();
        assert!(region.contains(2.0, 2.0));
        assert!(!region.contains(5.0, 8.0), "inside the notch");
        assert!(!region.contains(-1.0, 5.0));
    }

    #[test]
    fn a_lasso_needs_three_points() {
        assert!(Region::from_path(SelectShape::Lasso, &[p(0.0, 0.0), p(1.0, 1.0)]).is_none());
        assert!(Region::from_path(SelectShape::Click, &[p(0.0, 0.0)]).is_none());
    }

    #[test]
    fn atoms_inside_matches_a_serial_test_across_block_boundaries() {
        let n = 3 * BLOCK_BITS + 7;
        let positions: Vec<Vec3> = (0..n)
            .map(|i| Vec3::new(-1.0 + 2.0 * (i as f32 + 0.5) / n as f32, 0.0, 0.5))
            .collect();
        let region = Region::Box {
            min: p(25.0, 0.0),
            max: p(75.0, 100.0),
        };
        let mask = atoms_inside(&positions, Mat4::IDENTITY, (100.0, 100.0), &region);
        for (i, w) in positions.iter().enumerate() {
            let (x, y, _) = project_to_pixel(Mat4::IDENTITY, 100.0, 100.0, *w).unwrap();
            assert_eq!(mask.contains(i), region.contains(x, y), "atom {i}");
        }
        assert!(mask.count_ones(..) > 0);
    }

    #[test]
    fn a_region_never_selects_a_hidden_conformer() {
        use vv_core::altloc::AltlocPolicy;
        use vv_scene::{CommandHistory, Scene};
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/1AKE.pdb");
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(10);
        history
            .dispatch(&mut scene, Command::LoadStructure { path })
            .unwrap();
        let id = scene.structures().next().unwrap().0;
        let everything = Region::Box {
            min: p(-1e6, -1e6),
            max: p(1e6, 1e6),
        };
        let view = Mat4::orthographic_rh(-1e3, 1e3, -1e3, 1e3, -1e3, 1e3);
        let hits = |scene: &Scene| {
            let loaded = scene.structure(id).unwrap();
            visible_hits(loaded, view, (100.0, 100.0), &everything, SelectLevel::Atom)
        };
        let hidden: Vec<usize> = scene
            .structure(id)
            .unwrap()
            .shown_atoms()
            .expect("1AKE has conformers")
            .zeroes()
            .collect();
        assert!(!hidden.is_empty());
        let first = hits(&scene);
        assert!(hidden.iter().all(|&a| !first.contains(a)));
        history
            .dispatch(
                &mut scene,
                Command::SetAltloc {
                    id,
                    policy: AltlocPolicy::All,
                },
            )
            .unwrap();
        assert!(hidden.iter().all(|&a| hits(&scene).contains(a)));
    }

    /// Chain 0 holds residues 0 and 1 (atoms 0..2, 2..4); chain 1 holds
    /// residue 2 (atoms 4..6).
    fn three_residue_topology() -> Topology {
        let mut names = Interner::new();
        let comp = names.intern("ALA");
        let residue = |atoms, chain| ResidueRec {
            atoms,
            chain,
            comp,
            seq_id: 1,
            auth_seq_id: 1,
            ins_code: 0,
            ss: vv_core::SecondaryStructure::Unknown,
        };
        let chain = |residues| ChainRec {
            residues,
            label_asym: comp,
            auth_asym: comp,
            entity: 0,
        };
        Topology {
            residue_index: vec![0, 0, 1, 1, 2, 2],
            residues: vec![residue(0..2, 0), residue(2..4, 0), residue(4..6, 1)],
            chains: vec![chain(0..2), chain(2..3)],
            names,
            ..Topology::default()
        }
    }

    #[test]
    fn granularity_expands_to_residues_and_chains() {
        let topology = three_residue_topology();
        let selected = |level| {
            let mut mask = FixedBitSet::with_capacity(6);
            mask.insert(0);
            expand(&mut mask, &topology, level);
            mask.ones().collect::<Vec<_>>()
        };
        assert_eq!(selected(SelectLevel::Atom), [0]);
        assert_eq!(selected(SelectLevel::Residue), [0, 1]);
        assert_eq!(selected(SelectLevel::Chain), [0, 1, 2, 3]);
    }

    #[test]
    fn a_chain_pick_takes_the_whole_chain_the_atom_is_in() {
        use vv_scene::{Command, CommandHistory, Scene};
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/1AKE.pdb");
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(10);
        history
            .dispatch(&mut scene, Command::LoadStructure { path })
            .unwrap();
        let loaded = scene.structures().next().unwrap().1;
        let top = &loaded.structure.topology;
        let chain = top.chain_of_atom(0);
        let shown = loaded.shown_atoms().expect("1AKE has conformers");
        let expected: usize = top.chains[chain as usize]
            .residues
            .clone()
            .flat_map(|r| top.residues[r as usize].atoms.clone())
            .filter(|&a| shown.contains(a as usize))
            .count();
        let picked = picked_atoms(loaded, 0, SelectLevel::Chain);
        let residue = picked_atoms(loaded, 0, SelectLevel::Residue);
        assert_eq!(picked.count_ones(..), expected);
        assert!(picked.count_ones(..) > residue.count_ones(..));
        assert!(picked.ones().all(|a| top.chain_of_atom(a) == chain));
    }

    #[test]
    fn a_molecule_is_every_atom_of_a_touched_fragment() {
        let fragments = [0, 0, 1, 1, 1, 2];
        let mut mask = FixedBitSet::with_capacity(6);
        mask.insert(3);
        expand_fragments(&mut mask, &fragments);
        assert_eq!(mask.ones().collect::<Vec<_>>(), [2, 3, 4]);
        let mut none = FixedBitSet::with_capacity(6);
        expand_fragments(&mut none, &fragments);
        assert!(none.is_clear());
    }

    #[test]
    fn levels_and_shapes_cycle_through_all_and_wrap() {
        let mut seen = vec![SelectLevel::Atom];
        for _ in 0..4 {
            seen.push(seen.last().unwrap().next());
        }
        assert_eq!(
            seen,
            [
                SelectLevel::Atom,
                SelectLevel::Residue,
                SelectLevel::Chain,
                SelectLevel::Molecule,
                SelectLevel::Atom
            ]
        );
        assert_eq!(SelectShape::Lasso.next(), SelectShape::Click);
        for level in SelectLevel::ALL {
            assert_eq!(SelectLevel::from_word(level.word()), Some(level));
        }
    }

    #[test]
    fn the_tool_button_toggles_between_click_and_the_last_shape() {
        let mut tool = SelectTool::default();
        assert_eq!(tool.toggled_shape(), SelectShape::Box, "first use: box");
        tool.set("shape lasso").unwrap();
        assert_eq!(tool.toggled_shape(), SelectShape::Click);
        tool.set("shape toggle").unwrap();
        assert!(!tool.draws());
        tool.set("shape toggle").unwrap();
        assert_eq!(tool.shape, SelectShape::Lasso, "remembers the last shape");
        tool.set("shape next").unwrap();
        assert_eq!(tool.shape, SelectShape::Click);
        tool.set("level next").unwrap();
        assert_eq!(tool.level, SelectLevel::Residue);
    }

    #[test]
    fn a_corner_press_or_a_long_press_opens_the_flyout() {
        assert_eq!(classify_press(0.0, false), ToolPress::Activate);
        assert_eq!(classify_press(HOLD_SECS - 0.01, false), ToolPress::Activate);
        assert_eq!(classify_press(HOLD_SECS, false), ToolPress::OpenFlyout);
        assert_eq!(classify_press(0.0, true), ToolPress::OpenFlyout);
    }

    #[test]
    fn modifiers_pick_the_combine_mode() {
        let mods = |shift, alt, command| egui::Modifiers {
            shift,
            alt,
            command,
            ..Default::default()
        };
        assert_eq!(
            combine_for(mods(false, false, false), false),
            Combine::Replace
        );
        assert_eq!(combine_for(mods(true, false, false), false), Combine::Add);
        assert_eq!(
            combine_for(mods(false, true, false), false),
            Combine::Subtract
        );
        assert_eq!(
            combine_for(mods(false, false, true), false),
            Combine::Replace
        );
        assert_eq!(
            combine_for(mods(false, false, true), true),
            Combine::Subtract
        );
    }

    #[test]
    fn shift_adds_and_ctrl_subtracts() {
        let current: Mask = std::sync::Arc::new(FixedBitSet::with_capacity_and_blocks(4, [0b0011]));
        let hits = FixedBitSet::with_capacity_and_blocks(4, [0b0110]);
        let run = |how| {
            combine(Some(&current), hits.clone(), how)
                .ones()
                .collect::<Vec<_>>()
        };
        assert_eq!(run(Combine::Add), [0, 1, 2]);
        assert_eq!(run(Combine::Subtract), [0]);
        assert_eq!(run(Combine::Replace), [1, 2]);
    }

    #[test]
    fn the_verb_sets_shape_and_level_and_rejects_the_rest() {
        let mut tool = SelectTool::default();
        assert!(tool.set("shape lasso").is_ok());
        assert!(tool.set("level chain").is_ok());
        assert_eq!(
            (tool.shape, tool.level),
            (SelectShape::Lasso, SelectLevel::Chain)
        );
        assert!(tool.set("shape hexagon").is_err());
        assert!(tool.set("level").is_err());
        assert!(tool.set("shape box extra").is_err());
        assert_eq!((tool.shape_word(), tool.level_word()), ("lasso", "chain"));
        assert!(tool.set("level molecule").is_ok());
        assert_eq!(tool.level, SelectLevel::Molecule);
    }
}
