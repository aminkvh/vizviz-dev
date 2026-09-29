//! The `_atom_site` column-to-atom mapping shared by the mmCIF text reader
//! and the BinaryCIF reader: each supplies one row's values through
//! [`Fields`], [`AtomSink`] turns them into builder rows, and [`assemble`]
//! stitches the per-chunk builders back together.

use rayon::prelude::*;
use vv_core::glam::Vec3;
use vv_core::{AtomExtra, AtomRow, Element, Structure, TopologyBuilder};

use crate::mmcif_entity::EntityInfo;
use crate::ParseError;

/// The PDBx dictionary has no segment id. This local item carries a
/// PDB/MD segid so it survives a round trip through mmCIF.
pub(crate) const SEGID_ITEM: &str = "vizviz_segid";

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Role {
    Group,
    Id,
    Element,
    Name,
    AltLoc,
    Comp,
    Asym,
    Entity,
    SeqId,
    InsCode,
    X,
    Y,
    Z,
    Occupancy,
    BFactor,
    Charge,
    AuthSeqId,
    AuthAsym,
    AuthAtomName,
    Model,
    Segid,
}

impl Role {
    pub(crate) fn from_item(item: &str) -> Option<Role> {
        Some(match item {
            "group_PDB" => Role::Group,
            "id" => Role::Id,
            "type_symbol" => Role::Element,
            "label_atom_id" => Role::Name,
            "label_alt_id" => Role::AltLoc,
            "label_comp_id" => Role::Comp,
            "label_asym_id" => Role::Asym,
            "label_entity_id" => Role::Entity,
            "label_seq_id" => Role::SeqId,
            "pdbx_PDB_ins_code" => Role::InsCode,
            "Cartn_x" => Role::X,
            "Cartn_y" => Role::Y,
            "Cartn_z" => Role::Z,
            "occupancy" => Role::Occupancy,
            "B_iso_or_equiv" => Role::BFactor,
            "pdbx_formal_charge" => Role::Charge,
            "auth_seq_id" => Role::AuthSeqId,
            "auth_asym_id" => Role::AuthAsym,
            "auth_atom_id" => Role::AuthAtomName,
            "pdbx_PDB_model_num" => Role::Model,
            SEGID_ITEM => Role::Segid,
            _ => return None,
        })
    }
}

/// Column index of each role, or `usize::MAX` when absent.
pub(crate) struct Columns {
    pub(crate) index: [usize; 21],
    pub(crate) count: usize,
}

impl Columns {
    pub(crate) fn new(items: &[&str]) -> Result<Self, ParseError> {
        let mut index = [usize::MAX; 21];
        for (i, item) in items.iter().enumerate() {
            if let Some(role) = Role::from_item(item) {
                index[role as usize] = i;
            }
        }
        let cols = Self {
            index,
            count: items.len(),
        };
        for (role, name) in [
            (Role::X, "Cartn_x"),
            (Role::Y, "Cartn_y"),
            (Role::Z, "Cartn_z"),
            (Role::Comp, "label_comp_id"),
            (Role::Asym, "label_asym_id"),
        ] {
            if !cols.has(role) {
                return Err(ParseError::MissingColumn(name));
            }
        }
        if !cols.has(Role::Name) && !cols.has(Role::AuthAtomName) {
            return Err(ParseError::MissingColumn("label_atom_id"));
        }
        Ok(cols)
    }

    pub(crate) fn has(&self, role: Role) -> bool {
        self.index[role as usize] != usize::MAX
    }

    /// The column of `role`, when the table has one.
    pub(crate) fn of(&self, role: Role) -> Option<usize> {
        Some(self.index[role as usize]).filter(|&i| i != usize::MAX)
    }
}

/// One `_atom_site` row's values by role; `None` for an absent column or a
/// null (`.`/`?`) value.
pub(crate) trait Fields<'a> {
    fn text(&self, role: Role) -> Option<&'a [u8]>;
    fn int(&self, role: Role) -> Option<i32>;
    fn float(&self, role: Role) -> Option<f32>;
}

fn atom_name(bytes: &[u8]) -> [u8; 4] {
    let mut name = [b' '; 4];
    let n = bytes.len().min(4);
    name[..n].copy_from_slice(&bytes[..n]);
    name
}

fn as_str(v: Option<&[u8]>) -> &str {
    v.and_then(|b| std::str::from_utf8(b).ok()).unwrap_or("")
}

fn parse_charge(v: Option<i32>) -> i8 {
    v.map_or(0, |c| c.clamp(-9, 9) as i8)
}

pub(crate) struct ChunkResult {
    pub(crate) builder: TopologyBuilder,
    /// Coordinates of rows belonging to other models, keyed by model number.
    other_models: Vec<(i32, Vec<Vec3>)>,
    pub(crate) error: Option<ParseError>,
}

/// Rows of one chunk pushed into a private builder.
pub(crate) struct AtomSink<'a> {
    builder: TopologyBuilder,
    other_models: Vec<(i32, Vec<Vec3>)>,
    first_model: i32,
    last_segid: &'a str,
}

impl<'a> AtomSink<'a> {
    pub(crate) fn new(estimate: usize, first_model: i32) -> Self {
        Self {
            builder: TopologyBuilder::with_capacity(estimate),
            other_models: Vec::new(),
            first_model,
            last_segid: "",
        }
    }

    pub(crate) fn push<F: Fields<'a>>(&mut self, f: &F) -> Result<(), &'static str> {
        let (Some(x), Some(y), Some(z)) = (f.float(Role::X), f.float(Role::Y), f.float(Role::Z))
        else {
            return Err("unreadable coordinates");
        };
        let position = Vec3::new(x, y, z);
        let model = f.int(Role::Model).unwrap_or(self.first_model);
        if model != self.first_model {
            self.push_other_model(model, position);
            return Ok(());
        }
        let name_bytes = f
            .text(Role::Name)
            .or_else(|| f.text(Role::AuthAtomName))
            .unwrap_or(b"");
        let symbol = f.text(Role::Element);
        let element = symbol
            .map(Element::from_symbol)
            .filter(|e| !e.is_unknown())
            .unwrap_or_else(|| Element::from_atom_name(name_bytes));
        let segid = as_str(f.text(Role::Segid));
        if segid != self.last_segid {
            self.builder.end_chain();
            self.last_segid = segid;
        }
        let extra = AtomExtra {
            long_name: std::str::from_utf8(name_bytes).ok().filter(|n| n.len() > 4),
            segid,
            deuterium: matches!(symbol, Some(b"D" | b"d")),
        };
        let asym = as_str(f.text(Role::Asym));
        let auth_seq_id = f.int(Role::AuthSeqId);
        let seq_id = f.int(Role::SeqId);
        let row = AtomRow {
            element,
            name: atom_name(name_bytes),
            serial: f.int(Role::Id).unwrap_or(0).max(0) as u32,
            alt_loc: f
                .text(Role::AltLoc)
                .and_then(|a| a.first().copied())
                .unwrap_or(0),
            comp: as_str(f.text(Role::Comp)),
            asym,
            auth_asym: f.text(Role::AuthAsym).map_or(asym, |a| as_str(Some(a))),
            seq_id: seq_id.or(auth_seq_id).unwrap_or(0),
            auth_seq_id: auth_seq_id.or(seq_id).unwrap_or(0),
            ins_code: f
                .text(Role::InsCode)
                .and_then(|c| c.first().copied())
                .unwrap_or(0),
            entity: f.int(Role::Entity).unwrap_or(0) as u16,
            position,
            occupancy: f.float(Role::Occupancy).unwrap_or(1.0),
            b_factor: f.float(Role::BFactor).unwrap_or(0.0),
            charge: parse_charge(f.int(Role::Charge)),
            hetero: f.text(Role::Group).is_some_and(|g| g == b"HETATM"),
        };
        self.builder.push_with(&row, &extra);
        Ok(())
    }

    fn push_other_model(&mut self, model: i32, position: Vec3) {
        match self.other_models.iter_mut().find(|(m, _)| *m == model) {
            Some((_, v)) => v.push(position),
            None => self.other_models.push((model, vec![position])),
        }
    }

    pub(crate) fn finish(self, error: Option<ParseError>) -> ChunkResult {
        ChunkResult {
            builder: self.builder,
            other_models: self.other_models,
            error,
        }
    }
}

/// Joins the chunk results in file order into one structure; models other
/// than the first become extra frames when their atom count matches.
pub(crate) fn assemble(
    results: Vec<ChunkResult>,
    entity: &EntityInfo,
) -> Result<Structure, ParseError> {
    let mut builder = TopologyBuilder::new();
    let mut models: Vec<(i32, Vec<Vec3>)> = Vec::new();
    for result in results {
        if let Some(e) = result.error {
            return Err(e);
        }
        builder.append(result.builder);
        for (model, mut coords) in result.other_models {
            match models.iter_mut().find(|(m, _)| *m == model) {
                Some((_, v)) => v.append(&mut coords),
                None => models.push((model, coords)),
            }
        }
    }
    if builder.atom_count() == 0 {
        return Err(ParseError::NoAtoms);
    }
    let atom_count = builder.atom_count();
    models.sort_by_key(|(m, _)| *m);
    let extra: Vec<Vec<Vec3>> = models
        .into_iter()
        .filter(|(_, v)| v.len() == atom_count)
        .map(|(_, v)| v)
        .collect();
    builder.topology.polymer_hint = entity.hints(&builder.topology);
    builder.finish_with_frames(extra).map_err(ParseError::from)
}

/// Splits `0..rows` into chunks sized for the thread pool.
pub(crate) fn row_chunks(rows: usize) -> Vec<std::ops::Range<usize>> {
    let threads = rayon::current_num_threads().max(1);
    let size = (rows / (threads * 4)).max(16_384);
    (0..rows)
        .step_by(size)
        .map(|s| s..(s + size).min(rows))
        .collect()
}

/// Runs `f` over each chunk in parallel, keeping chunk order.
pub(crate) fn par_chunks<F>(rows: usize, f: F) -> Vec<ChunkResult>
where
    F: Fn(std::ops::Range<usize>) -> ChunkResult + Sync + Send,
{
    row_chunks(rows).into_par_iter().map(f).collect()
}
