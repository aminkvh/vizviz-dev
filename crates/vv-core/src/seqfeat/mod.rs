//! Per-residue sequence features: the computations behind the sequence
//! strip's annotation tracks. Each kernel returns plain per-residue data
//! (indices into `Topology::residues`), so a viewer, a script or the
//! Python bindings can share them.

pub mod align;
pub mod conservation;
pub mod contacts;
mod letters;
pub mod links;
pub mod missing;
pub mod motifs;
pub mod props;
pub mod sasa;

pub use conservation::{conservation, Column, Conservation};
pub use contacts::{residue_contacts, ResidueContact};
pub use letters::one_letter;
pub use links::{disulfides, glycosylated, Glycosylation};
pub use missing::{unobserved, unobserved_in_entity, EntityChain, Gap, Unobserved};
pub use motifs::{liabilities, sequons, Hit, Liability, Sequon, SequonKind};
pub use props::{properties, Props};
pub use sasa::relative_sasa;
