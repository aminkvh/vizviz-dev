//! Per-residue sequence features: the computations behind the sequence
//! strip's annotation tracks. Each kernel returns plain per-residue data
//! (indices into `Topology::residues`), so a viewer, a script or the
//! Python bindings can share them.

pub mod contacts;
pub mod links;
pub mod missing;
pub mod motifs;
pub mod sasa;

pub use contacts::{residue_contacts, ResidueContact};
pub use links::{disulfides, glycosylated, Glycosylation};
pub use missing::{unobserved, Gap, Unobserved};
pub use motifs::{liabilities, sequons, Hit, Liability, Sequon, SequonKind};
pub use sasa::relative_sasa;
