//! Core molecular data model.
//!
//! Columnar (SoA) atom storage, chain/residue hierarchy as index ranges,
//! bitset selections, and geometric bond perception. No GPU, no I/O.

pub mod altloc;
pub mod analysis;
pub mod annotations;
pub mod backbone;
pub mod bonds;
pub mod builder;
pub mod cartoon;
pub mod coords;
pub mod dssp;
pub mod element;
pub mod feasible_cell;
pub mod gaussian_mesh;
pub mod gaussian_surface;
pub mod glycan;
pub mod hbond;
pub mod intern;
pub mod material;
pub mod residue_class;
pub mod sasa;
pub mod select;
pub mod ses;
pub mod skin_surface;
pub mod spatial;
pub mod structure;
pub mod topology;
pub mod weighted_delaunay;

pub use analysis::distance;
pub use annotations::{AnnotationCategory, Annotations};
pub use backbone::{
    nucleic_ladder, putty_radius, trace as backbone_trace, tube as backbone_tube, Trace, Tube,
    TUBE_RADIUS,
};
pub use bonds::{
    adjacency, bond_strand_endpoints, bond_strands, strand_axis, Adjacency, BondOrder, BondStrand,
    BondTable,
};
pub use builder::{AtomExtra, AtomRow, TopologyBuilder};
pub use cartoon::{build as build_cartoon, CartoonMesh};
pub use coords::CoordSet;
pub use dssp::{assign as assign_dssp, DsspCode};
pub use element::Element;
pub use glycan::{
    build_linkage_mesh as build_glycan_linkage_mesh, build_mesh as build_glycan_mesh, Attachment,
    GlycanFrame, GlycanPlan, GlycanResidue, PolytopeMesh, Shape, SnfgColor,
};
pub use hbond::{hydrogen_bonds_into, HydrogenBond};
pub use intern::{InternId, Interner};
pub use material::MaterialPreset;
pub use residue_class::{ClassCounts, PolymerHint, ResidueClass, Roles};
pub use select::{select, Expr, SelectError};
pub use spatial::Grid;
pub use structure::{FrameSource, Structure, StructureError, DEFAULT_FRAME_BUDGET};
pub use topology::{
    flags, ChainRec, ExplicitBond, ExplicitBondKind, ResidueRec, SecondaryStructure, Topology,
    TopologyError,
};

pub use fixedbitset;
pub use glam;
