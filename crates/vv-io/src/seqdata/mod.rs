//! Sequence-level data the strip reads beside the coordinates: each
//! chain's full sequence from the structure file, and UniProt records
//! fetched over the network.

#[cfg(feature = "fetch")]
pub mod abnum;
pub mod entity;
#[cfg(feature = "fetch")]
pub mod uniprot;
