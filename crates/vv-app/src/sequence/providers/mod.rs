//! The built-in annotation tracks.

mod analysis;
mod antibody;
#[cfg(test)]
mod antibody_tests;
mod chemistry;
mod online;
mod site;
mod structure;
#[cfg(test)]
mod tests;

pub use analysis::{Burial, Conservation};
pub use antibody::Antibody;
pub use chemistry::{AltLocs, Disulfides, Glycans, Liabilities, Modified};
pub use online::{Uniprot, Variants};
pub use site::{Interface, LigandSite};
pub use structure::{Missing, Numbering, SecondaryStructure};

use vv_render::color::rgba;

pub(super) const fn hex(rgb: u32) -> u32 {
    rgba((rgb >> 16) as u8, (rgb >> 8) as u8, rgb as u8)
}
