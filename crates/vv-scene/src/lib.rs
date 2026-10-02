//! Document model and command bus.
//!
//! Everything the UI, CLI, and (later) Python or MCP can do to a loaded
//! set of structures is a [`Command`] applied to a [`Scene`]; undo is an
//! inverse command, never a snapshot (see [`CommandHistory`]). This crate
//! has no GPU dependency — `vv-app` is the only place that talks to both
//! this and `vv-render`.

/// The command language's version (semver, docs/VERSIONING.md): the
/// commands, their arguments and what they print. A client of a running
/// app (`vizviz --listen`) refuses a different major version.
pub const LANGUAGE_VERSION: &str = "2.1.0";

pub mod command;
pub mod history;
pub mod scene;
pub mod script;
mod select_ops;
pub mod selection;
pub mod session;
pub mod slotmap;
pub mod values;

pub use command::{Command, SceneError};
pub use history::CommandHistory;
pub use scene::{
    ActiveSelection, Caption, ColorScheme, LoadedStructure, Material, Measurement, Rep, RepId,
    Representation, Scene, StructureId,
};
pub use script::{
    parse_representation, run_line, run_script, split_script, ScriptError, Spec as CommandSpec,
    SPECS as COMMAND_SPECS,
};
pub use selection::{empty_mask, mask_of_one, Mask, SelectionSet};
pub use session::{SessionError, SessionFile};
pub use slotmap::Id;
pub use values::{ValueChannel, ValuesError};
