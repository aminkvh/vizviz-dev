//! Headless rendering for Python through `vv-cpu`.
//!
//! The CPU backend rather than the GPU one on purpose: it needs no
//! adapter, no window, and no driver, so `s.render()` works identically on
//! a laptop, an HPC login node, and a CI runner, and its output agrees
//! with the GPU renderer where the two overlap (spacefill; see
//! `crates/vv-cpu/tests/agreement.rs`). A GPU-backed option can be added
//! behind the same signature later.

use std::io::BufWriter;
use std::path::PathBuf;

use numpy::{
    PyArray1, PyArray3, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2, PyReadonlyArray3,
    PyUntypedArrayMethods,
};
use pyo3::exceptions::{PyOSError, PyValueError};
use pyo3::prelude::*;
use vv_render::{Camera, ColorScheme, StylePreset};

/// Extra room around the bounding sphere of atom *centers*, so the atoms
/// on the edge (drawn at their van der Waals radius) aren't clipped.
const FRAMING_PADDING: f32 = 2.0;

pub struct RenderOptions<'a> {
    pub width: u32,
    pub height: u32,
    pub style: &'a str,
    /// Radians, applied on top of the auto-framing camera.
    pub yaw: f32,
    pub pitch: f32,
    /// Multiplies the camera distance: `> 1` moves away, `< 1` closer.
    pub zoom: f32,
    /// RGBA override; `None` uses the style preset's own background.
    pub background: Option<(u8, u8, u8, u8)>,
}

pub fn parse_style(name: &str) -> PyResult<StylePreset> {
    Ok(match name {
        "dark_presentation" | "dark" => StylePreset::DarkPresentation,
        "publication_white" | "white" => StylePreset::PublicationWhite,
        "glossy" => StylePreset::Glossy,
        "flat_cel" | "flat" => StylePreset::FlatCel,
        other => {
            return Err(PyValueError::new_err(format!(
                "unknown style {other:?}; expected one of \"dark_presentation\", \
                 \"publication_white\", \"glossy\", \"flat_cel\""
            )))
        }
    })
}

pub fn parse_coloring(name: &str) -> PyResult<ColorScheme> {
    match name {
        "element" => Ok(ColorScheme::Element),
        "chain" => Ok(ColorScheme::Chain),
        "b_factor" => Ok(ColorScheme::BFactor),
        "structure" => Ok(ColorScheme::SecondaryStructure),
        "restype" => Ok(ColorScheme::ResidueType),
        "rainbow" => Ok(ColorScheme::Rainbow),
        "hetero" => Ok(ColorScheme::Hetero),
        other => Err(PyValueError::new_err(format!(
            "unknown coloring {other:?}; expected \"element\", \"chain\", \"structure\", \
             \"restype\", \"rainbow\", \"hetero\", or \"b_factor\""
        ))),
    }
}

/// One packed color per atom for a scene structure: its scheme, or its
/// value channel on the shared scalar ramp with the channel's
/// whole-trajectory range (same rule as the app's `gpu_cache::colors_of`).
pub fn colors_of(
    loaded: &vv_scene::LoadedStructure,
    coloring: &vv_scene::ColorScheme,
    frame: usize,
) -> Vec<u32> {
    let topology = &loaded.structure.topology;
    match coloring {
        vv_scene::ColorScheme::Element => vv_render::colors_for(ColorScheme::Element, topology),
        vv_scene::ColorScheme::Chain => vv_render::colors_for(ColorScheme::Chain, topology),
        vv_scene::ColorScheme::BFactor => vv_render::colors_for(ColorScheme::BFactor, topology),
        vv_scene::ColorScheme::SecondaryStructure => {
            let coords = loaded.structure.frame(frame);
            let positions = coords.positions();
            let single = loaded.structure.frame_count() == 1;
            let codes = vv_core::cartoon::secondary_structure(topology, positions, single);
            vv_render::colors_for_ss(topology, &codes)
        }
        vv_scene::ColorScheme::ResidueType => {
            vv_render::colors_for(ColorScheme::ResidueType, topology)
        }
        vv_scene::ColorScheme::Rainbow => vv_render::colors_for(ColorScheme::Rainbow, topology),
        vv_scene::ColorScheme::Hetero => vv_render::colors_for(ColorScheme::Hetero, topology),
        vv_scene::ColorScheme::ResidueName => {
            vv_render::colors_for(ColorScheme::ResidueName, topology)
        }
        vv_scene::ColorScheme::Occupancy => vv_render::colors_for(ColorScheme::Occupancy, topology),
        vv_scene::ColorScheme::Hydrophobicity => {
            vv_render::colors_for(ColorScheme::Hydrophobicity, topology)
        }
        vv_scene::ColorScheme::WimleyWhite => {
            vv_render::colors_for(ColorScheme::WimleyWhite, topology)
        }
        vv_scene::ColorScheme::SegmentName => {
            vv_render::colors_for(ColorScheme::SegmentName, topology)
        }
        vv_scene::ColorScheme::Zappo => vv_render::colors_for(ColorScheme::Zappo, topology),
        vv_scene::ColorScheme::Taylor => vv_render::colors_for(ColorScheme::Taylor, topology),
        vv_scene::ColorScheme::Clustal => vv_render::colors_for(ColorScheme::Clustal, topology),
        vv_scene::ColorScheme::HelixPropensity => {
            vv_render::colors_for(ColorScheme::HelixPropensity, topology)
        }
        vv_scene::ColorScheme::StrandPropensity => {
            vv_render::colors_for(ColorScheme::StrandPropensity, topology)
        }
        vv_scene::ColorScheme::TurnPropensity => {
            vv_render::colors_for(ColorScheme::TurnPropensity, topology)
        }
        vv_scene::ColorScheme::BuriedIndex => {
            vv_render::colors_for(ColorScheme::BuriedIndex, topology)
        }
        vv_scene::ColorScheme::Nucleotide => {
            vv_render::colors_for(ColorScheme::Nucleotide, topology)
        }
        vv_scene::ColorScheme::MoleculeClass => vv_render::colors_for(ColorScheme::Class, topology),
        vv_scene::ColorScheme::PurinePyrimidine => {
            vv_render::colors_for(ColorScheme::PurinePyrimidine, topology)
        }
        // No bond graph at hand here (unlike the app's `gpu_cache`, which
        // has `loaded.bonds`); falls back to element, like an unattached
        // `Values` channel does below.
        vv_scene::ColorScheme::Fragment => vv_render::colors_for(ColorScheme::Element, topology),
        vv_scene::ColorScheme::Constant([r, g, b]) => vv_render::colors_for(
            ColorScheme::Constant(vv_render::color::rgba(*r, *g, *b)),
            topology,
        ),
        vv_scene::ColorScheme::Values(name) => match loaded.values.get(name) {
            Some(channel) => {
                let (lo, hi) = channel.range();
                vv_render::colors_from_scalar_in(channel.frame(frame), lo, hi)
            }
            None => vv_render::colors_for(ColorScheme::Element, topology),
        },
    }
}

/// Colors for `Structure.render(coloring=, values=, colors=)`. `colors`
/// (an `(n, 3)` or `(n, 4)` uint8 array) wins over `values` (one number
/// per atom or per residue, drawn on the scalar ramp over `range` or its
/// own min..max), which wins over the named scheme.
pub fn colors_from_python(
    structure: &vv_core::Structure,
    coloring: &str,
    values: Option<&Bound<'_, PyAny>>,
    colors: Option<&Bound<'_, PyAny>>,
    range: Option<(f32, f32)>,
) -> PyResult<Vec<u32>> {
    let atoms = structure.topology.atom_count();
    if let Some(colors) = colors {
        let rgba = colors.extract::<PyReadonlyArray2<u8>>().map_err(|_| {
            PyValueError::new_err("colors must be an (atoms, 3) or (atoms, 4) uint8 array")
        })?;
        let shape = rgba.shape();
        if shape[0] != atoms || !(shape[1] == 3 || shape[1] == 4) {
            return Err(PyValueError::new_err(format!(
                "colors must be ({atoms}, 3) or ({atoms}, 4), got {shape:?}"
            )));
        }
        let a = rgba.as_array();
        return Ok((0..atoms)
            .map(|i| {
                let alpha = if shape[1] == 4 { a[[i, 3]] } else { 255 };
                vv_render::pack_rgba(a[[i, 0]], a[[i, 1]], a[[i, 2]], alpha)
            })
            .collect());
    }
    if let Some(values) = values {
        let mut column = scalar_column(values)?;
        let residues = structure.topology.residue_count();
        if column.len() == residues && residues != atoms {
            let index = &structure.topology.residue_index;
            column = index.iter().map(|&r| column[r as usize]).collect();
        }
        if column.len() != atoms {
            return Err(PyValueError::new_err(format!(
                "values has {} entries, the structure has {atoms} atoms ({residues} residues)",
                column.len()
            )));
        }
        return Ok(match range {
            Some((lo, hi)) => vv_render::colors_from_scalar_in(&column, lo, hi),
            None => vv_render::colors_from_scalar(&column),
        });
    }
    Ok(vv_render::colors_for(
        parse_coloring(coloring)?,
        &structure.topology,
    ))
}

/// Any 1-D numeric array or sequence as f32.
pub fn scalar_column(values: &Bound<'_, PyAny>) -> PyResult<Vec<f32>> {
    if let Ok(a) = values.extract::<PyReadonlyArray1<f32>>() {
        return Ok(a.as_array().iter().copied().collect());
    }
    if let Ok(a) = values.extract::<PyReadonlyArray1<f64>>() {
        return Ok(a.as_array().iter().map(|&v| v as f32).collect());
    }
    if let Ok(a) = values.extract::<PyReadonlyArray1<i64>>() {
        return Ok(a.as_array().iter().map(|&v| v as f32).collect());
    }
    if let Ok(a) = values.extract::<PyReadonlyArray1<i32>>() {
        return Ok(a.as_array().iter().map(|&v| v as f32).collect());
    }
    values.extract::<Vec<f32>>().map_err(|_| {
        PyValueError::new_err("values must be a 1-D numeric array (one number per atom)")
    })
}

/// The preset's background as RGBA8, the same conversion `vv-cpu`'s
/// agreement test uses.
pub fn background_bytes(style: StylePreset) -> [u8; 4] {
    let c = style.background();
    [
        (c.r * 255.0).round() as u8,
        (c.g * 255.0).round() as u8,
        (c.b * 255.0).round() as u8,
        (c.a * 255.0).round() as u8,
    ]
}

/// `render_to_array` with `representation: Spacefill` (`Structure.render()`
/// has no representation concept — see the module doc); `Session.render()`
/// calls [`render_to_array`] directly with the structure's actual one.
pub fn render_spacefill_to_array<'py>(
    py: Python<'py>,
    structure: &vv_core::Structure,
    colors: Vec<u32>,
    options: &RenderOptions<'_>,
) -> PyResult<Bound<'py, PyArray3<u8>>> {
    render_to_array(
        py,
        structure,
        colors,
        vv_scene::Representation::Spacefill,
        options,
    )
}

/// The one representation this backend draws exactly as the live GPU
/// viewport would (no bonds, no ribbon triangles -- see the module doc's
/// "spacefill only" framing, extended here to the two full-screen
/// surfaces, which have no such gap since neither draws bonds either).
/// `BallAndStick`, `Tube`, and `Cartoon` still fall back to spacefill,
/// with a warning; a CPU path for those needs real cylinder/mesh
/// geometry this backend doesn't have.
pub fn render_to_array<'py>(
    py: Python<'py>,
    structure: &vv_core::Structure,
    colors: Vec<u32>,
    representation: vv_scene::Representation,
    options: &RenderOptions<'_>,
) -> PyResult<Bound<'py, PyArray3<u8>>> {
    if options.width == 0 || options.height == 0 {
        return Err(PyValueError::new_err("width and height must be positive"));
    }
    if !(options.zoom.is_finite() && options.zoom > 0.0) {
        return Err(PyValueError::new_err("zoom must be a positive number"));
    }
    let style = parse_style(options.style)?;
    let (center, radius) = structure
        .frame(0)
        .bounding_sphere()
        .ok_or_else(|| PyValueError::new_err("cannot render a structure with no atoms"))?;
    let mut camera = Camera::framing(center, radius + FRAMING_PADDING);
    camera.orbit(options.yaw, options.pitch);
    camera.zoom(options.zoom);
    let background = options
        .background
        .map(|(r, g, b, a)| [r, g, b, a])
        .unwrap_or_else(|| background_bytes(style));
    let (width, height) = (options.width, options.height);
    let pixels = py.detach(|| match representation {
        vv_scene::Representation::GaussianSurface => {
            vv_cpu::GaussianSurfaceScene::from_structure_colored(
                structure,
                colors,
                vv_core::gaussian_surface::DEFAULT_BLOB_FACTOR,
            )
            .render(&camera, width, height, style, background)
        }
        vv_scene::Representation::SkinSurface => vv_cpu::SkinSurfaceScene::from_structure_colored(
            structure,
            0,
            colors,
            vv_core::skin_surface::DEFAULT_SHRINK,
        )
        .render(&camera, width, height, style, background),
        _ => vv_cpu::render_spacefill_colored(
            structure, colors, &camera, width, height, style, background,
        ),
    });
    PyArray1::from_vec(py, pixels).reshape([height as usize, width as usize, 4])
}

/// Writes an `(height, width, 4)` uint8 RGBA array (what `render()`
/// returns) to `path` as a PNG.
#[pyfunction]
pub fn write_png(path: PathBuf, image: PyReadonlyArray3<u8>) -> PyResult<()> {
    let shape = image.shape();
    let (height, width, channels) = (shape[0], shape[1], shape[2]);
    if channels != 4 || height == 0 || width == 0 {
        return Err(PyValueError::new_err(format!(
            "expected an (height, width, 4) uint8 RGBA image, got shape {shape:?}"
        )));
    }
    // `as_slice` fails on non-contiguous input (a slice or a transpose);
    // iterating handles any layout at the cost of a copy.
    let data: Vec<u8> = match image.as_slice() {
        Ok(slice) => slice.to_vec(),
        Err(_) => image.as_array().iter().copied().collect(),
    };
    let file = std::fs::File::create(&path)?;
    let mut encoder = png::Encoder::new(BufWriter::new(file), width as u32, height as u32);
    encoder.set_color(png::ColorType::Rgba);
    encoder.set_depth(png::BitDepth::Eight);
    let png_err = |e: png::EncodingError| PyOSError::new_err(format!("{}: {e}", path.display()));
    let mut writer = encoder.write_header().map_err(png_err)?;
    writer.write_image_data(&data).map_err(png_err)?;
    writer.finish().map_err(png_err)
}
