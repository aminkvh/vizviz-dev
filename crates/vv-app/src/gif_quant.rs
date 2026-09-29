//! Colour reduction for one GIF frame: a palette of at most 256 colours
//! and the pixels' indices into it.
//!
//! Flat regions must land on an exact palette entry in every frame, or
//! per-frame palettes make a plain background shimmer. So the palette
//! keeps the frame's dominant exact colours verbatim, fills the rest with
//! NeuQuant, and pixels that already match an entry skip dithering.

use std::collections::HashMap;

use color_quant::NeuQuant;

/// A colour covering this share of the frame is kept exactly.
const PIN_SHARE: f64 = 0.005;
const MAX_PINNED: usize = 32;
/// NeuQuant trains on about this many pixels; bigger samples are thinned.
const TRAIN_PIXELS: usize = 60_000;

pub type Rgb = [u8; 3];

fn key(c: Rgb) -> u32 {
    u32::from(c[0]) << 16 | u32::from(c[1]) << 8 | u32::from(c[2])
}

/// The colours of the frame's pixels where `active`, most frequent first.
fn histogram(rgb: &[u8], active: &[bool]) -> Vec<(Rgb, u32)> {
    let mut counts: HashMap<u32, u32> = HashMap::new();
    for (px, _) in rgb.chunks_exact(3).zip(active).filter(|(_, a)| **a) {
        *counts.entry(key([px[0], px[1], px[2]])).or_default() += 1;
    }
    let mut sorted: Vec<(Rgb, u32)> = counts
        .into_iter()
        .map(|(k, n)| ([(k >> 16) as u8, (k >> 8) as u8, k as u8], n))
        .collect();
    sorted.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    sorted
}

/// At most `slots` colours standing for the `active` pixels of `rgb`.
pub fn build_palette(rgb: &[u8], active: &[bool], slots: usize) -> Vec<Rgb> {
    let colours = histogram(rgb, active);
    if colours.len() <= slots {
        return colours.into_iter().map(|(c, _)| c).collect();
    }
    let total: u32 = colours.iter().map(|(_, n)| n).sum();
    let pinned: Vec<Rgb> = colours
        .iter()
        .take(MAX_PINNED.min(slots / 2))
        .take_while(|(_, n)| f64::from(*n) >= PIN_SHARE * f64::from(total))
        .map(|(c, _)| *c)
        .collect();
    let mut palette = pinned.clone();
    palette.extend(train(&colours, &pinned, slots - pinned.len()));
    palette
}

/// NeuQuant colours for every colour outside `pinned`, weighted by how
/// often it occurs.
fn train(colours: &[(Rgb, u32)], pinned: &[Rgb], slots: usize) -> Vec<Rgb> {
    let mut sample = Vec::new();
    for (c, n) in colours.iter().filter(|(c, _)| !pinned.contains(c)) {
        for _ in 0..*n.min(&64) {
            sample.extend_from_slice(&[c[0], c[1], c[2], 255]);
        }
    }
    let step = (sample.len() / 4 / TRAIN_PIXELS).clamp(1, 30) as i32;
    let map = NeuQuant::new(step, slots, &sample).color_map_rgb();
    map.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect()
}

/// Nearest-colour lookup over a palette: exact matches first, then a
/// cache over 6-bit colour cells.
struct Nearest<'a> {
    colours: &'a [Rgb],
    exact: HashMap<u32, u8>,
    cells: Vec<Option<u8>>,
}

impl<'a> Nearest<'a> {
    fn new(colours: &'a [Rgb]) -> Self {
        let exact = colours
            .iter()
            .enumerate()
            .map(|(i, c)| (key(*c), i as u8))
            .rev()
            .collect();
        Self {
            colours,
            exact,
            cells: vec![None; 1 << 18],
        }
    }

    fn exact(&self, c: Rgb) -> Option<u8> {
        self.exact.get(&key(c)).copied()
    }

    fn nearest(&mut self, c: Rgb) -> u8 {
        let cell =
            (usize::from(c[0] >> 2) << 12) | (usize::from(c[1] >> 2) << 6) | usize::from(c[2] >> 2);
        if let Some(i) = self.cells[cell] {
            return i;
        }
        let centre = [c[0] | 2, c[1] | 2, c[2] | 2];
        let best = self.closest(centre);
        self.cells[cell] = Some(best);
        best
    }

    fn closest(&self, c: Rgb) -> u8 {
        let dist = |p: &Rgb| -> i32 {
            (0..3)
                .map(|k| (i32::from(p[k]) - i32::from(c[k])).pow(2))
                .sum()
        };
        let (best, _) = self
            .colours
            .iter()
            .enumerate()
            .min_by_key(|(_, p)| dist(p))
            .unwrap_or((0, &[0; 3]));
        best as u8
    }
}

/// What to do with each pixel of a frame.
pub struct Pixels<'a> {
    pub rgb: &'a [u8],
    /// Pixels to draw; the rest become `transparent`.
    pub active: &'a [bool],
    pub width: usize,
}

/// Indices of `pixels` into `palette`; inactive pixels get `transparent`.
/// With `dither`, Floyd-Steinberg error diffusion (7, 3, 5, 1 sixteenths)
/// spreads what the palette could not reproduce.
pub fn map_to_palette(pixels: &Pixels, palette: &[Rgb], transparent: u8, dither: bool) -> Vec<u8> {
    let (w, h) = (pixels.width, pixels.active.len() / pixels.width);
    let mut lookup = Nearest::new(palette);
    let mut out = vec![transparent; w * h];
    let mut errors = [vec![[0.0f32; 3]; w + 2], vec![[0.0f32; 3]; w + 2]];
    for y in 0..h {
        let (this, next) = (y % 2, (y + 1) % 2);
        errors[next].fill([0.0; 3]);
        for x in 0..w {
            let i = y * w + x;
            let carried = std::mem::take(&mut errors[this][x + 1]);
            if !pixels.active[i] {
                continue;
            }
            let src = [
                pixels.rgb[i * 3],
                pixels.rgb[i * 3 + 1],
                pixels.rgb[i * 3 + 2],
            ];
            if let Some(index) = lookup.exact(src) {
                out[i] = index;
                continue;
            }
            let wanted: [f32; 3] =
                std::array::from_fn(|k| (f32::from(src[k]) + carried[k]).clamp(0.0, 255.0));
            let index = lookup.nearest(wanted.map(|v| v.round() as u8));
            out[i] = index;
            if dither {
                let got = palette[usize::from(index)];
                let miss: [f32; 3] = std::array::from_fn(|k| wanted[k] - f32::from(got[k]));
                spread(&mut errors, this, next, x + 1, miss);
            }
        }
    }
    out
}

/// Adds `miss` to the pixels right of and below column `x` (offset by
/// one so `x - 1` never underflows).
fn spread(errors: &mut [Vec<[f32; 3]>; 2], this: usize, next: usize, x: usize, miss: [f32; 3]) {
    let targets = [
        (this, x + 1, 7.0),
        (next, x - 1, 3.0),
        (next, x, 5.0),
        (next, x + 1, 1.0),
    ];
    for (row, col, share) in targets {
        for k in 0..3 {
            errors[row][col][k] += miss[k] * share / 16.0;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn frame(w: usize, h: usize, f: impl Fn(usize, usize) -> Rgb) -> Vec<u8> {
        let f = &f;
        (0..h)
            .flat_map(move |y| (0..w).flat_map(move |x| f(x, y)))
            .collect()
    }

    fn mean_error(rgb: &[u8], palette: &[Rgb], indices: &[u8]) -> f64 {
        let total: f64 = rgb
            .chunks_exact(3)
            .zip(indices)
            .map(|(p, &i)| {
                (0..3)
                    .map(|k| (f64::from(p[k]) - f64::from(palette[usize::from(i)][k])).abs())
                    .sum::<f64>()
                    / 3.0
            })
            .sum();
        total / indices.len() as f64
    }

    #[test]
    fn a_gradient_gets_at_most_256_colours_and_a_small_error() {
        let (w, h) = (128, 64);
        let rgb = frame(w, h, |x, y| {
            [(x * 2) as u8, (y * 4) as u8, ((x + y) * 2) as u8]
        });
        let active = vec![true; w * h];
        let palette = build_palette(&rgb, &active, 256);
        assert!(palette.len() <= 256);
        let pixels = Pixels {
            rgb: &rgb,
            active: &active,
            width: w,
        };
        let indices = map_to_palette(&pixels, &palette, 0, false);
        let error = mean_error(&rgb, &palette, &indices);
        assert!(error < 6.0, "mean error {error}");
    }

    #[test]
    fn dithering_keeps_the_average_colour_of_a_smooth_ramp() {
        let (w, h) = (256, 8);
        let rgb = frame(w, h, |x, _| [x as u8, x as u8, x as u8]);
        let active = vec![true; w * h];
        let palette: Vec<Rgb> = (0..16).map(|i| [i * 17; 3]).collect();
        let pixels = Pixels {
            rgb: &rgb,
            active: &active,
            width: w,
        };
        let indices = map_to_palette(&pixels, &palette, 0, true);
        let mean = |f: &dyn Fn(usize) -> f64| (0..w * h).map(f).sum::<f64>() / (w * h) as f64;
        let source = mean(&|i| f64::from(rgb[i * 3]));
        let shown = mean(&|i| f64::from(palette[usize::from(indices[i])][0]));
        assert!((source - shown).abs() < 1.0, "{source} vs {shown}");
    }

    #[test]
    fn a_flat_background_stays_exact_next_to_shading() {
        let (w, h) = (96, 96);
        let background = [243, 245, 249];
        let rgb = frame(w, h, |x, y| {
            let inside = (x as i32 - 48).pow(2) + (y as i32 - 48).pow(2) < 30 * 30;
            if inside {
                [(x * 2) as u8, (y * 2) as u8, (x + y) as u8]
            } else {
                background
            }
        });
        let active = vec![true; w * h];
        let palette = build_palette(&rgb, &active, 255);
        let pixels = Pixels {
            rgb: &rgb,
            active: &active,
            width: w,
        };
        let indices = map_to_palette(&pixels, &palette, 255, true);
        let corner = palette[usize::from(indices[0])];
        assert_eq!(corner, background);
        for i in (0..w * h).filter(|i| rgb[i * 3..i * 3 + 3] == background) {
            assert_eq!(palette[usize::from(indices[i])], background, "pixel {i}");
        }
    }

    #[test]
    fn inactive_pixels_take_the_transparent_index() {
        let rgb = frame(4, 1, |x, _| [x as u8 * 60, 0, 0]);
        let active = [true, false, true, false];
        let palette = build_palette(&rgb, &active, 255);
        let pixels = Pixels {
            rgb: &rgb,
            active: &active,
            width: 4,
        };
        let indices = map_to_palette(&pixels, &palette, 200, true);
        assert_eq!((indices[1], indices[3]), (200, 200));
        assert_eq!(palette.len(), 2);
    }
}
