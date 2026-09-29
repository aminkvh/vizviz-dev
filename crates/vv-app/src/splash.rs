//! The logo, for the title bar and for the splash: painted from the CPU
//! the moment a window opens, so the second the GPU takes to start
//! (mostly the Vulkan loader) shows vizviz instead of a blank window.
//!
//! Two windows exist before `State` is up: `Splash`, a small borderless
//! window shown immediately, and `Opening`, the real
//! (maximized or fixed-size) window, created hidden so it only appears
//! once `State::new` has it ready to draw.

use std::cell::RefCell;
use std::num::NonZeroU32;
use std::sync::Arc;

use winit::dpi::{LogicalSize, PhysicalPosition};
use winit::event_loop::ActiveEventLoop;
use winit::window::{Icon, Window, WindowId, WindowLevel};

use crate::layout::SavedLayout;

const SPLASH_SIZE: LogicalSize<f64> = LogicalSize::new(480.0, 360.0);
/// Logical px (DPI-scaled in `paint`) for the status and credit lines.
const STATUS_SIZE: f32 = 13.0;
const CREDIT_SIZE: f32 = 12.0;

/// The small borderless window shown while the GPU comes up and the real
/// window (`Opening`, hidden) is built.
pub struct Splash {
    pub window: Arc<Window>,
    theme: crate::theme::ThemeMode,
    /// What's loading right now, shown under the logo; empty until
    /// `State::new`'s progress callback (`set_status`) names the first
    /// stage.
    status: RefCell<String>,
}

impl Splash {
    pub fn new(event_loop: &ActiveEventLoop, theme: crate::theme::ThemeMode) -> Self {
        let mut attributes = Window::default_attributes()
            .with_title("vizviz")
            .with_window_icon(icon())
            .with_inner_size(SPLASH_SIZE)
            .with_resizable(false)
            .with_decorations(false)
            .with_window_level(WindowLevel::AlwaysOnTop)
            .with_position(centered_position(event_loop, SPLASH_SIZE));
        #[cfg(windows)]
        {
            use winit::platform::windows::WindowAttributesExtWindows;
            attributes = attributes.with_skip_taskbar(true);
        }
        Self {
            window: Arc::new(
                event_loop
                    .create_window(attributes)
                    .expect("create splash window"),
            ),
            theme,
            status: RefCell::new(String::new()),
        }
    }

    pub fn id(&self) -> WindowId {
        self.window.id()
    }

    pub fn paint(&self) {
        paint(&self.window, self.theme, &self.status.borrow());
    }

    /// Updates the status line and repaints immediately: `State::new`'s
    /// progress callback, called while it still blocks the event loop
    /// (softbuffer's `present` is synchronous, so this is safe from
    /// inside it).
    pub fn set_status(&self, status: &str) {
        *self.status.borrow_mut() = status.to_string();
        self.paint();
    }
}

/// The primary monitor's centre for a `size`-logical window, in physical
/// screen coordinates (falls back to the first available monitor, and to
/// the platform default placement if neither is reported).
fn centered_position(
    event_loop: &ActiveEventLoop,
    size: LogicalSize<f64>,
) -> PhysicalPosition<f64> {
    let monitor = event_loop
        .primary_monitor()
        .or_else(|| event_loop.available_monitors().next());
    let Some(monitor) = monitor else {
        return PhysicalPosition::new(0.0, 0.0);
    };
    let physical = size.to_physical::<f64>(monitor.scale_factor());
    let origin = monitor.position();
    let extent = monitor.size();
    PhysicalPosition::new(
        origin.x as f64 + (extent.width as f64 - physical.width) / 2.0,
        origin.y as f64 + (extent.height as f64 - physical.height) / 2.0,
    )
}

/// The real window: created hidden so it appears only once `State::new`
/// has brought up the GPU and is ready to draw the first frame.
pub struct Opening {
    pub window: Arc<Window>,
    pub layout: SavedLayout,
}

impl Opening {
    /// 1400x900 keeps a script's or live session's screenshots
    /// reproducible; anyone else's window is maximized once shown.
    pub fn new(event_loop: &ActiveEventLoop, layout: SavedLayout) -> Self {
        let attributes = Window::default_attributes()
            .with_title("vizviz")
            .with_window_icon(icon())
            .with_inner_size(LogicalSize::new(1400.0, 900.0))
            .with_min_inner_size(LogicalSize::new(1024.0, 640.0))
            .with_visible(false);
        Self {
            window: Arc::new(event_loop.create_window(attributes).expect("create window")),
            layout,
        }
    }
}

/// Decodes a straight-alpha RGBA8 PNG embedded at compile time: (pixels,
/// width, height). `None` for any other PNG color type/bit depth.
fn decode_png(bytes: &[u8]) -> Option<(Vec<u8>, u32, u32)> {
    let decoder = png::Decoder::new(bytes);
    let mut reader = decoder.read_info().ok()?;
    let mut rgba = vec![0; reader.output_buffer_size()];
    let info = reader.next_frame(&mut rgba).ok()?;
    if info.color_type != png::ColorType::Rgba || info.bit_depth != png::BitDepth::Eight {
        return None;
    }
    rgba.truncate(info.buffer_size());
    Some((rgba, info.width, info.height))
}

/// The mark and wordmark for `mode`, both transparent outside the
/// artwork with dark or light lettering to suit the theme
/// (`assets/logo-full-*.png`; also used by the viewport's start card,
/// `start.rs`).
pub fn full_logo(mode: crate::theme::ThemeMode) -> Option<(Vec<u8>, u32, u32)> {
    match mode {
        crate::theme::ThemeMode::Dark => decode_png(include_bytes!("../assets/logo-full-dark.png")),
        crate::theme::ThemeMode::Light => {
            decode_png(include_bytes!("../assets/logo-full-light.png"))
        }
    }
}

fn icon() -> Option<Icon> {
    let (rgba, w, h) = decode_png(include_bytes!("../assets/icon.png"))?;
    Icon::from_rgba(rgba, w, h).ok()
}

/// Fills `window` with `theme`'s surface colour, the full logo centred
/// above `status` (may be empty, before the first stage is known), and a
/// credit line at the bottom. Best effort: a platform without a CPU
/// surface just skips the splash.
fn paint(window: &Arc<Window>, theme: crate::theme::ThemeMode, status: &str) {
    let size = window.inner_size();
    let (Some(w), Some(h)) = (NonZeroU32::new(size.width), NonZeroU32::new(size.height)) else {
        return;
    };
    let Ok(context) = softbuffer::Context::new(window.clone()) else {
        return;
    };
    let Ok(mut surface) = softbuffer::Surface::new(&context, window.clone()) else {
        return;
    };
    if surface.resize(w, h).is_err() {
        return;
    }
    let Ok(mut buffer) = surface.buffer_mut() else {
        return;
    };
    let tokens = crate::theme::Tokens::of(theme);
    let [br, bg, bb, _] = tokens.surface.to_array().map(u32::from);
    buffer.fill(br << 16 | bg << 8 | bb);

    let scale = window.scale_factor() as f32;
    let [tr, tg, tb, _] = tokens.text_muted.to_array();
    let text_color = [tr, tg, tb];
    let logo = full_logo(theme);
    let logo_h = logo.as_ref().map_or(0, |&(_, _, lh)| lh) as f32;
    let status_px = STATUS_SIZE * scale;
    let gap = crate::theme::space::PAD * scale;
    let status_h = if status.is_empty() {
        0.0
    } else {
        gap + crate::text::line_height(status_px)
    };
    let y0 = ((size.height as f32 - logo_h - status_h) / 2.0).max(0.0) as u32;

    if let Some((rgba, lw, lh)) = &logo {
        let (lw, lh) = (*lw, *lh);
        let x0 = size.width.saturating_sub(lw) / 2;
        for y in 0..lh.min(size.height.saturating_sub(y0)) {
            for x in 0..lw.min(size.width) {
                let p = &rgba[((y * lw + x) * 4) as usize..][..4];
                let a = u32::from(p[3]);
                let mix = |fg: u8, bg: u32| (u32::from(fg) * a + bg * (255 - a)) / 255;
                buffer[((y0 + y) * size.width + x0 + x) as usize] =
                    mix(p[0], br) << 16 | mix(p[1], bg) << 8 | mix(p[2], bb);
            }
        }
    }

    if !status.is_empty() {
        let y = y0 as f32 + logo_h + gap;
        draw_centered(
            &mut buffer,
            size.width,
            size.height,
            y,
            status,
            status_px,
            text_color,
        );
    }

    let credit = format!(
        "vizviz {} \u{b7} \u{a9} 2026 {}",
        env!("CARGO_PKG_VERSION"),
        env!("CARGO_PKG_AUTHORS"),
    );
    let credit_px = CREDIT_SIZE * scale;
    let margin = crate::theme::space::WIDE * scale;
    let credit_y = size.height as f32 - margin - crate::text::line_height(credit_px);
    draw_centered(
        &mut buffer,
        size.width,
        size.height,
        credit_y,
        &credit,
        credit_px,
        text_color,
    );

    let _ = buffer.present();
}

/// Draws `text` horizontally centred with its top at `y`, alpha-blending
/// each glyph over the packed-0RGB `buffer` the same way the logo above
/// is blended.
fn draw_centered(
    buffer: &mut [u32],
    width: u32,
    height: u32,
    y: f32,
    text: &str,
    size_px: f32,
    color: [u8; 3],
) {
    let x = (width as f32 - crate::text::text_width(text, size_px)) / 2.0;
    crate::text::walk_glyphs(x, y, text, size_px, |px, py, coverage| {
        if px < 0 || py < 0 || px as u32 >= width || py as u32 >= height {
            return;
        }
        let idx = (py as u32 * width + px as u32) as usize;
        let bg = buffer[idx];
        let mix = |fg: u8, shift: u32| {
            let bg_c = ((bg >> shift) & 0xFF) as f32;
            (f32::from(fg) * coverage + bg_c * (1.0 - coverage)).round() as u32
        };
        buffer[idx] = mix(color[0], 16) << 16 | mix(color[1], 8) << 8 | mix(color[2], 0);
    });
}
