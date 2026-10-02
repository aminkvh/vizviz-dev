//! vizviz desktop application: entry point (CLI parsing, the winit event
//! loop) and glue to `State` (see state.rs), which is everything that
//! exists once a window and a GPU context do.
//!
//! The 3D scene renders into its own offscreen texture at the *viewport
//! tab's* size, not the window's; egui displays that texture as an
//! `Image`. Resizing the offscreen target is debounced (a few stable
//! frames) so dragging a dock divider doesn't rebuild the depth pyramid
//! every frame — see docs/RENDERING.md and state.rs.

// No console window when launched from Explorer or a shortcut. A process
// started from a shell attaches to that shell's console instead
// (`attach_parent_console`), so `--exec` output and `help` still show.
#![cfg_attr(windows, windows_subsystem = "windows")]

mod app_history;
mod commands;
mod gif_quant;
mod gpu_cache;
mod home;
mod home_interface;
mod keys;
mod layout;
mod live;
mod measure;
mod movie;
mod movie_cmd;
mod movie_export;
mod movie_gif;
mod movie_state;
mod movie_timeline;
mod movie_ui;
mod overlays;
mod prefetch;
mod render;
mod ribbon;
mod select_tool;
mod selection_outline;
mod selection_panel;
mod sequence;
mod splash;
mod start;
mod state;
mod studio;
mod text;
mod theme;
mod ui;
mod viewport;
mod widgets;

use std::path::PathBuf;
use std::time::Instant;

use state::State;
use winit::application::ApplicationHandler;
use winit::event::WindowEvent;
use winit::event_loop::{ActiveEventLoop, ControlFlow, EventLoop};
use winit::window::WindowId;

const USAGE: &str = "\
usage: vizviz [FILE] [options]
  FILE               structure to open on startup (.cif/.mmcif/.pdb, optionally .gz)
  --adapter NAME     use the first GPU whose name contains NAME
  --uncapped         render as fast as the GPU allows, no display-refresh cap
                     (default: capped — uncapped burns power/heat for no
                     visible benefit above the monitor's refresh rate)
  --exec SCRIPT      run commands once the window is up; repeatable, `;` separates
                     commands, e.g. --exec \"load 4HHB.cif; style white; screenshot fig.png; quit\"
  --exit-after SECS  quit automatically after SECS seconds (smoke tests)
  --listen [PORT]    let other programs drive this window (Python's
                     vizviz.connect(), the MCP server): 127.0.0.1 only, token
                     in a per-user file; PORT 0 or none lets the OS pick
Commands: `vizviz --exec help` lists them (also in docs/COMMANDS.md). Once running,
use File > Open, Ctrl+P for the command palette, or the Log tab's console line.";

pub(crate) struct Options {
    pub(crate) file: Option<PathBuf>,
    pub(crate) adapter: Option<String>,
    pub(crate) vsync: bool,
    pub(crate) exit_after: Option<f32>,
    /// Command lines to run after startup, already split and validated.
    pub(crate) exec: Vec<String>,
    /// `--listen`: the port (0 = any) to take requests on.
    pub(crate) listen: Option<u16>,
}

impl Options {
    /// Driven by a script or a live client rather than a person: a fixed
    /// window size, and the person's saved layout is left alone.
    pub(crate) fn scripted(&self) -> bool {
        !self.exec.is_empty() || self.listen.is_some()
    }
}

/// Saves the layout on the way out, unless a script or client drove the
/// run (its panel and theme commands are not the person's choices).
fn save_on_exit(options: &Options, state: &State) {
    if !options.scripted() {
        layout::save_layout(&state.saved_layout());
    }
}

struct App {
    options: Options,
    started: Instant,
    /// The small borderless splash, visible from startup until the real
    /// window is ready to take over.
    splash: Option<splash::Splash>,
    /// The real window, hidden until the first redraw builds `state`.
    opening: Option<splash::Opening>,
    state: Option<State>,
}

impl App {
    /// Events before the GPU is up. The splash's first redraw builds
    /// `state` (most of a second); its last painted frame stays on screen
    /// meanwhile.
    fn open(&mut self, event_loop: &ActiveEventLoop, id: WindowId, event: WindowEvent) {
        if self.opening.as_ref().is_some_and(|o| o.window.id() == id) {
            // Hidden: only a session end can close it this early.
            if let WindowEvent::CloseRequested = event {
                event_loop.exit();
            }
            return;
        }
        let Some(splash) = &self.splash else { return };
        if id != splash.id() {
            return;
        }
        match event {
            WindowEvent::Resized(_) => splash.paint(),
            WindowEvent::RedrawRequested => {
                splash.paint();
                let opening = self
                    .opening
                    .take()
                    .expect("main window created with splash");
                let fixed_size = self.options.scripted();
                let mut on_progress = |status: &str| splash.set_status(status);
                let state = State::new(opening, &self.options, &mut on_progress);
                if !fixed_size {
                    state.window.set_maximized(true);
                }
                state.window.set_visible(true);
                state.window.request_redraw();
                self.state = Some(state);
                self.splash = None;
            }
            _ => {}
        }
    }
}

impl ApplicationHandler for App {
    fn resumed(&mut self, event_loop: &ActiveEventLoop) {
        if self.state.is_none() && self.splash.is_none() {
            let layout = layout::load_layout();
            let splash = splash::Splash::new(event_loop, layout.prefs.theme);
            self.opening = Some(splash::Opening::new(event_loop, layout));
            // Not every platform sends an initial `RedrawRequested` for a
            // freshly created window; with `ControlFlow::Wait` (no
            // implicit ticking) this is the only guaranteed first frame.
            splash.window.request_redraw();
            self.splash = Some(splash);
        }
    }

    fn window_event(&mut self, event_loop: &ActiveEventLoop, id: WindowId, event: WindowEvent) {
        if self.state.is_none() {
            self.open(event_loop, id, event);
            return;
        }
        let Some(state) = self.state.as_mut() else {
            return;
        };
        if id != state.window.id() {
            // A stray event for the now-dropped splash window.
            return;
        }
        let response = state.egui_state.on_window_event(&state.window, &event);
        // `egui_winit` reports `repaint: true` for `RedrawRequested` itself
        // (it "may" need one, generically) - honoring that here would
        // self-perpetuate a redraw loop forever regardless of
        // `ControlFlow::Wait`. Whether the *next* frame is needed after a
        // redraw is `state.redraw()`'s own, more precise call to make.
        if response.repaint && !matches!(event, WindowEvent::RedrawRequested) {
            state.window.request_redraw();
        }
        match event {
            WindowEvent::CloseRequested => {
                if state.may_quit(self.options.scripted()) {
                    save_on_exit(&self.options, state);
                    event_loop.exit();
                }
            }
            WindowEvent::Resized(size) => {
                state.resize_window(size.width, size.height);
                state.window.request_redraw();
            }
            // Drag and drop from the desktop: one event per file. A
            // `.vviz`/`.json` is a session, anything else a structure.
            WindowEvent::DroppedFile(path) => state.drop_file(path),
            WindowEvent::HoveredFile(_) => state.set_drop_hover(true),
            WindowEvent::HoveredFileCancelled => state.set_drop_hover(false),
            WindowEvent::RedrawRequested => {
                state.redraw();
                if state.take_quit_request() && state.may_quit(self.options.scripted()) {
                    save_on_exit(&self.options, state);
                    event_loop.exit();
                    return;
                }
                if let Some(limit) = self.options.exit_after {
                    if self.started.elapsed().as_secs_f32() > limit {
                        eprintln!(
                            "exiting after {limit} s: {} frames, {:.1} fps",
                            state.frame_times.len(),
                            1e3 / (state.frame_times.iter().sum::<f32>()
                                / state.frame_times.len().max(1) as f32)
                                .max(1e-3)
                        );
                        save_on_exit(&self.options, state);
                        event_loop.exit();
                    } else {
                        // A benchmark/smoke-test run: keep ticking so the
                        // reported fps reflects sustained rendering, not
                        // whatever idle-detection `redraw()` decided.
                        state.window.request_redraw();
                    }
                }
            }
            _ if response.consumed => {}
            _ => {}
        }
    }
}

fn parse_options() -> Result<Options, String> {
    let mut options = Options {
        file: None,
        adapter: None,
        vsync: true,
        exit_after: None,
        exec: Vec::new(),
        listen: None,
    };
    let mut it = std::env::args().skip(1).peekable();
    while let Some(flag) = it.next() {
        match flag.as_str() {
            "--adapter" => options.adapter = Some(it.next().ok_or("--adapter needs a value")?),
            "--exit-after" => {
                let v = it.next().ok_or("--exit-after needs a value")?;
                options.exit_after = Some(v.parse().map_err(|e| format!("--exit-after: {e}"))?);
            }
            "--exec" => {
                let script = it.next().ok_or("--exec needs a script")?;
                for line in vv_scene::split_script(&script) {
                    commands::validate_unattended(&line).map_err(|e| format!("--exec: {e}"))?;
                    options.exec.push(line);
                }
            }
            "--uncapped" => options.vsync = false,
            "--listen" => {
                let port = match it.peek() {
                    Some(p) if !p.starts_with('-') && p.parse::<u16>().is_ok() => {
                        it.next().and_then(|p| p.parse().ok()).unwrap_or(0)
                    }
                    _ => 0,
                };
                options.listen = Some(port);
            }
            "-h" | "--help" => {
                println!("{USAGE}");
                std::process::exit(0);
            }
            other if other.starts_with('-') => {
                return Err(format!("unknown option `{other}`\n\n{USAGE}"))
            }
            file => options.file = Some(PathBuf::from(file)),
        }
    }
    Ok(options)
}

/// On Windows a GUI-subsystem process has no console; if the parent (a
/// shell) has one, borrow it so `eprintln!` reaches the terminal. Fails
/// harmlessly when there is no parent console (double-click launch).
/// Rust's stdio looks the handles up on every write, so no reopen is
/// needed. Known limit: `cmd.exe` and PowerShell do not wait for GUI
/// processes, so scripted runs there want `start /wait` or a pipe.
#[cfg(windows)]
fn attach_parent_console() {
    use windows_sys::Win32::System::Console::{AttachConsole, ATTACH_PARENT_PROCESS};
    // SAFETY: a plain Win32 call with no pointers; the result is ignored on
    // purpose (no parent console is the normal double-click case).
    unsafe {
        AttachConsole(ATTACH_PARENT_PROCESS);
    }
}

#[cfg(not(windows))]
fn attach_parent_console() {}

fn main() {
    attach_parent_console();
    let options = match parse_options() {
        Ok(o) => o,
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(2);
        }
    };
    let event_loop = EventLoop::new().expect("event loop");
    // Redraws are self-scheduled (state.rs's `redraw()` tail, plus explicit
    // `request_redraw()` calls on input that changes what's drawn): the
    // loop blocks between them instead of ticking continuously, so a
    // static viewport costs no CPU/GPU work at all.
    event_loop.set_control_flow(ControlFlow::Wait);
    let mut app = App {
        options,
        started: Instant::now(),
        splash: None,
        opening: None,
        state: None,
    };
    event_loop.run_app(&mut app).expect("run event loop");
    if app.state.as_ref().is_some_and(|s| s.exec_failed) {
        eprintln!("--exec: a command failed; the rest of the script was skipped");
        std::process::exit(1);
    }
}
