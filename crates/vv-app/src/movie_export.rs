//! Naming and encoding for a movie export: a numbered PNG sequence, an
//! MP4 of it when `ffmpeg` is on the PATH, and a GIF (`movie_gif`). The
//! frame rendering itself is `state`'s (it owns the GPU); this is only the
//! file side.

use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::OnceLock;

/// The MP4 written next to the frames.
pub const MP4_NAME: &str = "movie.mp4";
/// The GIF written next to the frames.
pub const GIF_NAME: &str = "movie.gif";

fn digits(total: usize) -> usize {
    total.to_string().len().max(4)
}

/// `frame_0001.png` for the first of `total` frames.
pub fn frame_path(dir: &Path, index: usize, total: usize) -> PathBuf {
    dir.join(format!(
        "frame_{:0width$}.png",
        index + 1,
        width = digits(total)
    ))
}

/// Whether `ffmpeg` runs, asked once.
pub fn ffmpeg_available() -> bool {
    static FOUND: OnceLock<bool> = OnceLock::new();
    *FOUND.get_or_init(|| {
        Command::new("ffmpeg")
            .arg("-version")
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .is_ok_and(|s| s.success())
    })
}

/// ffmpeg's arguments for encoding the sequence `frame_count` frames long
/// into `MP4_NAME`, both inside `dir`.
pub fn ffmpeg_args(dir: &Path, fps: f32, frame_count: usize) -> Vec<String> {
    let pattern = dir.join(format!("frame_%0{}d.png", digits(frame_count)));
    let out = dir.join(MP4_NAME);
    let text = |p: &Path| p.to_string_lossy().into_owned();
    [
        "-y",
        "-hide_banner",
        "-loglevel",
        "error",
        "-framerate",
        &fps.to_string(),
        "-start_number",
        "1",
        "-i",
        &text(&pattern),
        "-c:v",
        "libx264",
        "-pix_fmt",
        "yuv420p",
        &text(&out),
    ]
    .map(String::from)
    .to_vec()
}

/// Starts encoding in the background; poll the child for the end.
pub fn start_encoding(dir: &Path, fps: f32, frame_count: usize) -> std::io::Result<Child> {
    Command::new("ffmpeg")
        .args(ffmpeg_args(dir, fps, frame_count))
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frames_are_numbered_from_one_with_room_for_the_total() {
        let dir = Path::new("out");
        assert_eq!(frame_path(dir, 0, 300), dir.join("frame_0001.png"));
        assert_eq!(frame_path(dir, 299, 300), dir.join("frame_0300.png"));
        assert_eq!(frame_path(dir, 0, 12_345), dir.join("frame_00001.png"));
    }

    #[test]
    fn ffmpeg_reads_the_same_pattern_the_frames_are_written_with() {
        let dir = Path::new("out");
        let args = ffmpeg_args(dir, 30.0, 12_345);
        let input = &args[args.iter().position(|a| a == "-i").unwrap() + 1];
        assert!(input.ends_with("frame_%05d.png"), "{input}");
        assert!(args.iter().any(|a| a == "yuv420p"));
        assert!(args.last().unwrap().ends_with(MP4_NAME));
    }
}
