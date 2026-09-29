# Movies

The Movie panel lays camera moves and trajectory playback on a timeline and
writes them out as a numbered PNG sequence, an animated GIF, and an MP4 if
you have `ffmpeg`. Open it from File ▸ Export ▸ Movie…, View ▸ Panels, or `movie`.

## The idea

A movie has a length (seconds), a frame rate and a picture size. It holds
**tracks**; each track holds **clips** placed on one shared time axis.
Clips on different tracks run at the same time. Clips on one track can
touch but never overlap.

| Track | A clip on it |
|---|---|
| Rotate | turns the scene by some degrees about the x, y or z screen axis |
| Zoom | scales the camera distance from one factor to another (1 = as is, below 1 = closer) |
| Move to view | flies the camera to a view (the one on screen when you add the clip) |
| Play trajectory | steps a structure through a frame range at some frames per second, once or looping |

Every clip has an easing: `linear`, or `ease` (slow start and end).
Easing applies to the camera clips, not to trajectory playback.

## What a moment looks like

Camera clips apply in track order (top to bottom) to the camera you had at
time 0. A finished clip stays applied: a 90 degree rotation that ended at
2 s is still there at 5 s. A *Move to view* clip replaces what the tracks
above it did, so put it above any rotation you want on top of it. Before a
trajectory clip starts, its structure shows the clip's first frame; after
it ends, the last frame it reached.

Scrubbing or playing changes the live camera. If you then move the camera
yourself, that view becomes the new start (time 0) and the playhead goes
back to 0.

## In the panel

- The ruler is the scrubber: click or drag it. The transport buttons play
  and rewind. The row above sets fps, length and size.
- **Add track** adds a lane. Each lane's **+** adds a 2 s clip at the
  playhead, in the first free gap.
- Drag a bar to move it, drag its ends to resize it. It stops at its
  neighbours and at the ends.
- The **Clip** row edits the selected clip: start, length, easing and its own
  settings.
- **Export…** asks for a folder and writes `frame_0001.png`, `frame_0002.png`
  … there, one movie frame at a time, at the movie's size. If `ffmpeg` is on
  your PATH it also encodes `movie.mp4` (the panel says when it is not).
  "Smooth edges" renders at twice the size and averages down, as
  screenshots do. Switch off "PNG frames" to skip the sequence (and the
  MP4, which is built from it).
- **Also GIF** writes `movie.gif` too, built into vizviz (no ffmpeg). It is
  shrunk to at most the **GIF width** (800 px by default, never enlarged),
  averaging in linear light. Each frame gets its own 256-colour palette and
  only the part that changed is stored, so a plain background stays exact
  and the file stays small. **Dither** (on by default) spreads the colour
  rounding into a fine pattern; turn it off for flat colours. It loops
  forever. GIF times are whole hundredths of a second, so 30 fps plays as
  3, 4, 3, 3, 4 … hundredths (right on average). Browsers cannot play faster
  than 50 fps, so a higher movie rate drops frames to 50.

Your camera and frames are put back when the export ends. While a
trajectory plays, cartoons keep their secondary structure and only follow
the atoms, in the export as in the live view; a `screenshot` of one frame
recomputes it.

## From a script

Tracks and clips count from 1.

```
vizviz --exec "loadtrajectory 1CRN.pdb 1CRN_traj.dcd; \
  movie set duration 4; movie set fps 24; movie set size 1280x720; \
  movie addtrack rotate; movie clip 1 0 4 y 360 linear; \
  movie addtrack play;   movie clip 2 0 4 current 0 2 3 loop; \
  movie export frames; quit"
```

`movie clip TRACK START LENGTH ...` takes, by track kind:

| Kind | Arguments |
|---|---|
| rotate | `x\|y\|z DEGREES` |
| zoom | `FROM TO` |
| moveto | none (uses the view on screen) |
| play | `#N\|current FIRST LAST FPS [loop]` (`#N` as `structures` lists it) |

Add `linear` or `ease` (default) at the end. `movie time S` shows that
moment in the viewport, so `movie time 2; screenshot at2s.png` grabs a
still. `movie export DIR [nossaa] [nomp4] [nopng] [gif] [nodither]
[gifwidth=N]` runs in the background and a script waits for it. Plain
`export` keeps writing PNGs (and the MP4 if ffmpeg exists); add `gif` for
`movie.gif`. `gifwidth=N` (16 or more) and `nodither` only apply with
`gif`; `nopng` needs `gif` and drops the MP4 too.

## Saved and not saved

The movie is saved inside sessions. Not there yet: keyframe curves, audio,
per-clip undo, showing/hiding or restyling parts over time.
