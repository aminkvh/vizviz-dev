# Movie encoder spike: H.264 (openh264) vs AV1 (rav1e)

Question: can the movie maker write MP4 without an external `ffmpeg`, on
Windows, Linux and macOS, fast and with good quality? Measured on **Windows
11, x86_64-pc-windows-msvc, 16 logical cores, release build**. Linux and
macOS were **not** run; they are reasoned from crate docs and marked so.

Harness: `xtask/encoder-spike/` (its own workspace, not part of the app or
`cargo deny`). It reads the PNG frames of a real `movie export` and encodes,
decodes and scores them.

## Inputs

Two real exports from the release `vizviz.exe`, 1280x720, 30 fps, opaque
frames, "Smooth edges" on, light background:

| Clip | Content | Frames | Edge pixels |
|---|---|---|---|
| cartoon | 1AKE cartoon, coloured by chain, one turn in 5 s | 150 | 2.8 % |
| lines | 1AKE `lines`, rainbow, half a turn in 3 s (thin coloured lines) | 90 | 4.9 % |

## Method

- Colour: BT.709, limited range, 4:2:0. Chroma is the mean of each 2x2
  block going down; bilinear (centre-sited) going up. Both codecs get the
  same YUV and are tagged BT.709 (H.264 VUI, AV1 colour description).
- **PSNR** = 10 log10(255^2 / MSE) per frame, MSE over all samples (Y: luma
  plane; RGB: three channels), capped at 100 dB, averaged over frames.
- **SSIM**: standard, C1 = (0.01*255)^2, C2 = (0.03*255)^2, 8x8 uniform
  windows every 4 px; RGB is the mean of the three channels.
- **Edge PSNR**: the same MSE over source "edge pixels" only (any RGB
  channel differing by more than 32 from the right or lower neighbour, plus
  that neighbour). This is where thin coloured lines and chroma fringing show.
- **Subsampling floor** (row "4:2:0 round trip"): the source converted to
  4:2:0 and back with no codec. It is the best any 4:2:0 encoder can score
  on RGB metrics here.
- Speed is encode time only (RGB to YUV conversion and muxing excluded;
  mux took 4-8 ms). Each row is a single run, so fps carries run-to-run
  noise of tens of percent (openh264 camera mode at 4 Mbps happened to be
  run twice: 236 and 330 fps). Target bitrate is what was asked; **actual Mbps is what came out**.
- openh264: High profile, High complexity, BT.709 VUI, frame skipping on
  where noted (skipped 0 of every clip), threads left at the default.
  rav1e: `speed` preset 10 and 6, default threading, `asm` feature.
- Decoding back: H.264 was muxed with the `mp4` crate, read back through
  the same crate and decoded with openh264's own decoder, frame count
  asserted equal. rav1e **has no decoder**; its numbers are the encoder's
  own reconstruction (`Packet::rec`), which a conforming decoder reproduces,
  but the IVF was never independently decoded.
- An independent decoder for the H.264 file: Windows Media Foundation (via
  WPF `MediaPlayer`, Windows PowerShell 5.1) opened the 8 Mbps MP4, reported
  1280x720 and 5 s, and rendered a frame that scored 40.85 dB RGB PSNR
  against the nearest source frame (mean channel bias under 2 levels).
  That is a sanity check that the file plays and the colour tagging is
  honoured, not a quality measurement (frame alignment and MF's own chroma
  upsampling are unknown). ffmpeg/ffprobe are not installed here.

## Results: cartoon clip (150 frames)

| Encoder | Target | Enc fps | Actual Mbps | MB | PSNR Y | PSNR RGB | SSIM RGB | Edge PSNR |
|---|---|---|---|---|---|---|---|---|
| 4:2:0 round trip, no codec | - | - | - | - | 100 | 43.42 | 0.9977 | 29.50 |
| openh264 (camera mode) | 2 | 256 | 2.0 | 1.25 | 44.6 | 37.72 | 0.9902 | 24.43 |
| openh264 (camera mode) | 4 | 236 | 4.0 | 2.51 | 49.6 | 40.49 | 0.9941 | 27.05 |
| openh264 (camera mode) | 8 | 306 | 7.9 | 4.95 | 55.4 | 42.40 | 0.9961 | 28.81 |
| openh264 (camera mode) | 20 | 315 | 8.4 | 5.25 | 56.0 | 42.53 | 0.9962 | 28.92 |
| openh264 screen-content mode | 8 | 170 | 2.4 | 1.48 | 45.8 | 38.47 | 0.9919 | 25.15 |
| rav1e speed 10 | 4 | 8.6 | 4.0 | 2.50 | 45.1 | 39.25 | 0.9897 | 25.79 |
| rav1e speed 10 | 8 | 7.0 | 8.0 | 5.00 | 51.3 | 41.67 | 0.9943 | 28.18 |
| rav1e speed 10 | 20 | 5.4 | 19.4 | 12.11 | 65.3 | 43.27 | 0.9974 | 29.39 |
| rav1e speed 6 | 4 | 2.5 | 4.0 | 2.50 | 52.4 | 41.99 | 0.9966 | 28.02 |
| rav1e speed 6 | 8 | 2.0 | 8.0 | 4.99 | 60.4 | 43.04 | 0.9974 | 29.14 |
| rav1e speed 6 | 20 | 1.7 | 13.1 | 8.18 | 66.6 | 43.35 | 0.9976 | 29.44 |

## Results: thin coloured lines (90 frames)

| Encoder | Target | Enc fps | Actual Mbps | PSNR RGB | SSIM RGB | Edge PSNR |
|---|---|---|---|---|---|---|
| 4:2:0 round trip, no codec | - | - | - | 41.29 | 0.9945 | 28.95 |
| openh264 (camera mode) | 8 | 318 | 8.0 | 39.42 | 0.9885 | 27.30 |
| openh264 (camera mode) | 20 | 306 | 13.1 | 40.43 | 0.9912 | 28.28 |
| rav1e speed 10 | 8 | 7.3 | 8.2 | 39.42 | 0.9877 | 27.37 |
| rav1e speed 6 | 8 | 2.2 | 8.0 | 40.35 | 0.9923 | 28.04 |
| rav1e speed 6 | 20 | 1.6 | 18.3 | 41.21 | 0.9943 | 28.87 |

Full runs (camera, quality-mode and screen-content openh264 variants at
2/4/8/20 Mbps; rav1e speed 10 and 6) were run for both clips; the rows above
are representative.

## What the numbers say

- **openh264 is fast**: about 240-400 fps at 720p (nasm on PATH, so its
  assembly was used; the crate's README claims about 3x slower without
  nasm, not measured here). rav1e is 2-9 fps on the same 16 cores: a
  5 s clip takes 17-90 s. openh264 is roughly 30-150x faster.
- **Chroma subsampling dominates the error.** The no-codec 4:2:0 floor is
  43.4 dB RGB (29.5 dB on edge pixels); at 8 Mbps openh264 is within about
  1 dB of it and rav1e speed 6 within 0.4 dB. Past that, more bitrate cannot
  fix fringing on thin coloured lines: only 4:4:4 (or not subsampling) can.
  4:4:4 was not tested (openh264's encoder lists High444 but the crate's
  I420-only buffers were used; AV1 4:4:4 plays in few players).
- **openh264 stops spending bits** at about 8.4 Mbps on the cartoon clip
  (13 Mbps on lines): its rate control hits its lowest quantizer, so a
  "20 Mbps" request gives 8.4. It hit 2, 4 and 8 Mbps targets closely.
  Its screen-content mode was worse at any rate here (about 2.4 Mbps cap,
  lower PSNR) and half the speed: do not use it for rendered frames.
- **rav1e is more efficient per bit** at low rates (4 Mbps speed 6: 42.0 dB
  vs openh264 40.5 dB) but rav1e speed 10 is no better than openh264 at the
  same rate (39.25 vs 40.49 at 4 Mbps). At 8 Mbps both are near the floor.
- Build cost: `openh264` with the `source` feature compiled from C++ with
  the MSVC compiler found by `cc`, plus `nasm` if present; the whole spike
  (both codecs, 4 jobs) built in about 1.5 min. The two openh264 modes gave
  **bit-identical output** (same sizes and scores).

## Build and run needs, per mode

| | Windows (measured) | Linux / macOS (from the crate README, not run) |
|---|---|---|
| openh264 `source` | C++ compiler (cc finds MSVC); nasm optional (used here) | C++ compiler; nasm optional; README lists x86_64/aarch64 Linux and macOS as compiled, x86_64 as unit-tested |
| openh264 `libloading` | needs Cisco's `openh264-2.6.0-win64.dll` (a 452 KB bz2 download, sha256 checked by the crate against a whitelist; worked) next to the exe or by path | same with the `.so` / `.dylib` from Cisco |
| rav1e | pure Rust plus nasm for the `asm` feature; without it much slower (not measured) | same; nasm on x86 |

The `mp4` crate (MIT) muxed H.264 with no extra work (AVCC length-prefixed
samples, SPS/PPS from the first packet); the result opened in Media Foundation
and in the crate's own reader. It has no AV1 track writer; rav1e output went
to a hand-written IVF (32-byte header, 12-byte frame headers), not MP4.
`minimp4` was rejected on licence (MPL-2.0).

## Licence and patents (not legal advice)

- **openh264**: crate and Cisco code are BSD-2-Clause. H.264 is covered by
  a patent pool; Cisco pays the royalties only for **its own prebuilt binary
  downloaded at install/run time**, so a source-built openh264 shipped in our
  installers gets no such cover. Whoever distributes it carries that risk.
- **rav1e**: BSD-2-Clause with the AOMedia royalty-free patent grant; a
  third-party pool asserts claims on AV1 that AOMedia disputes. AV1-in-MP4
  is also less widely playable than H.264 (older players, some editors).
- `cargo deny` on the spike's tree: every crate in the built graph is
  MIT/Apache/BSD/Zlib/Unlicense/0BSD-family. It reports one rejection,
  `libfuzzer-sys` (NCSA, permissive) from rav1e's optional fuzzing feature,
  which is not compiled but is pulled in by the config's `all-features`.

## Recommendation

**openh264 + the `mp4` crate** for a built-in MP4: 30-150x faster than
rav1e, quality within about 1 dB of the 4:2:0 ceiling from 8 Mbps at 720p,
plays everywhere, MIT/BSD dependencies only. Pick the bitrate from the
frame size (about 8 Mbps at 720p for these scenes, scaled with pixel count)
and leave frame skipping off. Decide the patent question before shipping:
the runtime-loaded Cisco DLL gives royalty cover but needs a download and
a per-platform blob; the source build is self-contained but uncovered.
rav1e is worth offering only as a slow "small file / archival" option later;
it is not fast enough to be the default.

## Gaps

- Linux and macOS untested; openh264 not built without nasm; single timing runs only;
  only two clips; no 4:4:4 test; H.264 profile/level signalling
  was left to openh264 (no check that every player accepts High at any size).
- rav1e output not decoded by an independent AV1 decoder.
