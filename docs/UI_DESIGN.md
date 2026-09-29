# Interface design principles

Rules, not suggestions: retrofitting UI philosophy onto ad-hoc panels
later is expensive. Concrete tokens and components are in
[DESIGN.md](../DESIGN.md); recurring interaction patterns in
[UX_PATTERNS.md](UX_PATTERNS.md).

## Docking model: professional creative tool

Reference point: image, video and audio editors, whose panels dock,
float and save as workspaces.

- Dockable panels: left, right, top (bottom optional/collapsible), each
  resizable and collapsible.
- Any panel can float/undock (multi-monitor is common for this audience).
- Savable, switchable workspace layouts (editing vs. trajectory analysis
  vs. presentation shouldn't require manual rearranging every time).
- No modal dialogs except truly blocking actions (destructive confirms,
  save prompts). Everything else lives in a panel.
- The 3D viewport is itself a dock tab (rendered to a texture), so it
  docks, floats, and splits like any other panel; multiple viewports later.

## Interaction: select-to-analyze

- Click an atom/residue/chain → a contextual inspector panel updates
  immediately, reacting to whatever is active.
- Selections are composable, named, and persistent — not a transient
  highlight lost on the next click.
- Analysis from a selection docks next to it, not in a disconnected window.

## Visual design

- Minimal chrome by default; the 3D viewport is the primary surface.
  Toolbars/panels are opt-in, not maximal by default.
- One color palette, one spacing scale, one icon set — a real design
  system, not per-feature styling (concrete tokens defined when UI code
  starts).
- Light and dark themes both supported (light for figures/screenshots,
  dark for long sessions).
- Render style presets ("skins") switchable in one click: publication
  (white, flat lighting, outlines), dark presentation, cel/flat, glossy.
  A preset is a bundle of shading/background/outline parameters, so
  adding one is data, not code.

## Annotations & figures

Publication/figure quality is a first-class output, not a side effect of
a screenshot. Two distinct annotation types, both needed:

- **Structure-anchored** — labels, arrows, highlights, measurements bound
  to real 3D coordinates (an atom/residue/region). Must track camera
  rotation and zoom correctly (billboarded text, world-space lines); open
  question, to decide at implementation time, is whether an anchor tracks
  a moving atom during trajectory playback or freezes at its initial
  position — both are legitimate depending on use case, so this is
  probably a per-annotation toggle, not a global rule.
- **Screen-space** — free-floating text/shapes for general figure
  composition, not tied to any structure coordinate.
- Export to vector formats (SVG/PDF), not just raster — figures need to be
  editable in a vector editor afterward.

## Search / command access

Fast, fuzzy, keyboard-first command palette (Ctrl/Cmd+P) searching
commands, structures, residues, settings, and plugin actions from one
entry point — not a nested menu hunt.

## Architecture consequence

Panels, selection state, and the command palette all react to the same
underlying state, so the UI needs one shared application-state model that
panels subscribe to, not per-panel independent queries. Decide the concrete
mechanism (egui + a shared state struct is the natural fit) when `egui`
integration work starts.

## Panels

- **Viewport** — every loaded structure in one image; click an atom or bond
  to select it, Ctrl+click adds, right-drag pans, F resets the view.
- **Sequence** (above the viewport) — one-letter codes per chain; click a
  residue to select it. Virtualized, so huge structures cost nothing.
- **Scene** — loaded structures with per-structure representation/coloring.
- **Selection** — an expression box (docs/SELECTION.md), the active
  selection, and saved sets (which remember their expression).
- **Inspector** — one atom's details, or a distance/bond length for two.
- **Info** — the loaded file's annotations: title, method, resolution,
  organism, entities, citation.
- **Scene Settings** — style preset, background color, FXAA, outline,
  supersampling, culling, adaptive quality, point threshold.
- **Timeline** — opens itself when a multi-frame file loads: scrub, play,
  loop.
- **Log** — what the command bus did, and why something failed.
- Workspaces: Default, Trajectory (adds the timeline), Analysis (selection
  front and centre), Compact; plus any the user saves.
