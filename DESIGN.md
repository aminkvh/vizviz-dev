# vizviz design system

One system for every panel, dialog and overlay. If something needs a value
or component that is not here, add it here first. It is implemented in
`crates/vv-app/src/theme.rs` (tokens) and `crates/vv-app/src/widgets.rs`
(components); nothing else sets a colour, size or spacing.

vizviz is a desktop app (egui), not a web page: "tokens" are Rust
constants, and sizes are logical pixels (scaled by the OS).

## Spacing

A 4 px scale, nothing else: **4, 8, 12, 16, 24, 32**.

- 4: between an icon and its label; between a label and its field.
- 8: between controls in a row or column (default item spacing).
- 12: inside cards and dialogs, and between a control and its help text.
- 16: panel padding; between sections of a panel.
- 24 / 32: dialog padding; around empty states.

## Type

Six styles, one family (the egui default, Ubuntu-Light metrics):

| Style | Size | Weight / colour | Use |
|---|---|---|---|
| Caption | 11 | regular, muted | ribbon group names, hints under fields, viewport overlays |
| Body | 13 | regular, primary | everything by default |
| Label | 13 | strong (primary, brighter) | field labels, current item |
| Section | 15 | strong | a panel's section headings |
| Title | 18 | strong | dialog and empty-state titles |
| Mono | 12 | regular | selection expressions, the log, file paths |

Line height is 1.4 x size (egui row height). Icons (Phosphor) are 18 in
icon-only buttons and 24 over a large button's label. Sentence case
everywhere ("Open session", not "Open Session"), American spelling
("color"); an ellipsis (…) only on actions that ask for more input before
they act. Counts are words: "1 atom", "4,779 atoms", never "atom(s)". A
test holds every command title and help text to these rules.

## Color

Semantic tokens, dark and light. Every text/background pair meets WCAG AA
(4.5:1), checked by a test.

| Token | Use |
|---|---|
| `bg` | behind docked panels |
| `surface` | panels, dialogs |
| `surface_raised` | buttons, inputs, cards |
| `border` | control and panel outlines, dividers |
| `text` | primary text |
| `text_muted` | secondary text, captions, placeholders |
| `primary` | the one primary action, selected state, links, focus ring |
| `on_primary` | text on `primary` |
| `danger` | destructive actions, errors |
| `success` | done / confirmed (the logo's teal) |
| `warning` | slow or partial results (e.g. "left out of the render") |

States, derived from each token: **hover** = `surface_raised` → `hover`
fill; **active/pressed** = hover fill + `primary` border; **selected** =
`primary` fill, `on_primary` text; **disabled** = 40 % opacity, no hover;
**focus** = 2 px `primary` ring outside the control (never removed).

Colors inside the viewport (atoms, labels, measurements, the view cube's
red-green-blue axes) are data, not chrome, and are out of scope here.

## Radius

Two values: **4** (controls: buttons, inputs, chips) and **8** (containers:
cards, dialogs, menus, floating overlays).

## Components

Every one lives in `widgets.rs` with all of its states.

**Button**: height 28 (hit target 28 x at least 28; WCAG 2.2 asks 24).
Variants:
- *Primary*: `primary` fill. At most one per panel or dialog: the action
  the user came to do ("Render", "Open", "Save session").
- *Secondary*: `surface_raised` fill, `border` outline. Everything else.
- *Ghost*: no fill until hover. Toolbar and ribbon actions, icon buttons.
- *Danger*: `danger` text, danger fill on hover; always confirms (see below).
Sizes: regular (32), large (ribbon's main action per group: 56 tall, icon
over label, as wide as its label), icon-only (28 x 28). Icon-only
buttons carry a tooltip that says what they do; every ribbon action's
tooltip is its command's help, so labels stay short ("List" under
Measure, not "List measurements").

**Toggle** (on/off setting): a switch, never a highlighted button, so "on"
never looks like "selected". **Segmented choice** (one of a few): joined
buttons, the current one selected. **Select** (one of many, e.g.
materials): a drop-down with readable names ("Glass 1", not `glass1`). A
row with no room for even a short label (a Structures panel selection's
style/coloring/material, three to a line) shows a fixed icon instead of
its text — the same icon its ribbon equivalent uses, so the two read as
one system — with the current value in its tooltip; the open list still
spells every option out in full.

**Input**: height 28, `surface_raised`, label above or to the left, hint
below in Caption, error below in `danger` saying what to fix. Keeps its
text after a failure.

**Card** (the start card over an empty viewport): `surface`, `border`,
radius 8, padding 24. **Dialog**: modal, radius 8, padding 24, Title,
body, then buttons bottom-right: secondary "Cancel", then the verb
("Fetch", "Save"; Danger for a destructive one). Escape cancels, Enter
runs the verb.

**List row** (a structure, a rep, a set, a palette entry): one line cut
with "…" to the room its trailing icon buttons leave, the full text on
hover; the current one selected.

**Slider**: its label on its own line above, the value with two decimals
and its unit beside it.

**Table** (key/value, e.g. file info): two columns, labels in `text_muted`
right-aligned to a fixed column, values selectable.

**Toast** (the notice bar): bottom of the viewport, info (`primary`
icon) or error (`danger` icon and border); says what happened and what
to do next ("Could not read 4hhb.cif: not mmCIF or PDB. Check the file
or open another"). An undoable action's toast carries an Undo button
that runs `undo`. Info fades after 4 s; errors stay until dismissed. One
bar: a new toast replaces whatever was showing, never queues.

**Empty state**: icon, one line saying what goes here, one primary action
(e.g. Structures: "No structures yet" + "Open a structure").

**Loading state**: the notice bar says what is building ("building the
molecular surface for 4HHB; it appears when ready") and, when known, how
far it is ("rendering 43 %").

**Viewport overlays**: the mouse mode in a chip at the top left (frame
timings on hover), the view cube with its axes at the top right, the
scale bar at the bottom left, and a `primary` mark on every selected
atom, over whatever draws it.

## Layout

- **Ribbon** on top: exactly five tabs, left to right in the order of
  work — File ("get things in and out"), Represent ("how this rep is
  drawn"), Look ("how the whole picture is lit"), Analyze ("measure and
  annotate"), View ("where I'm looking from and what I see"). No
  contextual tabs: playback lives in the Timeline panel, surface
  materials in Represent. Each group's main action is the large button on
  its left. A fixed **quick-access bar** at the ribbon's top right —
  Undo, Redo, Save session, Find commands (Ctrl+P), Preferences (gear,
  Ctrl+,) — never moves with the tabs.
- **Popover** (a `NAME ▾` button): a settings group too small for its
  own panel (Look's Lights/Background/each effect, View's Clip/Layout/
  Panels). An on/off setting with a strength behind it is a switch
  *and* its own ▾, shown only while on — the switch alone is
  still a plain command; the ▾ never hides the only way to reach a
  setting. Every popover with more than one field gets its own Reset.
- **Form popover**: the same `NAME ▾` shape, for a "create" action that
  needs a couple of typed fields (View's Rotate, File's Add trajectory,
  File's Render) — labelled fields, then the primary verb and Cancel
  (`widgets::form_actions`). A full dialog is only for file choosers
  and anything that needs to block on a decision (Preferences' Danger
  group, deleting a workspace).
- **Viewport** in the centre, always visible, no header of its own (no
  ▼, no +; nothing docks into it) — the only dark-by-default surface in
  the light theme. Its top-left mode chip names the current mouse mode
  and has a "?" listing every mouse binding.
- **Left dock**: what is loaded and what is selected (Structures,
  Selection). **Right dock**: details of the current thing (Inspector,
  Info). **Bottom**: Timeline and Log.
- Each panel: 12 padding, 4 between rows, scrolling when it runs long; a
  Section heading per group; its primary action at the bottom of its
  form. One close button per panel, on its tab (the viewport has none);
  each tab bar's + reopens a closed panel there; panels collapse and
  resize. Launches open the default arrangement, or the one
  `startlayout` pinned.
- **Preferences** (gear, Ctrl+,, `preferences`): the one settings
  screen, three groups — Only me (theme, start card, quit reminder,
  start layout, clear recent files), Performance (adaptive quality and
  point size side by side, occlusion culling, supersample, FXAA),
  Danger (reset all preferences; delete saved layouts). Nothing a user
  needs mid-task lives here — lights, background, clip stay in Look/View.
- **Structures panel**: one compact row per loaded structure (name, atom
  count, an eye for the whole structure, ⋯ for rename/close/file info/
  export); its rows map 1:1 to that structure's reps, one row per
  selection (its expression, compact style/coloring/material selects, an
  eye, ⋯ for Select these atoms/Duplicate/Delete/Options) — the same rep
  model the Represent ribbon tab edits, so the two stay in sync.
- **Representation options** (atom size, bond and tube radius, probe,
  blobbiness, shrink) live in a selection row's own ⋯ ▸ Options popover,
  the first couple always shown and the rest (only Tube's three today)
  folded under "Advanced".
- **One home per setting.** A setting lives in one place (ribbon or
  panel); anything that appears in both uses the same control and wording.
- Window sizes: designed for 1280-2560 wide; the window cannot shrink
  below 1024 x 640. Below 1280 the ribbon folds to its tab row, and a tab
  click shows that tab until its next action runs. vizviz is a desktop
  tool; phone widths (375) are not supported.

## Interaction rules

- Every action gives visible feedback within 100 ms (pressed state, then
  a result or a loading state).
- Destructive or irreversible actions confirm first, naming what is lost
  ("Delete workspace 'Figures'? This can't be undone."). Undoable actions
  do not confirm; their toast carries an Undo button ("Closed 4HHB." +
  Undo).
- Quitting with unsaved changes asks to save.
- Slow actions show a loading state and ignore repeat clicks.
- Errors say what happened and what to do next.
- Forms validate as you type (e.g. a selection expression shows the
  count it matches or the error), and keep input on failure.
- Keyboard: every action reachable (palette, key tips, shortcuts); focus
  always visible; Tab order follows reading order.
- Selecting something shows it in the viewport.
