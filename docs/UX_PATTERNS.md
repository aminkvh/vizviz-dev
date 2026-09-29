# UX patterns

One way to do each recurring thing; every screen follows it. Components
are the ones in
`crates/vv-app/src/widgets.rs` and DESIGN.md.

## Create
- **Where**: the ribbon button for the kind of thing (Add rep, Label,
  Caption, Measure, Save set, Save layout), or the list's own **+**
  at the end of its section header.
- **How**: a small **form popover** anchored to the button: labelled
  fields, sensible defaults filled in, primary button named with the
  verb ("Add rep", "Save set"), Escape cancels. Never a half-typed
  command. A full dialog only for file choosers.
- **After**: the new thing becomes current and is highlighted in its
  list; a toast says what was made, with **Undo**.

## Edit
- **Inline** for one value on a row (rename a set, a rep's "Draws"
  field): click the text, type, Enter saves, Escape reverts.
- **Popover form** for a group of values (a rep's options, clip,
  a light): opens from the item's ▾ or ⋯.
- **Saving**: every edit applies immediately (live preview while
  dragging, one undo step when released). No Save/Apply buttons except
  in file dialogs.
- **Unsaved changes**: only at the session level (quit / open another
  session asks, with "Don't remind me").

## Delete / destructive
- **Undoable** deletes (rep, structure, set, annotation, label) happen
  at once, with a toast "Deleted X — Undo". No confirmation dialog.
- **Not undoable** (delete a saved layout, reset preferences, overwrite
  a file) confirm first. Title: "Delete layout 'Analysis'?"; body says
  what is lost; buttons: "Delete" (danger) and "Cancel". Never "OK".
- **Afterward**: the list selects the neighbouring item; nothing jumps.

## Bulk actions
- Lists that allow several (sets, annotations, reps): Ctrl/Shift+click
  selects rows; a slim bar appears at the list's top with the actions
  that apply to all ("Delete 3", "Hide 3"), and Escape clears.
- Atom selection is the viewport's bulk selection: its actions live in
  the Selection panel (Save set, Clear) and Analyze (Measure).

## Search / filter / sort
- One component: a search field at the top of the list it filters
  (structures, reps, sets, annotations, palette, layouts), filter as
  you type, x to clear, Escape clears.
- Active filters show as chips under the field, each with x.
- Lists sort in the order things were made; no sort controls.

## Lists and detail
- **List**: `list_row` rows, one line each: name, a short summary, an
  eye if it can hide, ⋯ for the rest. The current row is highlighted.
- **Detail**: under the list, in the same panel ("Options of Cartoon"
  under the rep list). No separate windows.
- **Back**: none needed, detail follows the current row.
- **Compact row** (a Structures panel selection: expression, style,
  coloring, material, eye, ⋯, all on one line): no room for a detail
  section below each row, so its own small controls (`select_icon`s, a
  fixed icon instead of the current value, matching their Represent-tab
  group's own) sit directly in the row, and its rare, bigger detail (rep
  options) opens from ⋯ as a popover anchored on that row instead — the
  same "no separate window" rule, just anchored to the row rather than
  living below the whole list.

## Settings screens
- Only Preferences. One scrolling page, groups with a heading and a
  one-line description each: Only me, Performance, Danger.
- Every change applies and saves at once; a caption under a setting
  says its effect when it isn't obvious.

## Feedback
- **Toast**: bottom of the viewport, one line, 4 s, an action when
  useful (Undo, Show). Every command that changes something toasts.
- **Error**: a toast in the danger colour: what failed, and why, in one
  sentence ("Couldn't open 1abc.pdb: not a PDB or mmCIF file"). Details
  go to the Log.
- **Empty state**: the `empty_state` component: icon, one line, one
  button that fixes it ("No structures yet · Open a structure").
- **Loading**: work over ~0.3 s shows progress in the item it affects
  (a rep row's spinner, a surface keeps its old shape until the new one
  is ready); nothing blocks the window.
- **Results** (SASA, contacts, render): a result card in the Inspector
  with Copy, plus a toast.

## Navigation
- Panels are the screens; the ribbon and quick-access bar are the
  global navigation. No breadcrumbs, no back button.
- **Deep links**: every panel, popover and action has a command
  (`panel selection`, `look lights`, `clip`), usable from the palette,
  the console, scripts and Python. The palette (Ctrl+P) finds all of
  them by name.

## Keyboard
- Escape closes the top popover/dialog, else clears the selection.
- Enter submits the focused form.
- One map, shown in the palette next to each command and in tooltips:

| Keys | Action |
|---|---|
| Ctrl+O / Ctrl+S / Ctrl+Shift+S | Open / Save session / Export image |
| Ctrl+Z / Ctrl+Shift+Z (Ctrl+Y) | Undo / Redo |
| Ctrl+P | Find commands |
| Ctrl+, | Preferences |
| Ctrl+1–8 | Show panel 1–8 |
| R T S C | Rotate / Translate / Scale / Centre mouse mode |
| Q | Cycle the selection shape (click, box, circle, lasso) |
| F | Reset view |
| Space | Play / pause |
| Alt or F10 | Ribbon key tips |
