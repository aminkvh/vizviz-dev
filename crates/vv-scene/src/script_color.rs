//! The `color` command: a scheme (with a property's scale, range and
//! ramp), a solid color, or per-target overrides over either.

use crate::coloring::{ColorOverride, ColorTarget, PropertyKind};
use crate::command::Command;
use crate::history::CommandHistory;
use crate::scene::{parse_solid, ColorScheme, Scene, StructureId};
use crate::script::{current_rep, resolve_structure, usage, ScriptError};

pub(crate) fn run_color(
    scene: &mut Scene,
    history: &mut CommandHistory,
    rest: &str,
) -> Result<String, ScriptError> {
    let words: Vec<&str> = rest.split_whitespace().collect();
    match words.first().copied() {
        None => Err(usage("color")),
        Some("set") => set_override(scene, history, &words[1..]),
        Some("unset") => unset_override(scene, history, &words[1..]),
        Some("overrides") => list_overrides(scene, words.get(1).copied().unwrap_or("")),
        Some(_) => set_scheme(scene, history, &words),
    }
}

fn set_scheme(
    scene: &mut Scene,
    history: &mut CommandHistory,
    words: &[&str],
) -> Result<String, ScriptError> {
    let (coloring, used) = ColorScheme::parse_words(words).ok_or_else(|| {
        ScriptError(format!(
            "unknown coloring or option in `{}`; see `help color`",
            words.join(" ")
        ))
    })?;
    let id_word = trailing_id(&words[used..])?;
    let id = resolve_structure(scene, id_word)?;
    let missing = match &coloring {
        ColorScheme::Property(p) => match &p.kind {
            PropertyKind::Values(name) => !scene
                .structure(id)
                .expect("resolved")
                .values
                .contains_key(name),
            _ => false,
        },
        _ => false,
    };
    let label = coloring.name();
    let rep = current_rep(scene, id);
    history.dispatch(scene, Command::SetColoring { id, rep, coloring })?;
    Ok(if missing {
        format!(
            "#{} colored by {label} (no such channel yet; drawn by element until `values` attaches it)",
            id.to_raw()
        )
    } else {
        format!("#{} colored by {label}", id.to_raw())
    })
}

/// The structure id word after a command's other words: none, or one.
fn trailing_id<'a>(words: &[&'a str]) -> Result<&'a str, ScriptError> {
    match words {
        [] => Ok(""),
        [id] if id.parse::<u32>().is_ok() => Ok(id),
        [word, ..] if matches!(*word, "scale" | "range" | "ramp") => Err(ScriptError(format!(
            "`{}` here takes no scale, range or ramp, or its value is missing",
            words.join(" ")
        ))),
        _ => Err(usage("color")),
    }
}

/// Where a `sel` target's expression ends in `args`: before a trailing
/// `COLOR [ID]` (`with_color`), or before a trailing structure id that
/// leaves a valid expression behind.
fn sel_end(args: &[&str], with_color: bool) -> usize {
    let n = args.len();
    let is_color = |i: usize| args.get(i).is_some_and(|w| parse_solid(w).is_some());
    if with_color {
        return match n {
            0 => 0,
            _ if is_color(n - 1) => n - 1,
            _ if n >= 2 && is_color(n - 2) => n - 2,
            _ => n,
        };
    }
    let id_last = n >= 2 && args[n - 1].parse::<u32>().is_ok();
    if id_last && vv_core::select::parse(&args[..n - 1].join(" ")).is_ok() {
        n - 1
    } else {
        n
    }
}

/// The target at the front of `words` (`KIND ARG`, or `sel EXPR...`) and
/// the words after it.
fn split_target<'a>(
    words: &'a [&'a str],
    with_color: bool,
) -> Result<(ColorTarget, &'a [&'a str]), ScriptError> {
    let [kind, args @ ..] = words else {
        return Err(usage("color"));
    };
    let (arg, tail) = match (*kind, args) {
        ("sel", _) => {
            let end = sel_end(args, with_color);
            (args[..end].join(" "), &args[end..])
        }
        (_, [arg, tail @ ..]) => ((*arg).to_string(), tail),
        (_, []) => (String::new(), args),
    };
    let target = ColorTarget::parse(kind, &arg).map_err(ScriptError)?;
    Ok((target, tail))
}

fn hex(rgb: [u8; 3]) -> String {
    format!("#{:02x}{:02x}{:02x}", rgb[0], rgb[1], rgb[2])
}

fn set_override(
    scene: &mut Scene,
    history: &mut CommandHistory,
    words: &[&str],
) -> Result<String, ScriptError> {
    let (target, tail) = split_target(words, true)?;
    let (color, id_word) = match tail {
        [color] => (*color, ""),
        [color, id] => (*color, *id),
        _ => return Err(usage("color")),
    };
    let color = parse_color(color)?;
    let id = resolve_structure(scene, id_word)?;
    let atoms = matched_atoms(scene, id, &target)?;
    let mut overrides = scene
        .structure(id)
        .expect("resolved")
        .color_overrides
        .clone();
    overrides.retain(|o| o.target != target);
    overrides.push(ColorOverride {
        target: target.clone(),
        color,
    });
    history.dispatch(scene, Command::SetColorOverrides { id, overrides })?;
    Ok(format!(
        "#{} {target} -> {} ({atoms} atoms)",
        id.to_raw(),
        hex(color)
    ))
}

fn parse_color(word: &str) -> Result<[u8; 3], ScriptError> {
    parse_solid(word)
        .ok_or_else(|| ScriptError(format!("`{word}` is not a color; use #RRGGBB or a name")))
}

fn matched_atoms(
    scene: &Scene,
    id: StructureId,
    target: &ColorTarget,
) -> Result<usize, ScriptError> {
    let loaded = scene.structure(id).expect("resolved");
    let bits = loaded
        .select(&target.expression(), loaded.frame)
        .map_err(|e| ScriptError(format!("{target}: {e}")))?;
    Ok(bits.count_ones(..))
}

fn unset_override(
    scene: &mut Scene,
    history: &mut CommandHistory,
    words: &[&str],
) -> Result<String, ScriptError> {
    if let ["all", tail @ ..] = words {
        let id = resolve_structure(scene, trailing_id(tail)?)?;
        history.dispatch(
            scene,
            Command::SetColorOverrides {
                id,
                overrides: Vec::new(),
            },
        )?;
        return Ok(format!("#{} color overrides cleared", id.to_raw()));
    }
    let (target, tail) = split_target(words, false)?;
    let id = resolve_structure(scene, trailing_id(tail)?)?;
    let mut overrides = scene
        .structure(id)
        .expect("resolved")
        .color_overrides
        .clone();
    let before = overrides.len();
    overrides.retain(|o| o.target != target);
    if overrides.len() == before {
        return Err(ScriptError(format!("no override for {target}")));
    }
    history.dispatch(scene, Command::SetColorOverrides { id, overrides })?;
    Ok(format!("#{} {target} override removed", id.to_raw()))
}

fn list_overrides(scene: &Scene, id_word: &str) -> Result<String, ScriptError> {
    let id = resolve_structure(scene, id_word)?;
    let overrides = &scene.structure(id).expect("resolved").color_overrides;
    if overrides.is_empty() {
        return Ok(format!("#{}: no color overrides", id.to_raw()));
    }
    let lines: Vec<String> = overrides
        .iter()
        .map(|o| format!("{} {}", o.target, hex(o.color)))
        .collect();
    Ok(lines.join("\n"))
}
