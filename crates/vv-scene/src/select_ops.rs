//! `select expand|grow|shrink|interface`: selection edits that start from
//! the active selection (or, for `interface`, from two expressions) and
//! each run as one undoable `Select` or `SelectExpr`.

use std::sync::Arc;

use vv_core::fixedbitset::FixedBitSet;
use vv_core::Topology;

use crate::command::Command;
use crate::history::CommandHistory;
use crate::scene::{LoadedStructure, Scene};
use crate::script::ScriptError;

/// Contact distance in angstroms when `interface` is given none.
const DEFAULT_REACH: f32 = 5.0;

#[derive(Debug, Clone, Copy, PartialEq)]
enum Reshape {
    Residues,
    Chains,
    Molecules,
    Within(f32),
    Grow(usize),
    Shrink(usize),
    ShrinkWithin(f32),
}

fn bad(usage: &str) -> ScriptError {
    ScriptError(format!("usage: select {usage}"))
}

fn radius(word: Option<&str>, usage: &str) -> Result<f32, ScriptError> {
    word.and_then(|w| w.parse::<f32>().ok())
        .filter(|r| r.is_finite() && *r > 0.0)
        .ok_or_else(|| bad(usage))
}

fn steps(word: Option<&str>, usage: &str) -> Result<usize, ScriptError> {
    word.map_or(Ok(1), |w| {
        w.parse().ok().filter(|&n| n > 0).ok_or_else(|| bad(usage))
    })
}

fn parse_reshape(verb: &str, args: &str) -> Result<Option<Reshape>, ScriptError> {
    let mut words = args.split_whitespace();
    let (first, second) = (words.next(), words.next());
    let usage = match verb {
        "expand" => "expand residue|chain|molecule|within N",
        "grow" => "grow [N]",
        "shrink" => "shrink [N] | shrink within N",
        _ => return Ok(None),
    };
    if words.next().is_some() {
        return Err(bad(usage));
    }
    let reshape = match (verb, first) {
        ("expand", Some("residue")) => Reshape::Residues,
        ("expand", Some("chain")) => Reshape::Chains,
        ("expand", Some("molecule")) => Reshape::Molecules,
        ("expand", Some("within")) => Reshape::Within(radius(second, usage)?),
        ("grow", n) if second.is_none() => Reshape::Grow(steps(n, usage)?),
        ("shrink", Some("within")) => Reshape::ShrinkWithin(radius(second, usage)?),
        ("shrink", n) if second.is_none() => Reshape::Shrink(steps(n, usage)?),
        _ => return Err(bad(usage)),
    };
    Ok(Some(reshape))
}

fn reshaped(loaded: &LoadedStructure, mask: &FixedBitSet, how: Reshape) -> FixedBitSet {
    let topology = &loaded.structure.topology;
    let coords = loaded.structure.frame(0);
    let positions = coords.positions();
    match how {
        Reshape::Residues => vv_core::select::byres(topology, mask),
        Reshape::Chains => whole_chains(topology, mask),
        Reshape::Molecules => whole_fragments(mask, &loaded.bonds.fragments(positions.len())),
        Reshape::Within(r) => vv_core::select::within(topology, positions, r, mask),
        Reshape::Grow(n) => (0..n).fold(vv_core::select::byres(topology, mask), |m, _| {
            step_along_chains(topology, &m, true)
        }),
        Reshape::Shrink(n) => (0..n).fold(vv_core::select::byres(topology, mask), |m, _| {
            step_along_chains(topology, &m, false)
        }),
        Reshape::ShrinkWithin(r) => {
            let mut outside = mask.clone();
            outside.toggle_range(..);
            let mut edge = vv_core::select::within(topology, positions, r, &outside);
            edge.intersect_with(mask);
            let mut kept = mask.clone();
            kept.difference_with(&edge);
            kept
        }
    }
}

fn whole_chains(topology: &Topology, mask: &FixedBitSet) -> FixedBitSet {
    let mut hit = vec![false; topology.chains.len()];
    for atom in mask.ones() {
        hit[topology.chain_of_atom(atom) as usize] = true;
    }
    fill_residues(topology, |_, res| hit[res.chain as usize])
}

fn whole_fragments(mask: &FixedBitSet, fragments: &[u32]) -> FixedBitSet {
    let mut hit = vec![false; fragments.iter().map(|&f| f as usize + 1).max().unwrap_or(0)];
    for atom in mask.ones() {
        hit[fragments[atom] as usize] = true;
    }
    let mut out = FixedBitSet::with_capacity(fragments.len());
    out.extend((0..fragments.len()).filter(|&a| hit[fragments[a] as usize]));
    out
}

fn fill_residues(
    topology: &Topology,
    pick: impl Fn(usize, &vv_core::ResidueRec) -> bool,
) -> FixedBitSet {
    let mut out = FixedBitSet::with_capacity(topology.atom_count());
    for (r, res) in topology.residues.iter().enumerate() {
        if pick(r, res) {
            out.insert_range(res.atoms.start as usize..res.atoms.end as usize);
        }
    }
    out
}

/// One residue step along each chain: growing adds the sequence neighbours
/// of every selected residue; shrinking keeps only residues whose two
/// neighbours are selected too, so each run loses a residue at both ends.
fn step_along_chains(topology: &Topology, mask: &FixedBitSet, grow: bool) -> FixedBitSet {
    let on: Vec<bool> = topology
        .residues
        .iter()
        .map(|res| mask.contains(res.atoms.start as usize))
        .collect();
    let neighbour_on = |r: usize, chain: u32, offset: isize| {
        r.checked_add_signed(offset)
            .and_then(|n| topology.residues.get(n).map(|res| (n, res)))
            .is_some_and(|(n, res)| res.chain == chain && on[n])
    };
    fill_residues(topology, |r, res| {
        let before = neighbour_on(r, res.chain, -1);
        let after = neighbour_on(r, res.chain, 1);
        if grow {
            on[r] || before || after
        } else {
            on[r] && before && after
        }
    })
}

/// The expression for `interface A [to B | with B] [within N]`: residues
/// of A near B, plus B's near A for `with`; B defaults to every polymer
/// atom outside A.
fn interface_expr(args: &str) -> Result<String, ScriptError> {
    const USAGE: &str = "interface A [to|with B] [within N]";
    let args = args.trim();
    let (args, reach) = match args.rsplit_once(" within ") {
        Some((head, n)) if n.trim().parse::<f32>().is_ok() => {
            (head, radius(Some(n.trim()), USAGE)?)
        }
        _ => (args, DEFAULT_REACH),
    };
    let split = [" to ", " with "]
        .into_iter()
        .filter_map(|sep| args.find(sep).map(|at| (at, sep)))
        .min();
    let (a, b, both) = match split {
        Some((at, sep)) => (&args[..at], args[at + sep.len()..].trim(), sep == " with "),
        None => (args, "", false),
    };
    let a = a.trim();
    if a.is_empty() {
        return Err(bad(USAGE));
    }
    let b = if b.is_empty() {
        format!("polymer and not ({a})")
    } else {
        b.to_owned()
    };
    let side = |x: &str, y: &str| format!("byres (({x}) and within {reach} of ({y}))");
    Ok(if both {
        format!("{} or {}", side(a, &b), side(&b, a))
    } else {
        side(a, &b)
    })
}

fn selected_count(scene: &Scene) -> usize {
    scene
        .active_selection()
        .map_or(0, |a| a.mask.count_ones(..))
}

/// Runs `rest` (the words after `select`) when it is one of this module's
/// forms; `None` leaves it to the expression forms.
pub(crate) fn run(
    scene: &mut Scene,
    history: &mut CommandHistory,
    rest: &str,
) -> Result<Option<String>, ScriptError> {
    let (verb, args) = rest.split_once(char::is_whitespace).unwrap_or((rest, ""));
    if verb == "interface" {
        let id = crate::script::current(scene)?;
        let expr = interface_expr(args)?;
        history.dispatch(scene, Command::SelectExpr { id, expr })?;
    } else if let Some(how) = parse_reshape(verb, args)? {
        let active = scene
            .active_selection()
            .filter(|a| !a.mask.is_clear())
            .ok_or_else(|| ScriptError("nothing is selected".into()))?;
        let (id, mask) = (active.structure, Arc::clone(&active.mask));
        let loaded = scene
            .structure(id)
            .ok_or_else(|| ScriptError("no structure loaded".into()))?;
        let out = Arc::new(reshaped(loaded, &mask, how));
        history.dispatch(scene, Command::Select { id, mask: out })?;
    } else {
        return Ok(None);
    }
    Ok(Some(format!("selected {} atom(s)", selected_count(scene))))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::script::run_line;

    fn loaded() -> (Scene, CommandHistory) {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(50);
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/1AKE.pdb");
        history
            .dispatch(&mut scene, Command::LoadStructure { path })
            .unwrap();
        (scene, history)
    }

    fn run_ok(scene: &mut Scene, history: &mut CommandHistory, line: &str) -> usize {
        run_line(scene, history, line).unwrap_or_else(|e| panic!("{line}: {e}"));
        selected_count(scene)
    }

    #[test]
    fn expand_grows_to_residue_chain_and_shell_and_undoes() {
        let (mut scene, mut history) = loaded();
        assert_eq!(run_ok(&mut scene, &mut history, "select index 0"), 1);
        let residue = run_ok(&mut scene, &mut history, "select expand residue");
        assert!(residue > 1);
        let shell = run_ok(&mut scene, &mut history, "select expand within 4");
        assert!(shell > residue);
        let chain = run_ok(&mut scene, &mut history, "select chain A");
        run_ok(&mut scene, &mut history, "select index 0");
        let whole = run_ok(&mut scene, &mut history, "select expand chain");
        assert!(
            whole > residue && whole <= chain,
            "one chain record, not its ligands"
        );
        history.undo(&mut scene).unwrap();
        assert_eq!(selected_count(&scene), 1);
    }

    #[test]
    fn grow_then_shrink_by_a_residue_round_trips_an_interior_run() {
        let (mut scene, mut history) = loaded();
        let base = run_ok(&mut scene, &mut history, "select chain A and resid 10-20");
        assert!(run_ok(&mut scene, &mut history, "select grow") > base);
        assert_eq!(run_ok(&mut scene, &mut history, "select shrink"), base);
        assert!(run_ok(&mut scene, &mut history, "select shrink within 3") < base);
    }

    #[test]
    fn invert_complements_and_undoes() {
        let (mut scene, mut history) = loaded();
        let a = run_ok(&mut scene, &mut history, "select chain A");
        let total = scene.structures().next().unwrap().1.structure.atom_count();
        assert_eq!(run_ok(&mut scene, &mut history, "select invert"), total - a);
        history.undo(&mut scene).unwrap();
        assert_eq!(selected_count(&scene), a);
    }

    #[test]
    fn remove_takes_an_expression_out_of_the_selection() {
        let (mut scene, mut history) = loaded();
        let all = run_ok(&mut scene, &mut history, "select chain A");
        let water = run_ok(&mut scene, &mut history, "select chain A and water");
        run_ok(&mut scene, &mut history, "select chain A");
        let rest = run_ok(&mut scene, &mut history, "select remove water");
        assert_eq!(rest, all - water);
    }

    #[test]
    fn interface_selects_residues_near_the_other_chain() {
        let (mut scene, mut history) = loaded();
        let chain_a = run_ok(&mut scene, &mut history, "select chain A");
        let one_side = run_ok(
            &mut scene,
            &mut history,
            "select interface chain A to chain B within 5",
        );
        assert!(one_side > 0 && one_side < chain_a);
        let both = run_ok(
            &mut scene,
            &mut history,
            "select interface chain A with chain B within 5",
        );
        assert!(both > one_side);
        assert!(run_ok(&mut scene, &mut history, "select interface chain A") > 0);
        history.undo(&mut scene).unwrap();
        assert_eq!(selected_count(&scene), both);
    }

    #[test]
    fn bad_forms_are_usage_errors() {
        let (mut scene, mut history) = loaded();
        assert!(run_line(&mut scene, &mut history, "select expand residue").is_err());
        run_ok(&mut scene, &mut history, "select chain A");
        for line in [
            "select expand",
            "select expand within",
            "select expand within -1",
            "select grow 0",
            "select shrink x",
            "select interface",
        ] {
            assert!(run_line(&mut scene, &mut history, line).is_err(), "{line}");
        }
    }

    #[test]
    fn interface_expression_text() {
        assert_eq!(
            interface_expr("chain A to chain B within 4").unwrap(),
            "byres ((chain A) and within 4 of (chain B))"
        );
        assert!(interface_expr("polymer with ligand")
            .unwrap()
            .contains(" or "));
    }
}
