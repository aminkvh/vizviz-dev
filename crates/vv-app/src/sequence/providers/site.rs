//! Tracks about neighbours in space: residues near a ligand, and residues
//! at a contact with another chain.

use std::collections::BTreeMap;

use vv_core::residue_class::Roles;
use vv_core::seqfeat::{residue_contacts, ResidueContact};

use super::hex;
use crate::sequence::rows::{describe_residue, is_polymer};
use crate::sequence::tracks::{legend, Glyph, TrackContext, TrackData, TrackProvider};

const LIGAND_CUTOFF: f32 = 4.0;
const INTERFACE_CUTOFF: f32 = 4.5;
const NOTED_PARTNERS: usize = 3;

pub struct LigandSite;

impl TrackProvider for LigandSite {
    fn id(&self) -> &'static str {
        "ligand"
    }

    fn label(&self) -> &'static str {
        "Ligand site"
    }

    fn per_frame(&self) -> bool {
        true
    }

    fn compute(&self, ctx: &TrackContext) -> TrackData {
        let top = ctx.top();
        let mut track = ctx.new_track(
            Glyph::Bar,
            vec![legend("Within 4 Å of a ligand", hex(0x1ABC9C))],
        );
        let ligand = |r: usize| top.residue_roles(r).contains(Roles::LIGAND);
        let polymer = |r: usize| is_polymer(top, r as u32);
        let contacts = residue_contacts(
            top,
            ctx.positions,
            LIGAND_CUTOFF,
            polymer,
            ligand,
            |_, _| true,
        );
        for group in contacts.chunk_by(|x, y| x.a == y.a) {
            let a = group[0].a;
            track.mark(a, 1);
            for c in group.iter().take(NOTED_PARTNERS) {
                let partner = describe_residue(top, c.b);
                track.note(a, format!("Contacts {partner}, {:.1} Å", c.distance));
            }
        }
        track.finish()
    }
}

pub struct Interface;

fn partner_chains(ctx: &TrackContext, group: &[ResidueContact]) -> String {
    let top = ctx.top();
    let mut per_chain: BTreeMap<&str, usize> = BTreeMap::new();
    for c in group {
        let chain = top.chain_name(top.residues[c.b as usize].chain as usize);
        *per_chain.entry(chain).or_default() += 1;
    }
    let parts: Vec<String> = per_chain
        .iter()
        .map(|(chain, n)| {
            format!(
                "chain {chain} ({n} residue{})",
                if *n == 1 { "" } else { "s" }
            )
        })
        .collect();
    format!("Interface with {}", parts.join(", "))
}

impl TrackProvider for Interface {
    fn id(&self) -> &'static str {
        "interface"
    }

    fn label(&self) -> &'static str {
        "Interface"
    }

    fn per_frame(&self) -> bool {
        true
    }

    fn compute(&self, ctx: &TrackContext) -> TrackData {
        let top = ctx.top();
        let mut track = ctx.new_track(
            Glyph::Bar,
            vec![legend("Within 4.5 Å of another chain", hex(0xF39C12))],
        );
        let polymer = |r: usize| is_polymer(top, r as u32);
        let chain_of = |r: u32| top.chain_name(top.residues[r as usize].chain as usize);
        let other_chain = |a: u32, b: u32| chain_of(a) != chain_of(b);
        let contacts = residue_contacts(
            top,
            ctx.positions,
            INTERFACE_CUTOFF,
            polymer,
            polymer,
            other_chain,
        );
        for group in contacts.chunk_by(|x, y| x.a == y.a) {
            track.mark(group[0].a, 1);
            track.note(group[0].a, partner_chains(ctx, group));
        }
        track.finish()
    }
}
