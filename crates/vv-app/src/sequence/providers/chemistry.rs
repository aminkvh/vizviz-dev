//! Tracks about what the residues are and what is bonded to them:
//! disulfides, glycans, sequence liabilities, modified residues, and
//! alternate locations.

use std::collections::HashSet;
use std::ops::Range;

use vv_core::residue_class::is_standard_polymer;
use vv_core::seqfeat::{self, Liability, SequonKind};

use super::hex;
use crate::sequence::rows::{describe_residue, is_protein, letters_and_breaks};
use crate::sequence::tracks::{legend, Glyph, TrackContext, TrackData, TrackProvider};

const PAIR_COLORS: [u32; 6] = [0xE67E22, 0x3498DB, 0x2ECC71, 0xE91E63, 0x9B59B6, 0x1ABC9C];

pub struct Disulfides;

impl TrackProvider for Disulfides {
    fn id(&self) -> &'static str {
        "disulfide"
    }

    fn label(&self) -> &'static str {
        "Disulfide"
    }

    /// Bonded cysteines share a color; the legend shows the first.
    fn compute(&self, ctx: &TrackContext) -> TrackData {
        let top = ctx.top();
        let entries = PAIR_COLORS.iter().map(|&c| legend("", hex(c))).collect();
        let mut track = ctx.new_track(Glyph::Bar, entries);
        track.legend = vec![legend(
            "Disulfide bridge (partners share a color)",
            hex(PAIR_COLORS[0]),
        )];
        for (i, [a, b]) in seqfeat::disulfides(top, &ctx.loaded.bonds)
            .into_iter()
            .enumerate()
        {
            let kind = (i % PAIR_COLORS.len()) as u8 + 1;
            for (here, there) in [(a, b), (b, a)] {
                track.mark(here, kind);
                track.note(
                    here,
                    format!("Disulfide to {}", describe_residue(top, there)),
                );
            }
        }
        track.finish()
    }
}

pub struct Glycans;

const SEQUON: u8 = 1;
const RARE_SEQUON: u8 = 2;
const N_GLYCAN: u8 = 3;
const O_GLYCAN: u8 = 4;

impl Glycans {
    fn mark_sequons(track: &mut TrackData, ctx: &TrackContext, residues: Range<u32>) {
        let top = ctx.top();
        let (letters, breaks) = letters_and_breaks(top, residues.clone());
        for s in seqfeat::sequons(&letters, &breaks) {
            let r = residues.start + s.asn as u32;
            let (kind, pattern) = match s.kind {
                SequonKind::Canonical => (SEQUON, "N-X-S/T"),
                SequonKind::Rare => (RARE_SEQUON, "N-X-C"),
            };
            let window = String::from_utf8_lossy(&letters[s.asn..s.asn + 3]).into_owned();
            track.mark(r, kind);
            track.note(r, format!("{pattern} sequon ({window}), no glycan modeled"));
        }
    }

    fn mark_attached(track: &mut TrackData, ctx: &TrackContext) {
        let top = ctx.top();
        let mut noted: HashSet<u32> = HashSet::new();
        for g in seqfeat::glycosylated(top, &ctx.loaded.bonds) {
            let n_linked = top.atom_name(g.atom as usize) == "ND2";
            track.mark(g.residue, if n_linked { N_GLYCAN } else { O_GLYCAN });
            if noted.insert(g.residue) {
                track.note(
                    g.residue,
                    if n_linked {
                        "N-glycosylated"
                    } else {
                        "Glycosylated"
                    },
                );
            }
            track.note(g.residue, format!("  {}", describe_residue(top, g.glycan)));
        }
    }
}

impl TrackProvider for Glycans {
    fn id(&self) -> &'static str {
        "glycan"
    }

    fn label(&self) -> &'static str {
        "Glycosylation"
    }

    fn compute(&self, ctx: &TrackContext) -> TrackData {
        let top = ctx.top();
        let mut track = ctx.new_track(
            Glyph::Bar,
            vec![
                legend("N-X-S/T sequon, no glycan modeled", hex(0x85C1E9)),
                legend("N-X-C sequon (rare), no glycan modeled", hex(0xD6EAF8)),
                legend("Glycosylated Asn", hex(0x1F618D)),
                legend("Glycosylated Ser, Thr or other", hex(0x117A65)),
            ],
        );
        for (_, residues) in ctx.rows {
            if residues.clone().any(|r| is_protein(top, r)) {
                Self::mark_sequons(&mut track, ctx, residues.clone());
            }
        }
        Self::mark_attached(&mut track, ctx);
        track.finish()
    }
}

pub struct Liabilities;

const LIABILITY_COLORS: [u32; 7] = [
    0xE67E22, 0x9B59B6, 0xA93226, 0xE91E8C, 0x27AE60, 0x16A085, 0xF1C40F,
];

fn liability_kind(kind: Liability) -> u8 {
    Liability::ALL.iter().position(|&k| k == kind).unwrap() as u8 + 1
}

impl Liabilities {
    fn mark_row(
        track: &mut TrackData,
        ctx: &TrackContext,
        residues: Range<u32>,
        paired: &HashSet<u32>,
    ) {
        let (letters, breaks) = letters_and_breaks(ctx.top(), residues.clone());
        for hit in seqfeat::liabilities(&letters, &breaks) {
            let first = residues.start + hit.start as u32;
            if hit.kind == Liability::FreeCysteine && paired.contains(&first) {
                continue;
            }
            let kind = liability_kind(hit.kind);
            for r in first..first + hit.len as u32 {
                if track.kind(r) == 0 || kind < track.kind(r) {
                    track.mark(r, kind);
                }
                track.note(r, format!("{}: {}", hit.motif, hit.kind.label()));
            }
        }
    }
}

impl TrackProvider for Liabilities {
    fn id(&self) -> &'static str {
        "liability"
    }

    fn label(&self) -> &'static str {
        "Liabilities"
    }

    fn compute(&self, ctx: &TrackContext) -> TrackData {
        let top = ctx.top();
        let paired: HashSet<u32> = seqfeat::disulfides(top, &ctx.loaded.bonds)
            .into_iter()
            .flatten()
            .collect();
        let entries = Liability::ALL
            .iter()
            .zip(LIABILITY_COLORS)
            .map(|(k, c)| legend(&liability_legend(*k), hex(c)))
            .collect();
        let mut track = ctx.new_track(Glyph::Bar, entries);
        for (_, residues) in ctx.rows {
            if residues.clone().any(|r| is_protein(top, r)) {
                Self::mark_row(&mut track, ctx, residues.clone(), &paired);
            }
        }
        track.finish()
    }
}

fn liability_legend(kind: Liability) -> String {
    let motifs = match kind {
        Liability::Deamidation => "NG NS NT",
        Liability::Isomerization => "DG DS DT",
        Liability::Fragmentation => "DP",
        Liability::Integrin => "RGD",
        Liability::FreeCysteine => "unpaired Cys",
        Liability::PyroGlutamate => "N-terminal Q/E",
        Liability::Oxidation => "Met, Trp",
    };
    format!("{} ({motifs})", kind.label())
}

pub struct Modified;

impl TrackProvider for Modified {
    fn id(&self) -> &'static str {
        "modified"
    }

    fn label(&self) -> &'static str {
        "Modified"
    }

    fn compute(&self, ctx: &TrackContext) -> TrackData {
        let top = ctx.top();
        let mut track = ctx.new_track(Glyph::Bar, vec![legend("Modified residue", hex(0xE84393))]);
        for r in 0..top.residue_count() as u32 {
            let name = top.residue_name(r as usize);
            if is_protein(top, r) && !is_standard_polymer(name) {
                track.mark(r, 1);
                track.note(r, format!("Modified residue {name}"));
            }
        }
        track.finish()
    }
}

pub struct AltLocs;

impl TrackProvider for AltLocs {
    fn id(&self) -> &'static str {
        "altloc"
    }

    fn label(&self) -> &'static str {
        "Alt. locations"
    }

    fn compute(&self, ctx: &TrackContext) -> TrackData {
        let top = ctx.top();
        let mut track = ctx.new_track(
            Glyph::Bar,
            vec![legend("Alternate locations", hex(0x7F8C8D))],
        );
        for (r, rec) in top.residues.iter().enumerate() {
            let mut labels: Vec<u8> = rec
                .atoms
                .clone()
                .map(|a| top.alt_loc[a as usize])
                .filter(|&l| l != 0)
                .collect();
            labels.sort_unstable();
            labels.dedup();
            if !labels.is_empty() {
                let names: Vec<String> = labels.iter().map(|&l| (l as char).to_string()).collect();
                track.mark(r as u32, 1);
                track.note(
                    r as u32,
                    format!("Alternate locations: {}", names.join(", ")),
                );
            }
        }
        track.finish()
    }
}
