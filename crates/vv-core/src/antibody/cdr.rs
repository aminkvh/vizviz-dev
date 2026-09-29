//! CDR boundary definitions, each stated in the numbering it was published
//! in (docs/ANTIBODY.md, "CDR definitions").

use super::numbering::{Label, Scheme};
use super::ChainType;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CdrDefinition {
    Kabat,
    Chothia,
    Imgt,
    /// MacCallum et al. 1996, antigen-contact residues.
    Contact,
    /// North et al. 2011.
    North,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Region {
    Fr1,
    Cdr1,
    Fr2,
    Cdr2,
    Fr3,
    Cdr3,
    Fr4,
}

type Range = (Label, Label);

fn num(n: u16) -> Label {
    Label::new(n)
}

fn ins(n: u16, letter: char) -> Label {
    Label::with_insertion(n, letter)
}

impl CdrDefinition {
    /// Numbering in which the definition's boundaries are stated.
    pub fn native_scheme(self) -> Scheme {
        match self {
            CdrDefinition::Kabat | CdrDefinition::Contact => Scheme::Kabat,
            CdrDefinition::Chothia | CdrDefinition::North => Scheme::Chothia,
            CdrDefinition::Imgt => Scheme::Imgt,
        }
    }

    /// Inclusive CDR1, CDR2, CDR3 ranges in the native numbering.
    fn ranges(self, chain: ChainType) -> [Range; 3] {
        let heavy = chain == ChainType::Heavy;
        let r = |a, b| (num(a), num(b));
        match (self, heavy) {
            (CdrDefinition::Kabat, true) => [(num(31), ins(35, 'B')), r(50, 65), r(95, 102)],
            (CdrDefinition::Chothia, true) => [r(26, 32), r(52, 56), r(95, 102)],
            (CdrDefinition::North, true) => [r(23, 35), r(50, 58), r(95, 102)],
            (CdrDefinition::Contact, true) => [(num(30), ins(35, 'B')), r(47, 58), r(93, 101)],
            (CdrDefinition::Contact, false) => [r(30, 36), r(46, 55), r(89, 96)],
            (CdrDefinition::Kabat | CdrDefinition::Chothia | CdrDefinition::North, false) => {
                [r(24, 34), r(50, 56), r(89, 97)]
            }
            (CdrDefinition::Imgt, _) => [
                (num(27), ins(38, '~')),
                (num(56), ins(65, '~')),
                (num(105), ins(117, '~')),
            ],
        }
    }

    /// Region of each residue given its labels in [`Self::native_scheme`].
    pub(super) fn regions(self, chain: ChainType, native: &[Label]) -> Vec<Region> {
        const ORDER: [(Region, Region); 3] = [
            (Region::Fr1, Region::Cdr1),
            (Region::Fr2, Region::Cdr2),
            (Region::Fr3, Region::Cdr3),
        ];
        let ranges = self.ranges(chain);
        let mut stage = 0;
        native
            .iter()
            .map(|&label| {
                while stage < 3 && past(ranges[stage], label) {
                    stage += 1;
                }
                match stage {
                    3 => Region::Fr4,
                    s if inside(ranges[s], label) => ORDER[s].1,
                    s => ORDER[s].0,
                }
            })
            .collect()
    }
}

fn inside((lo, hi): Range, label: Label) -> bool {
    lo <= label && label <= hi
}

fn past((_, hi): Range, label: Label) -> bool {
    label > hi
}
