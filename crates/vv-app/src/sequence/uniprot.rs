//! UniProt features for a structure's chains: SIFTS says which UniProt
//! residues a chain's residues are, the UniProt record says what is
//! annotated there. Parsing and mapping are pure; only [`load`] uses the
//! network (through the on-disk cache of `vv_io::seqdata::uniprot`).

use std::collections::HashMap;
use std::ops::Range;
use std::path::Path;

use serde_json::Value;
use vv_core::Topology;

/// One aligned stretch of a chain against a UniProt entry.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Segment {
    /// Chain id as the structure names it (see [`parse_sifts`]).
    pub chain: String,
    pub accession: String,
    pub name: String,
    pub unp: (u32, u32),
    /// The chain's residue numbers over the same stretch, in the
    /// numbering the structure's residues carry.
    pub numbers: (i32, i32),
}

impl Segment {
    /// The UniProt position of the residue numbered `seq`, when the
    /// stretch is one-to-one.
    fn position(&self, seq: i32) -> Option<u32> {
        let linear = i64::from(self.unp.1) - i64::from(self.unp.0)
            == i64::from(self.numbers.1) - i64::from(self.numbers.0);
        let inside = (self.numbers.0..=self.numbers.1).contains(&seq);
        (linear && inside).then(|| self.unp.0 + (seq - self.numbers.0) as u32)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Feature {
    pub kind: Kind,
    /// UniProt type name, as UniProt writes it.
    pub label: String,
    pub start: u32,
    pub end: u32,
    pub text: String,
}

/// How a feature is drawn; later kinds are painted over earlier ones.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Kind {
    Region,
    Topology,
    Variant,
    /// UniProt's generic "Site" (cleavage, glycation, ...): kept apart so
    /// red means an active or binding site only.
    OtherSite,
    Modification,
    Site,
}

impl Kind {
    pub fn legend(self) -> &'static str {
        match self {
            Kind::Region => "Domain, region, repeat, motif",
            Kind::Topology => "Signal, transmembrane, topology",
            Kind::Variant => "Natural variant",
            Kind::OtherSite => "Other site",
            Kind::Modification => "Modification",
            Kind::Site => "Active or binding site",
        }
    }

    fn of(uniprot_type: &str) -> Option<Kind> {
        Some(match uniprot_type {
            "Domain" | "Region" | "Repeat" | "Motif" | "Zinc finger" | "DNA binding"
            | "Coiled coil" | "Compositional bias" => Kind::Region,
            "Signal" | "Transit peptide" | "Propeptide" | "Transmembrane" | "Intramembrane"
            | "Topological domain" => Kind::Topology,
            "Natural variant" => Kind::Variant,
            "Modified residue" | "Glycosylation" | "Disulfide bond" | "Lipidation"
            | "Cross-link" => Kind::Modification,
            "Active site" | "Binding site" => Kind::Site,
            "Site" => Kind::OtherSite,
            _ => return None,
        })
    }
}

/// Everything fetched for one structure.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Data {
    pub segments: Vec<Segment>,
    /// Features by accession.
    pub features: HashMap<String, Vec<Feature>>,
}

fn number(v: &Value) -> Option<u32> {
    u32::try_from(v.as_u64()?).ok()
}

/// The chain-to-UniProt segments of a SIFTS mapping. Chains are named by
/// `struct_asym_id` (the mmCIF `label_asym_id`) when `by_label`, else by
/// the author chain id; residue numbers likewise.
pub fn parse_sifts(json: &str, by_label: bool) -> Vec<Segment> {
    let Ok(root) = serde_json::from_str::<Value>(json) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    let Some(entry) = root.as_object().and_then(|o| o.values().next()) else {
        return out;
    };
    let Some(entries) = entry["UniProt"].as_object() else {
        return out;
    };
    for (accession, record) in entries {
        let name = record["identifier"].as_str().unwrap_or(accession);
        for m in record["mappings"].as_array().into_iter().flatten() {
            if let Some(seg) = segment(m, accession, name, by_label) {
                out.push(seg);
            }
        }
    }
    out
}

fn segment(m: &Value, accession: &str, name: &str, by_label: bool) -> Option<Segment> {
    let (chain_key, number_key) = if by_label {
        ("struct_asym_id", "residue_number")
    } else {
        ("chain_id", "author_residue_number")
    };
    let seq = |end: &str| m[end][number_key].as_i64().map(|n| n as i32);
    Some(Segment {
        chain: m[chain_key].as_str()?.to_string(),
        accession: accession.to_string(),
        name: name.to_string(),
        unp: (number(&m["unp_start"])?, number(&m["unp_end"])?),
        numbers: (seq("start")?, seq("end")?),
    })
}

/// The drawable features of a UniProt record.
pub fn parse_features(json: &str) -> Vec<Feature> {
    let Ok(root) = serde_json::from_str::<Value>(json) else {
        return Vec::new();
    };
    root["features"]
        .as_array()
        .into_iter()
        .flatten()
        .filter_map(feature)
        .collect()
}

fn feature(f: &Value) -> Option<Feature> {
    let label = f["type"].as_str()?;
    Some(Feature {
        kind: Kind::of(label)?,
        label: label.to_string(),
        start: number(&f["location"]["start"]["value"])?,
        end: number(&f["location"]["end"]["value"])?,
        text: feature_text(f),
    })
}

fn feature_text(f: &Value) -> String {
    let mut parts: Vec<String> = Vec::new();
    if let Some(alt) = f["alternativeSequence"].as_object() {
        let original = alt
            .get("originalSequence")
            .and_then(Value::as_str)
            .unwrap_or("");
        let changed: Vec<&str> = alt
            .get("alternativeSequences")
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
            .filter_map(Value::as_str)
            .collect();
        if !original.is_empty() {
            parts.push(format!("{original} \u{2192} {}", changed.join("/")));
        }
    }
    if let Some(ligand) = f["ligand"]["name"].as_str() {
        parts.push(ligand.to_string());
    }
    if let Some(text) = f["description"].as_str().filter(|d| !d.is_empty()) {
        parts.push(text.to_string());
    }
    parts.join("; ")
}

/// A residue of a chain and the UniProt position it is.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Mapped {
    pub residue: u32,
    pub accession: String,
    pub position: u32,
}

/// The residues of chain `name` (`residues` in `top`) with a UniProt
/// position.
pub fn map_residues(data: &Data, name: &str, top: &Topology, residues: Range<u32>) -> Vec<Mapped> {
    let segments: Vec<&Segment> = data.segments.iter().filter(|s| s.chain == name).collect();
    residues
        .filter_map(|residue| {
            let seq = top.residues[residue as usize].seq_id;
            segments.iter().find_map(|s| {
                Some(Mapped {
                    residue,
                    accession: s.accession.clone(),
                    position: s.position(seq)?,
                })
            })
        })
        .collect()
}

/// Whether residue numbers in `path`'s structure are entity positions
/// (mmCIF) rather than author numbers (PDB format).
fn numbered_by_label(path: &Path) -> bool {
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    let name = name.strip_suffix(".gz").unwrap_or(&name);
    !(name.ends_with(".pdb") || name.ends_with(".ent"))
}

/// The PDB id of a structure: the file's own `entry.id`, else a file name
/// that is one.
pub fn pdb_id(top: &Topology, path: Option<&Path>) -> Option<String> {
    let named = top.annotations.get("entry", "id").map(str::to_string);
    let stem = || {
        let name = path?.file_name()?.to_str()?;
        let stem = name.split('.').next()?;
        vv_io::fetch::normalize_id(stem).ok()
    };
    named.or_else(stem).filter(|id| id.len() == 4)
}

/// Fetches (or reads from the cache) the mapping and the features of every
/// UniProt entry it names. `None` when the mapping is unavailable, e.g.
/// offline with nothing cached; an entry whose features fail is left out.
pub fn load(id: &str, path: &Path) -> Option<Data> {
    let by_label = numbered_by_label(path);
    let sifts = vv_io::seqdata::uniprot::fetch_sifts(id).ok()?;
    let segments = parse_sifts(&sifts, by_label);
    if segments.is_empty() {
        return None;
    }
    let mut features = HashMap::new();
    for seg in &segments {
        if features.contains_key(&seg.accession) {
            continue;
        }
        if let Ok(json) = vv_io::seqdata::uniprot::fetch_features(&seg.accession) {
            features.insert(seg.accession.clone(), parse_features(&json));
        }
    }
    Some(Data { segments, features })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str) -> String {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/uniprot")
            .join(name);
        std::fs::read_to_string(path).unwrap()
    }

    #[test]
    fn sifts_gives_one_segment_per_chain_and_entry() {
        let segs = parse_sifts(&fixture("sifts_4hhb.json"), true);
        assert_eq!(segs.len(), 4);
        let a = segs.iter().find(|s| s.chain == "A").unwrap();
        assert_eq!(
            (a.accession.as_str(), a.name.as_str()),
            ("P69905", "HBA_HUMAN")
        );
        assert_eq!((a.unp, a.numbers), ((2, 142), (1, 141)));
        let b = segs.iter().find(|s| s.chain == "D").unwrap();
        assert_eq!(b.accession, "P68871");
    }

    #[test]
    fn a_segment_maps_only_inside_a_one_to_one_stretch() {
        let seg = Segment {
            chain: "A".into(),
            accession: "P".into(),
            name: "N".into(),
            unp: (2, 142),
            numbers: (1, 141),
        };
        assert_eq!(seg.position(1), Some(2));
        assert_eq!(seg.position(141), Some(142));
        assert_eq!(seg.position(0), None);
        assert_eq!(seg.position(142), None);
        let stretched = Segment {
            numbers: (1, 150),
            ..seg
        };
        assert_eq!(stretched.position(5), None);
    }

    #[test]
    fn features_keep_kind_range_and_words() {
        let feats = parse_features(&fixture("P69905.json"));
        let heme = feats
            .iter()
            .find(|f| f.kind == Kind::Site && f.start == 88)
            .unwrap();
        assert_eq!(heme.text, "heme b; proximal binding residue");
        let variant = feats.iter().find(|f| f.kind == Kind::Variant).unwrap();
        assert!(variant.text.starts_with("V \u{2192} E"));
        assert!(feats
            .iter()
            .any(|f| f.kind == Kind::Region && (f.start, f.end) == (2, 142)));
    }

    #[test]
    fn a_variant_without_an_original_residue_does_not_break_parsing() {
        let json = r#"{"features":[{"type":"Natural variant",
            "location":{"start":{"value":5},"end":{"value":6}},
            "alternativeSequence":{},"description":"deletion"}]}"#;
        let feats = parse_features(json);
        assert_eq!(feats.len(), 1);
        assert_eq!(feats[0].text, "deletion");
    }

    #[test]
    fn unparseable_text_gives_nothing() {
        assert!(parse_sifts("not json", true).is_empty());
        assert!(parse_features("{}").is_empty());
    }

    /// Hemoglobin's proximal and distal histidines: UniProt counts the
    /// initiator Met, so its 59 and 88 are the PDB's 58 and 87 (alpha),
    /// 64 and 93 the PDB's 63 and 92 (beta).
    #[test]
    fn hemoglobin_heme_histidines_land_on_the_right_residues() {
        let top = vv_io::load(
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/4HHB.cif"),
        )
        .unwrap()
        .topology
        .clone();
        let mut data = Data {
            segments: parse_sifts(&fixture("sifts_4hhb.json"), true),
            features: HashMap::new(),
        };
        data.features
            .insert("P69905".into(), parse_features(&fixture("P69905.json")));
        data.features
            .insert("P68871".into(), parse_features(&fixture("P68871.json")));
        for (chain, accession, unp_positions) in
            [("A", "P69905", [59, 88]), ("B", "P68871", [64, 93])]
        {
            let ci = (0..top.chain_count())
                .find(|&c| top.chain_name(c) == chain)
                .unwrap();
            let range = top.chains[ci].residues.clone();
            let mapped = map_residues(&data, chain, &top, range);
            for position in unp_positions {
                let m = mapped
                    .iter()
                    .find(|m| m.position == position && m.accession == accession)
                    .unwrap();
                assert_eq!(top.residue_name(m.residue as usize), "HIS");
                let site = data.features[accession]
                    .iter()
                    .any(|f| f.kind == Kind::Site && f.start == position);
                assert!(site, "{accession} {position}");
            }
        }
    }
}
