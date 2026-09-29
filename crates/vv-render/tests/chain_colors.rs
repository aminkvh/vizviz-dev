//! Per-chain colouring keys on the chain name, not the chain record.

use vv_render::{colors_for, ColorScheme};

fn record(kind: &str, serial: u32, res: &str, chain: char, seq: u32) -> String {
    format!(
        "{kind:<6}{serial:>5}  CA  {res:<4}{chain}{seq:>4}    {:8.3}{:8.3}{:8.3}  1.00  0.00           C\n",
        seq as f32 * 3.8,
        0.0,
        0.0
    )
}

#[test]
fn ter_split_chain_records_share_a_chain_colour() {
    let text = record("ATOM", 1, "ALA", 'A', 1)
        + "TER\n"
        + &record("HETATM", 3, "HOH", 'A', 2)
        + &record("ATOM", 4, "GLY", 'B', 1);
    let s = vv_io::pdb::parse(text.as_bytes()).unwrap();
    assert_eq!(s.topology.chain_count(), 3);
    for scheme in [ColorScheme::Chain, ColorScheme::Hetero] {
        let c = colors_for(scheme, &s.topology);
        assert_eq!(c[0], c[1], "{scheme:?}: both records of chain A");
        assert_ne!(c[0], c[2], "{scheme:?}: chain B differs");
    }
}
