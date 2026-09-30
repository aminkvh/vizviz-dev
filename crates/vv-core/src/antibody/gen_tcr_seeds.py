#!/usr/bin/env python3
"""Regenerates `tcr_seeds.rs` from public wwPDB coordinates (CC0).

Not part of the build: vv-core ships the generated table. Run by hand:

    python3 crates/vv-core/src/antibody/gen_tcr_seeds.py

Reads the entries in `tcr_seed_ids.txt` (downloaded to the git-ignored
`fixtures/real/tcr/seeds/`) and keeps every alpha or beta chain whose
depositor numbering is IMGT-framed: Cys 23, Trp 41 and Cys 104 at those
numbers. The chain type comes from the entry's COMPND record. Only the
framework columns (FR1 1-26, FR2 39-55, FR3 66-104, FR4 118-128) are kept,
`-` marking a column with no residue. Each chain joins the family whose
empty FR3 columns it shares (`FAMILIES`, mirrored in `profile.rs`); chains
that fit none are dropped. CDRs are discarded, but their lengths are
printed, because they set the loop priors in `profile.rs`. Chains are
clustered greedily at 90% identity.

The entry list: X-ray entries whose polymer description names a T cell
receptor, from an RCSB search; every fourth entry in sorted-ID order (and
1AO7, 1MI5, 2BNR, 3HG1, 5HHO) is held out for `tests/antibody_tcr.rs` and
is absent from `tcr_seed_ids.txt`.
"""
import collections
import os
import re
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
PDB_DIR = os.path.join(REPO, "fixtures", "real", "tcr", "seeds")

FR_RANGES = ((1, 26), (39, 55), (66, 104), (118, 128))
ALPHA_EXTRA = [(84, "A"), (84, "B"), (84, "C")]
# name, chain type, empty columns, columns 84A-84C filled. Chains of one
# type differ in which FR3 columns they leave empty; the depositor
# numberings measured in the seed entries fall into these groups.
FAMILIES = [
    ("ALPHA_69_73", "alpha", [69, 70, 71, 72, 73], False),
    ("ALPHA_71_77", "alpha", [71, 72, 73, 74, 75, 76, 77], True),
    ("BETA_73_82", "beta", [73, 82], False),
    ("BETA_82", "beta", [82], False),
]
# Columns whose emptiness tells the families apart.
DECIDING = {
    "alpha": [(n, "") for n in range(69, 78)] + ALPHA_EXTRA,
    "beta": [(73, ""), (82, "")],
}
THREE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M",
}
LOOPS = {"cdr1": (27, 38), "cdr2": (56, 65), "cdr3": (105, 117)}
IDENTITY = 0.90


def fetch(pid):
    os.makedirs(PDB_DIR, exist_ok=True)
    path = os.path.join(PDB_DIR, pid + ".pdb")
    if not (os.path.exists(path) and os.path.getsize(path) > 1000):
        url = "https://files.rcsb.org/download/%s.pdb" % pid
        with open(path, "wb") as f:
            f.write(urllib.request.urlopen(url, timeout=60).read())
    return path


def kind(name):
    n = name.lower()
    for k in ("alpha", "beta"):
        if k in n or re.search(r"\btcr[ -]?%s\b" % k[0], n):
            return k
    return None


def molecules(lines):
    mol, cur = {}, None
    for line in lines:
        if line.startswith("ATOM"):
            break
        if not line.startswith("COMPND"):
            continue
        text = line[10:].strip()
        m = re.match(r"MOLECULE:\s*(.*?);?$", text)
        if m:
            cur = m.group(1)
        m = re.match(r"CHAIN:\s*(.*?);?$", text)
        if m and cur is not None:
            for c in m.group(1).replace(" ", "").split(","):
                mol[c] = cur
    return mol


def chains(path):
    lines = open(path, errors="ignore").read().split("\n")
    mol = molecules(lines)
    out = collections.defaultdict(dict)
    for line in lines:
        if line.startswith("ATOM") and line[12:16] == " CA " and line[16] in " A":
            try:
                n = int(line[22:26])
            except ValueError:
                continue
            out[line[21]][(n, line[26].strip())] = THREE.get(line[17:20], "X")
    return [(c, kind(mol.get(c, "")), res) for c, res in out.items()]


def anchored(res):
    return (res.get((23, "")) == "C" and res.get((41, "")) == "W"
            and res.get((104, "")) == "C")


def family_of(res, k):
    empty = {c for c in DECIDING[k] if c not in res}
    for name, kind_, cols, extras in FAMILIES:
        if kind_ != k:
            continue
        want = {(n, "") for n in cols} | (set() if extras else set(ALPHA_EXTRA if k == "alpha" else []))
        if empty == want:
            return name
    return None


def framework(res, name):
    _, _, cols, extras = next(f for f in FAMILIES if f[0] == name)
    blocks = []
    for lo, hi in FR_RANGES:
        keys = [(n, "") for n in range(lo, hi + 1) if n not in cols]
        if (lo, hi) == (66, 104) and extras:
            at = keys.index((84, ""))
            keys[at + 1:at + 1] = ALPHA_EXTRA
        blocks.append("".join(res.get(c, "-") for c in keys))
    return blocks


def loop_lengths(res):
    return {name: sum(1 for (n, _) in res if lo <= n <= hi) for name, (lo, hi) in LOOPS.items()}


def identity(a, b):
    pairs = [(x, y) for x, y in zip(a, b) if x != "-" and y != "-"]
    return sum(x == y for x, y in pairs) / max(1, len(pairs))


def cluster(counts):
    reps = []
    for seg, count in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
        joined = "".join(seg)
        for rep in reps:
            if identity(joined, "".join(rep["seg"])) >= IDENTITY:
                rep["weight"] += count
                break
        else:
            reps.append({"seg": seg, "weight": count})
    return reps


def main():
    with open(os.path.join(HERE, "tcr_seed_ids.txt")) as f:
        ids = [l.strip() for l in f if l.strip()]
    found = {f[0]: collections.Counter() for f in FAMILIES}
    dropped = collections.Counter()
    lengths = {"alpha": collections.defaultdict(list), "beta": collections.defaultdict(list)}
    for pid in ids:
        try:
            path = fetch(pid)
        except Exception:
            continue
        for _, k, res in chains(path):
            if k in lengths and anchored(res):
                fam = family_of(res, k)
                if fam is None:
                    dropped[k] += 1
                    continue
                found[fam][tuple(framework(res, fam))] += 1
                for name, n in loop_lengths(res).items():
                    lengths[k][name].append(n)
    print("dropped (fit no family):", dict(dropped))
    for k in lengths:
        for name, v in lengths[k].items():
            v.sort()
            q = lambda p: v[min(len(v) - 1, int(p * len(v)))]
            print(k, name, "n=%d min=%d p10=%d p90=%d max=%d" % (len(v), v[0], q(0.1), q(0.9), v[-1]))
    write({k: cluster(c) for k, c in found.items()})


def write(clusters):
    lines = [
        "// Generated by gen_tcr_seeds.py from wwPDB coordinates (CC0); do not edit.",
        "",
        "use super::profile::Seed;",
        "",
    ]
    for name, *_ in FAMILIES:
        lines.append("#[rustfmt::skip]")
        lines.append("pub(super) static %s: &[Seed] = &[" % name)
        for rep in clusters[name]:
            fr = ", ".join('"%s"' % x for x in rep["seg"])
            lines.append("    Seed { weight: %d, fr: [%s] }," % (rep["weight"], fr))
        lines += ["];", ""]
    with open(os.path.join(HERE, "tcr_seeds.rs"), "w", newline="\n") as f:
        f.write("\n".join(lines))
    for k in clusters:
        print(k, len(clusters[k]), "clusters,", sum(r["weight"] for r in clusters[k]), "chains")


if __name__ == "__main__":
    main()
