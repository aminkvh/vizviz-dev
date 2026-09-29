"""Writes membrane.pdb, a hand-built classification fixture: residue names
and atom names as the force fields spell them, geometry only as good as
bond perception needs (1.5 A zig-zag chains, planar peptide links).

    py -3 fixtures/membrane/make.py

Chains: P peptide with an unknown-named modified residue, L CHARMM POPC
and CHL1, A Amber Lipid21-style PC + PA + OL split lipid, U an unknown
lipid-like residue, W TIP3 waters and an unknown-named water, I SOD/CLA,
X unknown ligands (a phosphate ion and a benzene-like ring).
"""

from math import cos, radians, sin
from pathlib import Path

lines = []
serial = 0


def atom(chain, resseq, resname, name, x, y, z, element, hetero=True):
    global serial
    serial += 1
    record = "HETATM" if hetero else "ATOM  "
    shown = f" {name:<3}" if len(name) < 4 else name
    lines.append(
        f"{record}{serial:5d} {shown} {resname:<4}{chain}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}"
    )


def residue(chain, resseq, resname, atoms, origin, hetero=True, zsign=1):
    ox, oy, oz = origin
    for name, element, (x, y, z) in atoms:
        atom(chain, resseq, resname, name, ox + x, oy + y, oz + zsign * z, element, hetero)


def chain_atoms(prefix, first, count, element="C"):
    """`count` carbons named prefix+first.. along +z, 1.53 A apart."""
    return [
        (f"{prefix}{first + i}", element, (0.88 if i % 2 else 0.0, 0.0, 1.25 * (i + 1)))
        for i in range(count)
    ]


def water(chain, resseq, resname, names, origin):
    a = radians(104.5) / 2
    residue(
        chain,
        resseq,
        resname,
        [
            (names[0], "O", (0.0, 0.0, 0.0)),
            (names[1], "H", (0.9572 * sin(a), 0.9572 * cos(a), 0.0)),
            (names[2], "H", (-0.9572 * sin(a), 0.9572 * cos(a), 0.0)),
        ],
        origin,
    )


def popc(resseq, origin, chain="L", zsign=1):
    head = [
        ("N", "N", (0.0, 0.0, 0.0)),
        ("C13", "C", (1.3, 0.0, 0.5)),
        ("C14", "C", (-0.7, 1.2, 0.5)),
        ("C15", "C", (-0.7, -1.2, 0.5)),
        ("C12", "C", (0.1, 0.0, -1.5)),
        ("C11", "C", (0.1, 0.0, -3.0)),
        ("P", "P", (0.1, 0.0, -4.6)),
        ("O13", "O", (1.5, 0.0, -4.6)),
        ("O14", "O", (-0.7, 1.2, -4.6)),
        ("O11", "O", (0.1, 0.0, -6.0)),
        ("O12", "O", (0.1, -1.3, -4.0)),
        ("C1", "C", (0.1, 0.0, -7.4)),
        ("C2", "C", (0.1, 0.0, -8.9)),
        ("C3", "C", (0.1, 0.0, -10.4)),
        ("O21", "O", (1.4, 0.0, -9.4)),
        ("C21", "C", (2.5, 0.0, -9.0)),
        ("O22", "O", (2.5, 0.0, -7.8)),
        ("O31", "O", (0.1, 0.0, -11.7)),
        ("C31", "C", (0.1, 0.0, -13.0)),
        ("O32", "O", (1.3, 0.0, -13.0)),
    ]
    sn2 = [(f"C2{i}", "C", (3.0 + (0.88 if i % 2 else 0.0), 0.0, -9.0 - 1.25 * (i - 1))) for i in range(2, 19)]
    sn1 = [(f"C3{i}", "C", (0.1 + (0.88 if i % 2 else 0.0), 0.0, -13.0 - 1.25 * (i - 1))) for i in range(2, 17)]
    residue(chain, resseq, "POPC", head + sn2 + sn1, origin, zsign=zsign)


def chl1(resseq, origin, zsign=1):
    core = [(f"C{i}", "C", ((0.88 if i % 2 else 0.0), (i // 5) * 1.4, 1.25 * (i % 5))) for i in range(1, 28)]
    core.append(("O3", "O", (-1.3, 0.0, -0.3)))
    residue("L", resseq, "CHL1", core, origin, zsign=zsign)


def peptide():
    def res(k, name, n, cb=True):
        residue(
            "P",
            k,
            name,
            [
                ("N", "N", (0.0, 0.0, 0.0)),
                ("CA", "C", (1.46, 0.0, 0.0)),
                ("C", "C", (2.22, 1.32, 0.0)),
                ("O", "O", (2.22, 2.55, 0.0)),
            ]
            + ([("CB", "C", (1.46, -1.53, 0.0))] if cb else []),
            (3.42 * n, 0.74 * n, 0.0),
            hetero=False,
        )

    res(1, "ALA", 0)
    res(2, "XYZ", 1)
    res(3, "GLY", 2, cb=False)


def lipid21(origin):
    ox, oy, oz = origin
    pc = [
        ("N31", "N", (0.0, 0.0, 0.0)),
        ("C32", "C", (1.3, 0.0, 0.5)),
        ("C33", "C", (-0.7, 1.2, 0.5)),
        ("C34", "C", (-0.7, -1.2, 0.5)),
        ("C35", "C", (0.1, 0.0, -1.5)),
        ("C36", "C", (0.1, 0.0, -3.0)),
        ("P31", "P", (0.1, 0.0, -4.6)),
        ("O33", "O", (1.5, 0.0, -4.6)),
        ("O34", "O", (-0.7, 1.2, -4.6)),
        ("O32", "O", (0.1, 0.0, -6.0)),
        ("C31", "C", (0.1, 0.0, -7.4)),
        ("C2", "C", (0.1, 0.0, -8.9)),
        ("C1", "C", (0.1, 0.0, -10.4)),
    ]
    residue("A", 1, "PC", pc, origin)
    tail = lambda first: [
        ("C21", "C", (0.0, 0.0, 0.0)),
        ("O21", "O", (1.2, 0.0, 0.0)),
        ("O22", "O", (-0.5, -1.1, 0.0)),
    ] + chain_atoms("C", first, 15)
    residue("A", 2, "PA", [(n, e, (x, y, -z)) for n, e, (x, y, z) in tail(22)], (ox + 3.0, oy, oz - 8.5))
    residue("A", 3, "OL", tail(32), (ox - 3.0, oy, oz - 10.4))


def unknown_lipid(origin):
    head = [
        ("P", "P", (0.0, 0.0, 0.0)),
        ("O1", "O", (1.4, 0.0, 0.0)),
        ("O2", "O", (-0.7, 1.2, 0.0)),
        ("O3", "O", (-0.7, -1.2, 0.0)),
        ("O4", "O", (0.0, 0.0, 1.5)),
        ("C1", "C", (0.0, 0.0, 2.9)),
        ("C2", "C", (0.0, 0.0, 4.4)),
        ("O5", "O", (1.3, 0.0, 4.8)),
        ("C3", "C", (2.5, 0.0, 4.3)),
        ("O6", "O", (2.5, 0.0, 3.1)),
    ]
    tail = [(f"CA{i}", "C", (3.0 + (0.88 if i % 2 == 0 else 0.0), 0.0, 4.3 + 1.25 * i)) for i in range(1, 15)]
    residue("U", 1, "LPX", head + tail, origin)


def ligands(origin):
    ox, oy, oz = origin
    po4 = [
        ("P", "P", (0.0, 0.0, 0.0)),
        ("O1", "O", (1.5, 0.0, 0.0)),
        ("O2", "O", (-0.75, 1.3, 0.0)),
        ("O3", "O", (-0.75, -1.3, 0.0)),
        ("O4", "O", (0.0, 0.0, 1.5)),
    ]
    residue("X", 1, "PO4", po4, origin)
    ring = [
        (f"C{i}", "C", (1.4 * cos(radians(60 * i)), 1.4 * sin(radians(60 * i)), 0.0))
        for i in range(6)
    ]
    residue("X", 2, "BZX", ring, (ox + 10.0, oy, oz))


peptide()
popc(1, (20.0, 0.0, 0.0))
popc(2, (26.0, 0.0, 0.0))
chl1(3, (32.0, 0.0, 0.0))
lipid21((50.0, 0.0, 0.0))
unknown_lipid((70.0, 0.0, 0.0))
for k in range(3):
    water("W", k + 1, "TIP3", ("OH2", "H1", "H2"), (0.0, 20.0 + 4.0 * k, 0.0))
water("W", 4, "SW9", ("OW", "HW1", "HW2"), (0.0, 34.0, 0.0))
atom("I", 1, "SOD", "SOD", 10.0, 40.0, 0.0, "NA")
atom("I", 2, "CLA", "CLA", 14.0, 40.0, 0.0, "CL")
ligands((0.0, 50.0, 0.0))
lines.append("END")

out = Path(__file__).with_name("membrane.pdb")
out.write_text("\n".join(lines) + "\n", newline="\n")
print(f"wrote {out} ({serial} atoms)")


def bilayer():
    """bilayer.pdb: two 6x6 POPC leaflets (tails toward the middle), a
    cholesterol in each, and a lone POPC (chain B) far from them, the way
    a lipid bound to a protein sits away from any membrane."""
    global serial
    lines.clear()
    serial = 0
    pitch, seq = 6.5, 0
    for zsign, z, shift in ((1, 20.0, 0.0), (-1, -20.0, 3.0)):
        for ix in range(6):
            for iy in range(6):
                seq += 1
                popc(seq, (pitch * ix + shift, pitch * iy + shift, z), zsign=zsign)
        seq += 1
        chl1(seq, (36.5 + shift, shift, z), zsign=zsign)
    popc(1, (300.0, 0.0, 0.0), chain="B")
    lines.append("END")
    out = Path(__file__).with_name("bilayer.pdb")
    out.write_text("\n".join(lines) + "\n", newline="\n")
    print(f"wrote {out} ({serial} atoms)")


bilayer()
