"""A structure through FastSASA and back into vizviz, with no change to
FastSASA.

    pip install vizviz            # or the local wheel, see README.md
    pip install git+https://github.com/aminkvh/FastSASA
    python examples/fastsasa_pipe.py [PDB_ID]

The contract is NumPy in, NumPy out. vizviz hands FastSASA positions and
the names it needs for radii; FastSASA returns one SASA per atom; that
column comes back as a value channel the renderers color by. Prints a
few results, writes `sasa_<ID>.png` (CPU render, no GPU needed) and
`sasa_<ID>.npy`, which the desktop app takes with

    vizviz --exec "fetch 4HHB; values sasa sasa_4HHB.npy; color values sasa"
"""

import sys
import time

import fastsasa
import fastsasa_adapters as fa
import numpy as np
import vizviz

pdb_id = sys.argv[1] if len(sys.argv) > 1 else "4HHB"
s = vizviz.load(vizviz.fetch(pdb_id))
print(s)

# --- out: the columns FastSASA wants --------------------------------------
t0 = time.perf_counter()
positions = s.positions.astype(np.float64)  # (atoms, 3), Å
atom_names = np.char.strip(s.name.astype("U4"))
residue_names = np.asarray(s.residue_names)[s.residue_index]  # per residue -> per atom
symbols = np.array([vizviz.element_symbol(int(z)) or "X" for z in range(256)])[s.element]
radii = fa.load_radius_config().radii(residue_names, atom_names, elements=symbols, default=1.70)

# --- compute (FastSASA's own API, untouched) -------------------------------
total, atom = fastsasa.sasa(positions, radii, probe_radius=1.4, n_points=100, atom_sasa=True)
per_atom = atom[0]  # FastSASA returns (frames, atoms)
t1 = time.perf_counter()
print(f"total SASA {total[0]:.0f} Å², {len(per_atom)} atoms, {1e3 * (t1 - t0):.0f} ms end to end")

# Per residue, using vizviz's residue index as FastSASA's grouping.
grouped = fastsasa.sasa(
    positions,
    radii,
    atom_sasa=True,
    residue_ids=s.residue_index.astype(np.int32),
    n_residues=s.residue_count,
)
per_residue = grouped["residue"][0]
top = np.argsort(per_residue)[::-1][:5]
print("most exposed residues:")
for r in top:
    chain = s.chain_names[s.residue_chain[r]]
    print(
        f"  {s.residue_name(int(r))} {s.residue_seq_id[r]} chain {chain}: {per_residue[r]:.0f} Å²"
    )

# Interface buried area from two selection-language selections.
a, b = s.select("chain A and protein"), s.select("chain B and protein")
if a.size and b.size:

    def sasa_of(idx):
        return float(fastsasa.sasa(positions[idx], radii[idx])[0])

    buried = (sasa_of(a) + sasa_of(b) - sasa_of(np.concatenate([a, b]))) / 2
    print(f"A/B interface: {buried:.0f} Å² buried per side")

# --- back in: color by the result --------------------------------------------
if vizviz.HAS_RENDER:
    img = s.render(1200, 900, values=per_atom, style="publication_white")
    vizviz.write_png(f"sasa_{pdb_id}.png", img)
    print(f"wrote sasa_{pdb_id}.png (blue buried .. red exposed)")
np.save(f"sasa_{pdb_id}.npy", per_atom)
print(f"wrote sasa_{pdb_id}.npy; in the app: values sasa sasa_{pdb_id}.npy; color values sasa")

# The same channel in a Session (the desktop app's command bus).
session = vizviz.Session()
sid = session.load(vizviz.fetch(pdb_id))
session.set_values(sid, "sasa", per_atom)
print(session.exec("color values sasa; values"))
