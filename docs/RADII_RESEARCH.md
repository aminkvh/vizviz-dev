# Van der Waals display radii

`Element::vdw_radius` (`crates/vv-core/src/element.rs`) gives the sphere
size spacefill and the other atom styles scale (spacefill 1.0x, ball and
stick 0.25x). SASA is separate: it uses ProtOr radii (Tsai et al. 1999),
matching FreeSASA and FastSASA.

## Sources

| Source | Coverage | Method |
|---|---|---|
| Bondi, *J. Phys. Chem.* **68**, 441 (1964) (+ metals follow-up, **70**, 3006, 1966) | about 28 main-group elements and noble gases, plus Ni, Cu, Zn, Pd, Ag, Cd, Pt, Au, Hg, U | literature survey of contact distances |
| Rowland & Taylor, *J. Phys. Chem.* **100**, 7384 (1996) | revises 9 of Bondi's elements; H is the one large change (1.20 to 1.10), the other 8 shift by 0.09 A or less | more than 25,000 CSD contacts |
| Mantina et al., *J. Phys. Chem. A* **113**, 5806 (2009) | fills the remaining main-group gaps on Bondi's scale | scale-matched extension of Bondi |
| Batsanov, *Inorg. Mater.* **37**, 871 (2001) | about 65 elements including most transition metals, almost no lanthanides | crystallochemical contact-distance survey |
| Alvarez, *Dalton Trans.* **42**, 8617 (2013) | Z = 1-103, the most complete table | statistics over more than 5,000,000 CSD/ICSD contacts |

## The cascade

Per element, the first source that has it wins:

1. Bondi 1964, with Rowland and Taylor's hydrogen (1.10 A) and Bondi's own
   metals and Li (1.82).
2. Mantina 2009 for the main-group gaps.
3. Batsanov 2001 for the remaining d-block elements.
4. Alvarez 2013 for the f-block (La-Lu, Ac-Lr), except U, which stays at
   Bondi's 1.86 A (the only actinide common in real structures).
5. A flat 2.0 A placeholder where no source has data (Pm, Z >= 100).

Transition metals are where the sources disagree by more than rounding:
Cu is 1.40 A (Bondi), 2.00 A (Batsanov) or 2.38 A (Alvarez). Alvarez's Fe
(2.44 A) and Zn (2.39 A) exceed the Fe-N(heme) (about 2.0 A) and Zn-S(Cys)
(about 2.3 A) bond lengths, so in spacefill the metal's sphere would
swallow its ligand atom's centre. That is why Batsanov, not Alvarez, is used for the
d-block.

All radii above are neutral-atom radii, while metals in PDB entries are
almost always ions: a neutral Ca (2.31 A) nearly reaches a coordinating
water oxygen (Ca-O about 2.4 A). Ionic radii keyed by formal charge
(Shannon 1976) would be the general fix.

## Why force-field radii are not display radii

CHARMM and AMBER define one Lennard-Jones radius per atom type, not per
element. For hydrogen alone the CHARMM36 Rmin/2 values span about 6x
(0.2245 A for an amide H, 1.34 A for an aliphatic H), because a polar
hydrogen's bulk is folded into its parent atom and the electrostatics.
Some ion parameters give Mg2+ a radius under 1 A, smaller than hydrogen,
because they are fitted to hydration thermodynamics. Choosing the right type
needs the full topology, and the force fields disagree with each other, so
there is no single value to draw. They earn their keep in electrostatics:
`pdb2pqr` writes charges and radii from the same fitted set into a PQR
file for Poisson-Boltzmann solvers.

## Per-element table (Z 1-100)

Generated from the `mendeleev` package's Bondi, Mantina, Batsanov and
Alvarez columns, patched with Bondi's own metals and Li (`mendeleev`'s
Bondi field lacks them):

```python
from mendeleev import element
for z in range(1, 101):
    e = element(z)
    print(e.symbol, e.vdw_radius_bondi, e.vdw_radius_truhlar,
          e.vdw_radius_batsanov, e.vdw_radius_alvarez)
```

| Z | El | A | Source | Z | El | A | Source |
|---|----|-----|--------|---|----|-----|--------|
| 1 | H | 1.10 | Rowland-Taylor 1996 | 51 | Sb | 2.06 | Mantina 2009 |
| 2 | He | 1.40 | Bondi 1964 | 52 | Te | 2.06 | Bondi 1964 |
| 3 | Li | 1.82 | Bondi 1964 | 53 | I | 1.98 | Bondi 1964 |
| 4 | Be | 1.53 | Mantina 2009 | 54 | Xe | 2.16 | Bondi 1964 |
| 5 | B | 1.92 | Mantina 2009 | 55 | Cs | 3.43 | Mantina 2009 |
| 6 | C | 1.70 | Bondi 1964 | 56 | Ba | 2.68 | Mantina 2009 |
| 7 | N | 1.55 | Bondi 1964 | 57 | La | 2.98 | Alvarez 2013 |
| 8 | O | 1.52 | Bondi 1964 | 58 | Ce | 2.88 | Alvarez 2013 |
| 9 | F | 1.47 | Bondi 1964 | 59 | Pr | 2.92 | Alvarez 2013 |
| 10 | Ne | 1.54 | Bondi 1964 | 60 | Nd | 2.95 | Alvarez 2013 |
| 11 | Na | 2.27 | Bondi 1964 | 61 | Pm | 2.00 | fallback (no data anywhere) |
| 12 | Mg | 1.73 | Bondi 1964 | 62 | Sm | 2.90 | Alvarez 2013 |
| 13 | Al | 1.84 | Mantina 2009 | 63 | Eu | 2.87 | Alvarez 2013 |
| 14 | Si | 2.10 | Bondi 1964 | 64 | Gd | 2.83 | Alvarez 2013 |
| 15 | P | 1.80 | Bondi 1964 | 65 | Tb | 2.79 | Alvarez 2013 |
| 16 | S | 1.80 | Bondi 1964 | 66 | Dy | 2.87 | Alvarez 2013 |
| 17 | Cl | 1.75 | Bondi 1964 | 67 | Ho | 2.81 | Alvarez 2013 |
| 18 | Ar | 1.88 | Bondi 1964 | 68 | Er | 2.83 | Alvarez 2013 |
| 19 | K | 2.75 | Bondi 1964 | 69 | Tm | 2.79 | Alvarez 2013 |
| 20 | Ca | 2.31 | Mantina 2009 | 70 | Yb | 2.80 | Alvarez 2013 |
| 21 | Sc | 2.30 | Batsanov 2001 | 71 | Lu | 2.74 | Alvarez 2013 |
| 22 | Ti | 2.15 | Batsanov 2001 | 72 | Hf | 2.25 | Batsanov 2001 |
| 23 | V | 2.05 | Batsanov 2001 | 73 | Ta | 2.20 | Batsanov 2001 |
| 24 | Cr | 2.05 | Batsanov 2001 | 74 | W | 2.10 | Batsanov 2001 |
| 25 | Mn | 2.05 | Batsanov 2001 | 75 | Re | 2.05 | Batsanov 2001 |
| 26 | Fe | 2.05 | Batsanov 2001 | 76 | Os | 2.00 | Batsanov 2001 |
| 27 | Co | 2.00 | Batsanov 2001 | 77 | Ir | 2.00 | Batsanov 2001 |
| 28 | Ni | 1.63 | Bondi 1964 | 78 | Pt | 1.72 | Bondi 1964 |
| 29 | Cu | 1.40 | Bondi 1964 | 79 | Au | 1.66 | Bondi 1964 |
| 30 | Zn | 1.39 | Bondi 1964 | 80 | Hg | 1.55 | Bondi 1964 |
| 31 | Ga | 1.87 | Bondi 1964 | 81 | Tl | 1.96 | Bondi 1964 |
| 32 | Ge | 2.11 | Mantina 2009 | 82 | Pb | 2.02 | Bondi 1964 |
| 33 | As | 1.85 | Bondi 1964 | 83 | Bi | 2.07 | Mantina 2009 |
| 34 | Se | 1.90 | Bondi 1964 | 84 | Po | 1.97 | Mantina 2009 |
| 35 | Br | 1.85 | Bondi 1964 | 85 | At | 2.02 | Mantina 2009 |
| 36 | Kr | 2.02 | Bondi 1964 | 86 | Rn | 2.20 | Mantina 2009 |
| 37 | Rb | 3.03 | Mantina 2009 | 87 | Fr | 3.48 | Mantina 2009 |
| 38 | Sr | 2.49 | Mantina 2009 | 88 | Ra | 2.83 | Mantina 2009 |
| 39 | Y | 2.40 | Batsanov 2001 | 89 | Ac | 2.80 | Alvarez 2013 |
| 40 | Zr | 2.30 | Batsanov 2001 | 90 | Th | 2.93 | Alvarez 2013 |
| 41 | Nb | 2.15 | Batsanov 2001 | 91 | Pa | 2.88 | Alvarez 2013 |
| 42 | Mo | 2.10 | Batsanov 2001 | 92 | U | 1.86 | Bondi 1964 |
| 43 | Tc | 2.05 | Batsanov 2001 | 93 | Np | 2.82 | Alvarez 2013 |
| 44 | Ru | 2.05 | Batsanov 2001 | 94 | Pu | 2.81 | Alvarez 2013 |
| 45 | Rh | 2.00 | Batsanov 2001 | 95 | Am | 2.83 | Alvarez 2013 |
| 46 | Pd | 1.63 | Bondi 1964 | 96 | Cm | 3.05 | Alvarez 2013 |
| 47 | Ag | 1.72 | Bondi 1964 | 97 | Bk | 3.40 | Alvarez 2013 |
| 48 | Cd | 1.58 | Bondi 1964 | 98 | Cf | 3.05 | Alvarez 2013 |
| 49 | In | 1.93 | Bondi 1964 | 99 | Es | 2.70 | Alvarez 2013 |
| 50 | Sn | 2.17 | Bondi 1964 | 100 | Fm | 2.00 | fallback (no data anywhere) |

Pm (Z = 61) and Z >= 100 have no contact data to fit a radius from; they
are rare or synthetic and do not occur in real structures, so the flat
2.0 A placeholder stands.

## References

Bondi 1964 and 1966, Rowland and Taylor 1996, Mantina et al. 2009,
Batsanov 2001, Alvarez 2013 (full citations above); Tsai et al., *J. Mol.
Biol.* **290**, 253 (1999) for ProtOr; Shannon, *Acta Cryst. A* **32**, 751
(1976) for ionic radii. Cross-checked with the `mendeleev` package and the
`ase` package's Alvarez table.
