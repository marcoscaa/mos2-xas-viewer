#!/usr/bin/env python3
"""Build data.json for the XAS viewer from calculations/systems/*.

For each system with at least one edge that has spectrum.dat:
  - parse the relaxed structure (01_relax/pw.relax.out [+ .restart.out]) into a CIF string
  - collect every edge's spectrum.dat (energy, avg, soc, soc+gauss)
and write it all into data.json (structures + spectra, no external files needed).
"""
import json
import re
import tempfile
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import write as ase_write

CALC_ROOT = Path(__file__).resolve().parent.parent / "calculations"
SYSTEMS_ROOT = CALC_ROOT / "systems"

READY_SYSTEMS = [
    "pristine",
    "pristine_2H",
    "adsorbates/Co_hollow",
    "adsorbates/Co_Motop",
    "adsorbates/Co_Stop",
    "adsorbates/Cr_hollow",
    "adsorbates/Cr_Motop",
    "adsorbates/Cu_hollow",
    "adsorbates/Cu_Motop",
    "adsorbates/Cu_Stop",
    "adsorbates/Fe_hollow",
    "adsorbates/Fe_Motop",
    "adsorbates/Mn_hollow",
    "adsorbates/Mn_Motop",
    "adsorbates/Mn_Stop",
    "adsorbates/Ti_hollow",
    "adsorbates/Ti_Stop",
    "adsorbates/V_hollow",
    "adsorbates/V_Stop",
    "adsorbates/Zn_hollow",
    "adsorbates/Zn_Motop",
    "adsorbates/Zn_Stop",
    "substitutional/Co_at_Mo",
    "substitutional/Cr_at_Mo",
    "substitutional/Cr_at_S",
    "substitutional/Cu_at_Mo",
    "substitutional/Cu_at_S",
    "substitutional/Fe_at_S",
    "substitutional/Mn_at_S",
    "substitutional/Ti_at_Mo",
    "substitutional/Ti_at_S",
    "substitutional/V_at_Mo",
    "substitutional/V_at_S",
    "substitutional/Zn_at_Mo",
    "substitutional/Zn_at_S",
]

FINAL_BLOCK_RE = re.compile(
    r"Begin final coordinates(.*?)End final coordinates", re.DOTALL
)


def parse_cell_from_input(relax_in: Path):
    text = relax_in.read_text()
    m = re.search(r"CELL_PARAMETERS\s*\S*\s*\n(.*?)\n(?:\n|[A-Z])", text, re.DOTALL)
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.strip().startswith("CELL_PARAMETERS"):
            vecs = [list(map(float, lines[i + 1 + k].split())) for k in range(3)]
            return np.array(vecs)
    raise ValueError(f"no CELL_PARAMETERS found in {relax_in}")


def last_final_block(out_texts):
    blocks = []
    for text in out_texts:
        blocks.extend(FINAL_BLOCK_RE.findall(text))
    if not blocks:
        raise ValueError("no 'Begin final coordinates' block found")
    return blocks[-1]


def parse_final_block(block, fallback_cell):
    lines = [l for l in block.splitlines() if l.strip()]
    lines = [l for l in lines if not l.strip().startswith(("new unit-cell", "density"))]
    cell = fallback_cell
    i = 0
    if lines[0].strip().startswith("CELL_PARAMETERS"):
        cell = np.array([list(map(float, lines[i + 1 + k].split())) for k in range(3)])
        i += 4
    assert lines[i].strip().startswith("ATOMIC_POSITIONS"), lines[i]
    units = "crystal" if "crystal" in lines[i] else "angstrom"
    i += 1
    symbols, coords = [], []
    for line in lines[i:]:
        parts = line.split()
        if len(parts) != 4:
            continue
        symbols.append(parts[0])
        coords.append(list(map(float, parts[1:4])))
    coords = np.array(coords)
    if units == "crystal":
        atoms = Atoms(symbols=symbols, scaled_positions=coords, cell=cell, pbc=True)
    else:
        atoms = Atoms(symbols=symbols, positions=coords, cell=cell, pbc=True)
    return atoms


def get_relaxed_structure(sys_dir: Path) -> Atoms:
    relax_dir = sys_dir / "01_relax"
    relax_in = relax_dir / "pw.relax.in"
    fallback_cell = parse_cell_from_input(relax_in)

    out_texts = []
    main_out = relax_dir / "pw.relax.out"
    if main_out.exists():
        out_texts.append(main_out.read_text())
    restart_out = relax_dir / "pw.relax.restart.out"
    if restart_out.exists():
        out_texts.append(restart_out.read_text())

    block = last_final_block(out_texts)
    return parse_final_block(block, fallback_cell)


def atoms_to_cif_string(atoms: Atoms) -> str:
    with tempfile.NamedTemporaryFile(suffix=".cif") as tmp:
        ase_write(tmp.name, atoms, format="cif")
        return Path(tmp.name).read_text()


def is_stale(edge_dir: Path) -> bool:
    """True if any xspectra polarization predates the edge's own
    pw.scf.coreh.out -- i.e. it was computed against a since-superseded
    core-hole SCF (e.g. a relaxed-coordinate fix triggered a coreh rerun,
    but not all of x/y/z have been resubmitted yet). Averaging polarizations
    from two different coreh states is not physically meaningful, so such an
    edge must be excluded until all three are refreshed."""
    coreh_out = edge_dir / "pw.scf.coreh.out"
    if not coreh_out.exists():
        return False
    coreh_mtime = coreh_out.stat().st_mtime
    for c in "xyz":
        out = edge_dir / f"xanes_{c}.out"
        if not out.exists() or out.stat().st_mtime < coreh_mtime:
            return True
    return False


def load_edges(sys_dir: Path):
    edges = {}
    xanes_dir = sys_dir / "03_xanes"
    for d in sorted(xanes_dir.glob("*/")):
        spec = d / "spectrum.dat"
        if not spec.exists():
            continue
        if is_stale(d):
            print(f"  [skip] {d.name}: stale vs. current pw.scf.coreh.out")
            continue
        data = np.loadtxt(spec)
        edges[d.name] = {
            "energy": data[:, 0].round(4).tolist(),
            "avg": data[:, 1].tolist(),
            "soc": data[:, 2].tolist(),
            "soc_gauss": data[:, 3].tolist(),
        }
    return edges


def main():
    out = {"systems": {}}
    for key in READY_SYSTEMS:
        sys_dir = SYSTEMS_ROOT / key
        print(f"processing {key} ...")
        edges = load_edges(sys_dir)
        if not edges:
            print(f"  skipping {key}: no self-consistent edges right now")
            continue
        try:
            atoms = get_relaxed_structure(sys_dir)
        except ValueError as e:
            print(f"  skipping {key}: edges are self-consistent but ion relaxation "
                  f"hasn't converged yet ({e})")
            continue
        cif = atoms_to_cif_string(atoms)
        out["systems"][key] = {
            "formula": atoms.get_chemical_formula(),
            "natoms": len(atoms),
            "cif": cif,
            "edges": edges,
        }
        print(f"  {atoms.get_chemical_formula()}  edges={list(edges)}")

    out_path = Path(__file__).resolve().parent / "data.json"
    out_path.write_text(json.dumps(out))
    print(f"\nwrote {out_path} ({out_path.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
