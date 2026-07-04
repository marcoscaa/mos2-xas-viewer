# MoS2 XAS Viewer

Static, self-contained viewer for relaxed structures + simulated XAS spectra
of pristine and doped MoS2 systems. Pick a system and an absorption edge to see
the 3D structure (3Dmol.js) and spectrum (Plotly.js) side by side.

Live at: https://marcoscaa.github.io/mos2-xas-viewer/

## Regenerating data.json

`data.json` is built from the QE calculations in `../calculations/systems/`.
Whenever more systems finish post-processing (`scripts/run_xas.sh post <sys>`
in `calculations/`), add them to `READY_SYSTEMS` in `build_data.py` and rerun:

```
python3 build_data.py
```

Requires `ase` and `numpy`. Then commit and push `data.json`.
