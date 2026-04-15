# Atomistic Shearing of Al Sheet with W Punch and W L-Die

[![LAMMPS](https://img.shields.io/badge/LAMMPS-Molecular%20Dynamics-blue)](https://www.lammps.org/)
[![Language](https://img.shields.io/badge/Language-LAMMPS%20Script-orange)](https://docs.lammps.org/)
[![Material](https://img.shields.io/badge/Material-Al%2FW%20System-darkred)](https://en.wikipedia.org/wiki/Sheet_metal)
[![Potential](https://img.shields.io/badge/Potential-EAM%2FAlloy-green)](https://docs.lammps.org/pair_eam.html)

A LAMMPS input script for performing **uniaxial shearing** of an aluminum sheet using a tungsten punch and a tungsten L-die at the atomistic scale. The script models a periodic Al FCC sheet clamped between a rigid W die (left/bottom) and a downward-driven W punch (right), capturing full force–displacement response, per-atom von Mises stress evolution, and punch force partitioning against the Al mobile region.
<img width="1600" height="1200" alt="shearing2" src="https://github.com/user-attachments/assets/4ce63ae3-7046-4d7e-b94a-d8baa7c928ed" />

---

## Table of Contents

- [Overview](#overview)
- [Key Fixes and Design Decisions](#key-fixes-and-design-decisions)
- [Simulation Stages](#simulation-stages)
- [System Architecture](#system-architecture)
- [Parameters](#parameters)
- [Requirements](#requirements)
- [Usage](#usage)
- [Outputs](#outputs)
- [Visualization](#visualization)
- [Post-processing](#post-processing-python-snippet)

---

## Overview

This simulation builds an FCC Aluminum sheet (40 Å thick, 300 Å long) with two embedded rigid Tungsten bodies: a fixed L-shaped die on the left and a displacement-driven punch descending from the right. Shearing is applied by moving the W punch linearly downward at 0.5 Å/ps until its bottom face clears 5 Å below the top of the die. Per-atom von Mises stress and punch-on-Al group forces are written at regular intervals.

The script is designed to **avoid common MD shearing artifacts** — including grip drift under NVT, null-force reads before the punch engages, premature halt condition triggers, and rigid body integration conflicts with thermostat groups.

---

## Key Fixes and Design Decisions

| # | Issue | Fix Applied |
|---|-------|-------------|
| 1 | `fix rigid/nve` conflicting with `fix nvt` over shared atoms | `fix nvt` applied exclusively to `al_mobile`; die and punch groups are fully disjoint |
| 2 | Al support atoms drifting downward under punch load | `fix setforce 0 0 0` applied to `al_fixed` (bottom 5 Å of sheet over die) every step |
| 3 | Punch drifting laterally or in Z during descent | `fix drive_punch` uses `move linear 0 -0.5 0` — only Y displacement; `freeze_die` pins the die |
| 4 | `halt` condition triggering before punch engages sheet | `v_punch_bottom` evaluated every 100 steps; condition is `< punch_stop` (−6.0 Å), well below sheet top |
| 5 | Overlap between freshly created BCC and FCC lattices | `delete_atoms overlap 2.3 all all` removes any atoms closer than 2.3 Å across phase boundaries |

---

## Simulation Stages

```
┌─────────────────────────────────────────────────────────────────────┐
│  Stage 1  │  Energy Minimization (CG, implicit via initial overlap   │
│           │  removal + velocity assignment)                          │
├─────────────────────────────────────────────────────────────────────┤
│  Stage 2  │  NVT Equilibration — 10,000 steps @ 300 K (5 ps)       │
│           │  Punch held fixed; Al mobile thermalised                 │
├─────────────────────────────────────────────────────────────────────┤
│  Stage 3  │  Shearing Run — 100,000 steps (100 ps)                  │
│           │  Punch descends at 0.5 Å/ps; halts at punch_bottom < −6 │
└─────────────────────────────────────────────────────────────────────┘
```

---

## System Architecture

The simulation box is divided into functional zones along Y:

```
  Y ▲
    │   ┌───────────────────────────────────────────┐
    │   │         W Punch  (w_punch)                │  descends at −0.5 Å/ps
    │   │         x: [die_right+6 → lx−3]           │  rigid/nve + drive_punch
    │   │         y: [sheet_y1+1 → sheet_y1+100]    │
    │   ├───────────────────────────────────────────┤  y = +41 Å  (sheet top + 1)
    │   │   Al Sheet — mobile region (al_mobile)    │
    │   │   NVT thermostat @ 300 K                  │  y: [0 → 40 Å]
    │   │   Al Sheet — fixed region (al_fixed)      │
    │   │   x: [0 → die_right], y: [0 → 5 Å]       │
    │   ├───────────────────────────────────────────┤  y = 0 Å
    │   │                                           │
    │   │         W L-Die  (w_die)                  │  fully frozen
    │   │         x: [0 → 150 Å], y: [−120 → −1]   │  rigid/nve + setforce 0 0 0
    │   └───────────────────────────────────────────┘  y = −122 Å (box floor)
```

**Clearance gap:** A 6 Å gap between the die right edge (x = 150 Å) and punch left edge (x = 156 Å) defines the shear zone where the Al sheet is cut. The sheet material in this gap is the primary deformation region.

---

## Parameters

### Physical

| Parameter | Value | Description |
|-----------|-------|-------------|
| `lx` | 300.0 Å | Lateral (X) box dimension |
| `lz` | 60.0 Å | Depth (Z, periodic) dimension |
| `t` | 40.0 Å | Al sheet thickness (Y) |
| `die_right` | 150.0 Å | X position of die right edge |
| `c` | 6.0 Å | Clearance gap between die and punch |
| `T` | 300 K | Simulation temperature |
| `dt` (equil) | 0.0005 ps | Timestep during equilibration |
| `dt` (shear) | 0.001 ps | Timestep during shearing |
| Punch velocity | −0.5 Å/ps | Downward (−Y) punch speed |
| `punch_stop` | −6.0 Å | Y halt threshold for punch bottom face |

### Geometry Summary

| Component | Lattice | Constant | Region |
|-----------|---------|----------|--------|
| Al sheet | FCC | 4.05 Å | x:[0,300], y:[0,40], z:[0,60] |
| W die | BCC | 3.165 Å | x:[0,150], y:[−120,−1], z:[0,60] |
| W punch | BCC | 3.165 Å | x:[156,297], y:[41,141], z:[0,60] |

### Simulation Length

| Stage | Steps | Timestep | Wall Time (approx.) |
|-------|-------|----------|---------------------|
| Equilibration | 10,000 | 0.5 fs | 5 ps |
| Shearing | up to 100,000 | 1.0 fs | up to 100 ps |

---

## Requirements

- **LAMMPS** (any recent stable release, compiled with `MANYBODY` and `RIGID` packages)
- **Potential file:** `CuAlW.txt` — EAM/alloy potential covering Cu, Al, and W
  - Confirm the element ordering in the file header matches the `pair_coeff * * CuAlW.txt Cu Al W` call
  - Type 1 = Cu (unused here), Type 2 = Al, Type 3 = W
- **OVITO** (recommended for visualization and CNA/von Mises analysis)

### LAMMPS Packages Required

```
MANYBODY    # for pair_style eam/alloy
RIGID       # for fix rigid/nve
```

---

## Usage

1. **Clone the repository**
   ```bash
   git clone https://github.com/<your-username>/<repo-name>.git
   cd <repo-name>
   ```

2. **Place the potential file** in the same directory as the input script:
   ```
   CuAlW.txt
   ```

3. **Run the simulation**
   ```bash
   # Serial
   lammps -in shear_AlW.lammps

   # Parallel (recommended for large box)
   mpirun -np 8 lammps -in shear_AlW.lammps
   ```

4. **Monitor thermo output** during shearing:
   ```
   Step   Temp   PE   KE   v_Fy_punch   v_punch_bottom
   ```
   - `v_Fy_punch` becomes increasingly negative as the punch loads the Al sheet
   - `v_punch_bottom` descends from ~41 Å toward the halt threshold of −6.0 Å
   - Watch for a force drop in `v_Fy_punch` signalling crack initiation / shear fracture

---

## Outputs

| File | Description |
|------|-------------|
| `dump.shear.*.cfg` | Per-atom dump every 500 steps (id, type, xs, ys, zs, von Mises stress) |
| `final_shear.data` | Final atomic configuration in LAMMPS data format |
| Thermo log | Step, Temp, PE, KE, punch Y-force, punch bottom Y — written to `log.lammps` |

### Per-atom quantities in shear dump

| Column | Quantity |
|--------|----------|
| `id` | Atom ID |
| `type` | Atom type (1 = Cu unused, 2 = Al, 3 = W) |
| `xs ys zs` | Scaled atomic coordinates |
| `v_vm` | Per-atom von Mises stress (bar·Å³) |

---

## Visualization

Open trajectory files in **OVITO**:

1. Load `dump.shear.*.cfg` (OVITO auto-detects the wildcard sequence)
2. Apply **Color Coding** by `v_vm` (von Mises stress) to visualize stress concentration in the shear zone and beneath the punch
3. Apply **Common Neighbor Analysis (CNA)** to track FCC → disordered transitions in the Al sheet as shearing progresses
4. Use the **Expression Select** modifier to isolate Al (Type == 2) or W (Type == 3) for phase-specific inspection
5. Plot the thermo log data in Python/MATLAB to extract:
   - Punch force–displacement curve (engineering load vs. stroke)
   - Onset of plastic deformation via von Mises stress threshold in Al
   - Fracture point from the force drop

---

## Post-processing (Python snippet)

```python
import numpy as np
import matplotlib.pyplot as plt

# Parse thermo log — adjust column indices to match your thermo_style
# Columns: Step Temp PE KE Fy_punch punch_bottom
log = []
with open("log.lammps") as f:
    capture = False
    for line in f:
        if line.startswith("Step"):
            capture = True
            continue
        if capture:
            try:
                log.append([float(x) for x in line.split()])
            except ValueError:
                capture = False

data = np.array(log)
step         = data[:, 0]
fy_punch     = data[:, 4]          # punch force on Al (eV/Å → convert to nN if needed)
punch_bottom = data[:, 5]          # punch bottom Y position (Å)

# Punch displacement = initial_punch_bottom - current_punch_bottom
punch_disp = punch_bottom[0] - punch_bottom  # positive = downward stroke

fig, ax = plt.subplots()
ax.plot(punch_disp, np.abs(fy_punch), color="steelblue")
ax.set_xlabel("Punch Displacement (Å)")
ax.set_ylabel("|F_y| (eV/Å)")
ax.set_title("Shear Force–Displacement Curve — Al/W System")
plt.tight_layout()
plt.savefig("shear_force_displacement.png", dpi=150)
```

---

## Citation

If you use this script or build upon it in your research, please cite:

```bibtex
@misc{mishra2026shearing,
  author    = {Mishra, A.},
  title     = {Atomistic Shearing of Al Sheet with W Punch and W L-Die},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.19591374},
  url       = {https://doi.org/10.5281/zenodo.19591374}
}
```

> Mishra, A. (2026). *Atomistic Shearing of Al Sheet with W Punch and W L-Die*. Zenodo. https://doi.org/10.5281/zenodo.19591374

---

## License

This project is open-source. Feel free to use and adapt the script for your own sheet-forming or fracture MD simulations. Attribution appreciated.
