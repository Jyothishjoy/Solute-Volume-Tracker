# Solute Volume Tracker

A Python script for computing the **Solvent-Excluded Surface (SES) volume and surface area** of a solute molecule as a function of simulation time, across multiple XYZ trajectory files simultaneously. Designed for quasiclassical trajectory analysis where small differences in molecular volume expansion between reaction pathways are physically meaningful.

---

## Background and Motivation

In direct dynamics trajectory simulations, explicit solvent molecules can impose **dynamical steric barriers** that alter product selectivity. To quantify this effect, one needs a high-precision measure of how much space the solute occupies — and how that changes over time — for each product-forming pathway.

This script computes the **Solvent-Excluded Surface (SES) volume** (also called the Connolly surface volume) using the MSMS program (Sanner, 1996). The SES is defined by rolling a probe sphere of radius 1.4 Å (representing the solvent) over the van der Waals surface of the molecule. The enclosed volume represents the space the solvent physically cannot penetrate — making it the correct descriptor for solvent-barrier hypotheses.

MSMS computes this surface **analytically**, not on a grid, giving reproducibility of ~10⁻⁴ Å³ per frame. This precision is necessary when comparing small volume differences between reaction pathways.

### Atomic radii

Van der Waals radii are taken from two sources, deliberately mixed for consistency:

- **Bondi (1964)** — *J. Phys. Chem.*, 68, 441 — for H, C, N, O
- **Alvarez (2013)** — *Dalton Trans.*, 42, 8617 — for Fe and other transition metals

---

## What the script does

1. Discovers all `.xyz` trajectory files in a specified directory automatically.
2. For each trajectory, iterates through every frame using **MDAnalysis**.
3. Extracts the Cartesian coordinates of the first `N_SOLUTE` atoms (the solute) from each frame.
4. Calls **MSMS** analytically to compute the SES volume (Å³) and SES area (Å²) for each frame.
5. Assigns a simulation time to each frame based on a fixed timestep.
6. Writes two output CSV files:
   - `all_volumes.csv` — SES volume vs. time for all trajectories
   - `all_areas.csv`  — SES surface area vs. time for all trajectories

Each CSV has `time_fs` as the first column, followed by one column per trajectory file. This layout is directly ready for plotting in Excel, Origin, or Python.

### Output format

```
time_fs, pathway_A.xyz, pathway_B.xyz, pathway_C.xyz
0.00,    514.1632,      508.2241,      521.4403
0.75,    514.8821,      507.9934,      522.1107
1.50,    515.3301,      508.4412,      523.8821
...
```

---

## Requirements

| Software | Purpose |
|---|---|
| Python ≥ 3.10 | Runtime |
| MSMS 2.6.1 | Analytical SES computation |
| MDAnalysis | XYZ trajectory I/O and frame iteration |
| NumPy | Array operations |

---

## Installation

### 1. Install Miniconda (if not already installed)

Download from [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html) and follow the installer instructions for your operating system.

### 2. Create the conda environment

```bash
conda create -n vol_analysis python=3.11 -y
conda activate vol_analysis
```

### 3. Install MSMS and MDAnalysis via conda

Install these together so conda can resolve dependencies correctly:

```bash
conda install -c bioconda -c conda-forge msms mdanalysis -y
```

### 4. Install remaining packages via pip

```bash
pip install numpy
```

### 5. Verify MSMS is accessible

```bash
msms -h
```

You should see the MSMS usage message. If you get `command not found`, the bioconda install did not complete correctly — retry Step 3.

---

## Usage

### 1. Prepare your trajectory files

Place all `.xyz` trajectory files in a single directory. Each file must be a standard multi-frame XYZ file where:

- Line 1 of each frame: total atom count
- Line 2 of each frame: comment line (can be empty)
- Lines 3 onward: element symbol followed by X Y Z coordinates

The solute atoms must be the **first N atoms** in every frame, in a consistent order across all files.

### 2. Configure the script

Open `traj_vol_extract.py` and edit the Settings block at the top:

```python
TRAJ_DIR     = '.'      # path to the directory containing your .xyz files
                        # '.' means the current working directory
N_SOLUTE     = 74       # number of solute atoms (first N atoms in each frame)
PROBE_RADIUS = 1.4      # solvent probe radius in Å (1.4 Å = water/ACN)
MSMS_DENSITY = 3.0      # triangulation density; increase to 5.0 for higher precision
TIMESTEP_FS  = 0.75     # time between frames in femtoseconds
OUT_VOL_CSV  = 'all_volumes.csv'
OUT_AREA_CSV = 'all_areas.csv'
```

If your system contains elements other than H, C, N, O, Fe, add their van der Waals radii (from Alvarez 2013) to the `BONDI` dictionary:

```python
BONDI = {
    'H' : 1.20,
    'C' : 1.70,
    'N' : 1.55,
    'O' : 1.52,
    'Fe': 2.05,
    'Cu': 1.96,   # example: add any additional elements here
}
```

### 3. Run the script

```bash
conda activate vol_analysis
cd /path/to/your/trajectory/directory
python traj_vol_extract.py
```

Progress is printed to the terminal frame by frame. On completion, two CSV files are written to the current directory.

---

## Adapting to other systems

**Different solute size:** change `N_SOLUTE` to match the number of atoms in your solute. The script always reads the first `N_SOLUTE` atoms from each frame, so trajectories containing explicit solvent atoms are handled automatically — the solvent atoms are simply ignored.

**Solute-only trajectories:** these work identically. The script does not require solvent atoms to be present.

**Different probe radius:** 1.4 Å is appropriate for water and acetonitrile. For larger solvents, use their corresponding van der Waals radius.

**Different timestep:** change `TIMESTEP_FS` to match your simulation timestep. The script computes time as `frame_index × TIMESTEP_FS`.

**Atom label format:** MDAnalysis reads element symbols from column 1 of the XYZ file. If your XYZ files use labels like `C1`, `H12` instead of bare `C`, `H`, add these lines before the trajectory loop:

```python
from MDAnalysis.topology.guessers import guess_types
elems = guess_types(u.atoms.names)
u.add_TopologyAttr('elements', elems)
# then replace solute.names with solute.elements in the loop
```

---

## References

- Sanner, M. F., Olson, A. J., Spehner, J. C. (1996). Reduced surface: an efficient way to compute molecular surfaces. *Biopolymers*, 38, 305–320.
- Bondi, A. (1964). Van der Waals volumes and radii. *J. Phys. Chem.*, 68, 441–451.
- Alvarez, S. (2013). A cartography of the van der Waals territories. *Dalton Trans.*, 42, 8617–8636.
- Michaud-Agrawal, N. et al. (2011). MDAnalysis: A toolkit for the analysis of molecular dynamics simulations. *J. Comput. Chem.*, 32, 2319–2327.
