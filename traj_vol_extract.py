import subprocess, tempfile, os, csv, glob
import numpy as np
import MDAnalysis as mda

# ── Bondi (1964) for organics + Alvarez (2013) for Fe ────────────────────────
BONDI = {
    'H' : 1.20,
    'C' : 1.70,
    'N' : 1.55,
    'O' : 1.52,
    'Fe': 2.05,
}

# ── Settings ──────────────────────────────────────────────────────────────────
TRAJ_DIR      = '.'               # directory containing all XYZ trajectories
                                  # '.' means the current working directory
                                  # change to e.g. r'C:\data\trajectories' if needed
N_SOLUTE      = 74
PROBE_RADIUS  = 1.4
MSMS_DENSITY  = 3.0
TIMESTEP_FS   = 0.75
OUT_VOL_CSV   = 'all_volumes.csv'
OUT_AREA_CSV  = 'all_areas.csv'

# ── MSMS runner ───────────────────────────────────────────────────────────────
def run_msms(coords, radii, probe_radius=1.4):
    """
    Write a temporary .xyzr file, call MSMS, parse stdout.
    Returns (ses_area_A2, ses_volume_A3).
    """
    with tempfile.TemporaryDirectory() as tmp:
        xyzr = os.path.join(tmp, 'mol.xyzr')
        out  = os.path.join(tmp, 'out')

        with open(xyzr, 'w') as f:
            for (x, y, z), r in zip(coords, radii):
                f.write(f'{x:.6f} {y:.6f} {z:.6f} {r:.4f}\n')

        cmd = ['msms',
               '-if', xyzr,
               '-probe_radius', str(probe_radius),
               '-density',      str(MSMS_DENSITY),
               '-of', out]

        result = subprocess.run(cmd, capture_output=True, text=True)

        area = vol = None

        for line in result.stdout.splitlines():
            stripped = line.strip()

            # Volume: "Total ses_volume:   514.163"
            if stripped.startswith('Total ses_volume:'):
                try:
                    vol = float(stripped.split()[-1])
                except ValueError:
                    pass

            # Area: numerical data row "   0      1.50      514.163    396.964"
            parts = stripped.split()
            if len(parts) == 4:
                try:
                    comp  = int(parts[0])
                    probe = float(parts[1])
                    s_area = float(parts[3])
                    if comp == 0 and abs(probe - probe_radius) < 0.01:
                        area = s_area
                except ValueError:
                    pass

        if vol is None:
            raise RuntimeError(
                f'MSMS failed at probe={probe_radius}\n'
                f'stdout:\n{result.stdout}\nstderr:\n{result.stderr}'
            )

        return area, vol


# ── Discover all XYZ files ────────────────────────────────────────────────────
xyz_files = sorted(glob.glob(os.path.join(TRAJ_DIR, '*.xyz')))

if not xyz_files:
    raise FileNotFoundError(
        f'No .xyz files found in directory: {os.path.abspath(TRAJ_DIR)}'
    )

traj_names = [os.path.basename(f) for f in xyz_files]
print(f'Found {len(xyz_files)} trajectory file(s):')
for name in traj_names:
    print(f'  {name}')
print()

# ── Process each trajectory ───────────────────────────────────────────────────
# vol_data[frame_index] = {'time_fs': t, 'traj1.xyz': v1, 'traj2.xyz': v2, ...}
# area_data mirrors the same structure.
vol_data  = {}
area_data = {}

for xyz_path, traj_name in zip(xyz_files, traj_names):
    print(f'Processing: {traj_name}')

    u        = mda.Universe(xyz_path)
    solute   = u.atoms[:N_SOLUTE]
    n_frames = len(u.trajectory)

    print(f'  Frames: {n_frames}  |  Atom names (first 5): {u.atoms.names[:5]}')

    for ts in u.trajectory:
        coords  = solute.positions
        syms    = solute.names
        time_fs = round(ts.frame * TIMESTEP_FS, 4)

        try:
            radii = np.array([BONDI[s] for s in syms])
        except KeyError as e:
            raise KeyError(
                f'Element {e} in {traj_name} not in BONDI table. '
                f'Add its Alvarez (2013) radius and rerun.'
            )

        area_ses, vol_ses = run_msms(coords, radii, probe_radius=PROBE_RADIUS)

        # Initialise the row dict the first time we see this frame index
        if ts.frame not in vol_data:
            vol_data[ts.frame]  = {'time_fs': time_fs}
            area_data[ts.frame] = {'time_fs': time_fs}

        vol_data[ts.frame][traj_name]  = round(vol_ses,  4)
        area_data[ts.frame][traj_name] = round(area_ses, 4) if area_ses is not None else ''

        print(f'  frame {ts.frame+1:4d}/{n_frames}  '
              f'time={time_fs:.2f} fs  '
              f'V_SES={vol_ses:.4f} Å³  '
              f'SA_SES={area_ses:.4f} Å²')

    print(f'  Done: {traj_name}\n')

# ── Write CSVs ────────────────────────────────────────────────────────────────
# Column order: time_fs, then each trajectory name in the order they were found
fields = ['time_fs'] + traj_names

# Sort rows by frame index so time runs 0 → end
sorted_frames = sorted(vol_data.keys())

with open(OUT_VOL_CSV, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
    w.writeheader()
    w.writerows([vol_data[i] for i in sorted_frames])

with open(OUT_AREA_CSV, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
    w.writeheader()
    w.writerows([area_data[i] for i in sorted_frames])

print(f'Volume data  → {OUT_VOL_CSV}')
print(f'Area data    → {OUT_AREA_CSV}')
print('All done.')