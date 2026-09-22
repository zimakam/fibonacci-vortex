import numpy as np
import os

def parse_massmodels(path):
    """Parse MassModels_Lelli2016c.mrt — per-radius rotation curves."""
    galaxies = {}
    with open(path, 'r') as f:
        for line in f:
            # Skip headers and empty lines
            if not line or line[0] == '#' or len(line) < 40:
                continue
            # Try fixed-width parse
            try:
                gid = line[0:11].strip()
                if not gid or gid.startswith('Title') or gid.startswith('Author'):
                    continue
                D = float(line[12:18].strip())
                R = float(line[19:25].strip())
                Vobs = float(line[26:32].strip())
                eV = float(line[33:38].strip())
                Vg = float(line[39:45].strip())
                Vd = float(line[46:52].strip())
                Vb = float(line[53:59].strip())
                SBd = float(line[60:67].strip()) if len(line) > 67 else 0.0
                SBb = float(line[68:76].strip()) if len(line) > 76 else 0.0
            except (ValueError, IndexError):
                continue
            galaxies.setdefault(gid, []).append((R, Vobs, eV, Vg, Vd, Vb, SBd, SBb))
    return {k: np.array(v) for k, v in galaxies.items()}

path = os.path.expanduser('~/sparc-vortex-nfw/data/MassModels_Lelli2016c.mrt')
gals = parse_massmodels(path)
print(f"Loaded {len(gals)} galaxies")
print()
# Show 3 examples
for name in list(gals.keys())[:3]:
    print(f"--- {name} ---")
    print(f"  N_points = {len(gals[name])}")
    print(f"  First row: R={gals[name][0,0]:.2f}, Vobs={gals[name][0,1]:.2f}, eV={gals[name][0,2]:.2f}")
    print(f"  Last  row: R={gals[name][-1,0]:.2f}, Vobs={gals[name][-1,1]:.2f}, eV={gals[name][-1,2]:.2f}")
print()
# Show list of galaxies
print("First 10 galaxies:", list(gals.keys())[:10])
