#!/usr/bin/env bash
set -u

# Align all mARC-selected TSRE_S conformers (3a-3n).
#
# Alignment convention (0-based atom indices):
#   atom 37 = Ni  -> origin
#   atom 42 = N
#   atom 43 = N
#   vector N43 -> N42 -> +y
#   N42 and N43 lie in the xy plane with x < 0
#
# Run from:
# ~/HT_TSs_Opt/Step_1_Templates_Generation/Reisman_JACS_2014
#
# Expected per substrate:
#   TSRE_S_3a_*/
#     mARC_input/
#       conf_XXXX_accepted.xyz
#
# Output:
#   conf_XXXX_accepted_aligned.xyz
#
# No separate .py file is created.

python - <<'PY'
import numpy as np
from pathlib import Path

NI = 37
N42 = 42
N43 = 43
ROOT = Path(".")

def read_xyz(path):
    lines = path.read_text().splitlines()
    nat = int(lines[0])
    comment = lines[1] if len(lines) > 1 else ""
    atoms = []
    coords = []

    for line in lines[2:2+nat]:
        p = line.split()
        atoms.append(p[0])
        coords.append([float(p[1]), float(p[2]), float(p[3])])

    return nat, comment, atoms, np.array(coords, dtype=float)

def write_xyz(path, nat, comment, atoms, xyz):
    with open(path, "w") as f:
        f.write(f"{nat}\n")
        f.write(f"{comment}\n")
        for atom, r in zip(atoms, xyz):
            f.write(f"{atom:2s} {r[0]:16.10f} {r[1]:16.10f} {r[2]:16.10f}\n")

def rotation_from_a_to_b(a, b):
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)

    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)

    if s < 1e-12:
        if c > 0:
            return np.eye(3)

        trial = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(a, trial)) > 0.9:
            trial = np.array([0.0, 0.0, 1.0])

        axis = np.cross(a, trial)
        axis /= np.linalg.norm(axis)

        K = np.array([
            [0.0, -axis[2], axis[1]],
            [axis[2], 0.0, -axis[0]],
            [-axis[1], axis[0], 0.0]
        ])
        return np.eye(3) + 2.0 * (K @ K)

    axis = v / s
    K = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0]
    ])

    return np.eye(3) + s * K + (1.0 - c) * (K @ K)

for substrate in sorted(ROOT.glob("TSRE_S_3?_*/")):
    marc_dir = substrate / "mARC_input"
    if not marc_dir.is_dir():
        continue

    files = sorted(marc_dir.glob("conf_[0-9][0-9][0-9][0-9]_accepted.xyz"))

    print()
    print("=" * 60)
    print(f"Substrate: {substrate.name}")
    print(f"Accepted conformers found: {len(files)}")
    print("=" * 60)

    for path in files:
        nat, comment, atoms, xyz = read_xyz(path)

        if nat <= max(NI, N42, N43):
            print(f"CHECK: {path.name} has only {nat} atoms")
            continue

        # Step 1: Ni37 -> origin
        xyz = xyz - xyz[NI]

        # Step 2: N43 -> N42 -> +y
        vec = xyz[N42] - xyz[N43]
        R = rotation_from_a_to_b(vec, np.array([0.0, 1.0, 0.0]))
        xyz = xyz @ R.T

        # Step 3: rotate around y so N42/N43 are in xy plane
        # and on the negative-x side.
        xN = 0.5 * (xyz[N42, 0] + xyz[N43, 0])
        zN = 0.5 * (xyz[N42, 2] + xyz[N43, 2])

        rho = np.hypot(xN, zN)

        if rho > 1e-12:
            c = -xN / rho
            s = -zN / rho

            Ry = np.array([
                [ c, 0.0, s],
                [0.0, 1.0, 0.0],
                [-s, 0.0, c]
            ])

            xyz = xyz @ Ry.T

        out = path.with_name(path.stem + "_aligned.xyz")
        write_xyz(out, nat, comment, atoms, xyz)

        vfinal = xyz[N42] - xyz[N43]

        ok = (
            np.linalg.norm(xyz[NI]) < 1e-6
            and abs(xyz[N42, 2]) < 1e-6
            and abs(xyz[N43, 2]) < 1e-6
            and xyz[N42, 0] < 0
            and xyz[N43, 0] < 0
            and abs(vfinal[0]) < 1e-6
            and vfinal[1] > 0
            and abs(vfinal[2]) < 1e-6
        )

        status = "OK" if ok else "CHECK"
        print(f"{status}: {path.name} -> {out.name}")
        print(f"  Ni37       = {xyz[NI]}")
        print(f"  N42        = {xyz[N42]}")
        print(f"  N43        = {xyz[N43]}")
        print(f"  N43->N42   = {vfinal}")

print()
print("=" * 60)
print("All available TSRE_S accepted conformers processed.")
print("=" * 60)
PY
