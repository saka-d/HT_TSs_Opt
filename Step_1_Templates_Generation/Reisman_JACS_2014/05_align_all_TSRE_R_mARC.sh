#!/usr/bin/env bash
set -u

# Align all accepted mARC conformers for TSRE_R substrates (3a-3n).
# 0-based atom indices:
#   Ni  = 37
#   N42 = 42
#   N43 = 43
# Alignment conditions:
#   1) Ni37 -> origin
#   2) vector N43 -> N42 -> +y
#   3) N42/N43 -> (-x, y, 0)
#
# Run from:
#   ~/HT_TSs_Opt/Step_1_Templates_Generation/Reisman_JACS_2014
#
# Requires:
#   conda activate crest_env
#   numpy installed in the active environment

command -v python >/dev/null 2>&1 || {
    echo "ERROR: python not found in PATH."
    exit 1
}

python - <<'PY'
import numpy as np
from pathlib import Path

NI = 37
N42 = 42
N43 = 43
ROOT = Path.cwd()


def read_xyz(path):
    lines = path.read_text().splitlines()
    nat = int(lines[0])
    comment = lines[1] if len(lines) > 1 else ""

    atoms = []
    coords = []
    for line in lines[2:2 + nat]:
        p = line.split()
        atoms.append(p[0])
        coords.append([float(p[1]), float(p[2]), float(p[3])])

    xyz = np.array(coords, dtype=float)
    if len(xyz) != nat:
        raise ValueError(f"Expected {nat} atoms but read {len(xyz)} atoms")

    return nat, comment, atoms, xyz


def write_xyz(path, nat, comment, atoms, xyz):
    with open(path, "w") as f:
        f.write(f"{nat}\n")
        f.write(f"{comment}\n")
        for atom, r in zip(atoms, xyz):
            f.write(f"{atom:2s} {r[0]:16.10f} {r[1]:16.10f} {r[2]:16.10f}\n")


def rotation_from_a_to_b(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)

    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    if na < 1e-12 or nb < 1e-12:
        raise ValueError("Cannot align a zero-length vector")

    a = a / na
    b = b / nb

    v = np.cross(a, b)
    c = np.clip(np.dot(a, b), -1.0, 1.0)
    s = np.linalg.norm(v)

    if s < 1e-12:
        if c > 0.0:
            return np.eye(3)

        # 180-degree rotation around any axis perpendicular to a
        trial = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(a, trial)) > 0.9:
            trial = np.array([0.0, 0.0, 1.0])

        axis = np.cross(a, trial)
        axis /= np.linalg.norm(axis)
        K = np.array([
            [0.0, -axis[2], axis[1]],
            [axis[2], 0.0, -axis[0]],
            [-axis[1], axis[0], 0.0],
        ])
        return np.eye(3) + 2.0 * (K @ K)

    axis = v / s
    K = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ])

    return np.eye(3) + s * K + (1.0 - c) * (K @ K)


def align_xyz(xyz):
    if xyz.shape[0] <= max(NI, N42, N43):
        raise ValueError("Structure does not contain the required atom indices")

    xyz = xyz.copy()

    # Step 1: Ni37 -> origin
    xyz -= xyz[NI]

    # Step 2: N43 -> N42 -> +y
    vec = xyz[N42] - xyz[N43]
    R = rotation_from_a_to_b(vec, np.array([0.0, 1.0, 0.0]))
    xyz = xyz @ R.T

    # Step 3: rotate around y so N42/N43 lie in xy plane on negative-x side
    xN = 0.5 * (xyz[N42, 0] + xyz[N43, 0])
    zN = 0.5 * (xyz[N42, 2] + xyz[N43, 2])
    rho = np.hypot(xN, zN)

    if rho > 1e-12:
        cosphi = -xN / rho
        sinphi = -zN / rho

        Ry = np.array([
            [ cosphi, 0.0, sinphi],
            [ 0.0,    1.0, 0.0   ],
            [-sinphi, 0.0, cosphi],
        ])
        xyz = xyz @ Ry.T

    return xyz


substrate_dirs = sorted(ROOT.glob("TSRE_R_3[abcdefghijklmn]_*/"))

if not substrate_dirs:
    raise SystemExit(
        "ERROR: No TSRE_R_3a-3n substrate directories found. "
        "Run this script from Reisman_JACS_2014."
    )

n_substrates = 0
n_structures = 0
n_failed = 0

tol = 1e-6

for substrate in substrate_dirs:
    marc_dir = substrate / "mARC_input"
    if not marc_dir.is_dir():
        print(f"SKIP: {substrate.name} (mARC_input not found)")
        continue

    accepted = sorted(marc_dir.glob("*_accepted.xyz"))
    if not accepted:
        print(f"SKIP: {substrate.name} (no *_accepted.xyz files)")
        continue

    n_substrates += 1
    print()
    print("============================================================")
    print(f"Substrate: {substrate.name}")
    print(f"Accepted conformers: {len(accepted)}")
    print("============================================================")

    for path in accepted:
        try:
            nat, comment, atoms, xyz = read_xyz(path)
            aligned = align_xyz(xyz)

            out = path.with_name(path.stem + "_aligned.xyz")
            write_xyz(out, nat, comment, atoms, aligned)

            # Validation
            ni_norm = np.linalg.norm(aligned[NI])
            v = aligned[N42] - aligned[N43]
            x_common = 0.5 * (aligned[N42, 0] + aligned[N43, 0])
            z42 = aligned[N42, 2]
            z43 = aligned[N43, 2]

            ok = (
                ni_norm < tol
                and abs(v[0]) < tol
                and v[1] > 0.0
                and abs(v[2]) < tol
                and abs(z42) < tol
                and abs(z43) < tol
                and x_common < 0.0
            )

            status = "OK" if ok else "CHECK"
            print(f"{status}: {path.name} -> {out.name}")
            print(
                f"      Ni37=({aligned[NI,0]: .6f}, {aligned[NI,1]: .6f}, {aligned[NI,2]: .6f})"
            )
            print(
                f"      N42 =({aligned[N42,0]: .6f}, {aligned[N42,1]: .6f}, {aligned[N42,2]: .6f})"
            )
            print(
                f"      N43 =({aligned[N43,0]: .6f}, {aligned[N43,1]: .6f}, {aligned[N43,2]: .6f})"
            )
            print(
                f"      43->42=({v[0]: .6f}, {v[1]: .6f}, {v[2]: .6f})"
            )

            n_structures += 1
            if not ok:
                n_failed += 1

        except Exception as e:
            n_failed += 1
            print(f"ERROR: {path}: {e}")

print()
print("============================================================")
print(f"Substrates processed: {n_substrates}")
print(f"Structures aligned:   {n_structures}")
print(f"Checks/errors:        {n_failed}")
print("============================================================")
PY
