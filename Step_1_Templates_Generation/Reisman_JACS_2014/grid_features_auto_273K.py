#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

KB_HARTREE_PER_K = 3.166811563e-6
DEFAULT_T = 273.15
GRID_STEP = 2.0
ELECTRONIC_DENSITY_CLIP = 1e-3
ESP_DENSITY_MASK_THRESHOLD = 1e-2

SCF_RE = re.compile(
    r"SCF Done:\s+E\([^)]+\)\s*=\s*([-+]?\d+\.\d+(?:[DEde][-+]?\d+)?)"
)


def read_sp_energy_hartree(log_path: Path) -> float:
    text = log_path.read_text(errors="replace")
    matches = SCF_RE.findall(text)
    if not matches:
        raise ValueError(f"No SCF Done energy found in {log_path}")
    return float(matches[-1].replace("D", "E").replace("d", "e"))


def read_gaussian_cube(cube_path: Path):
    with cube_path.open("r", encoding="utf-8", errors="replace") as f:
        f.readline()
        f.readline()

        third = f.readline().split()
        if len(third) < 4:
            raise ValueError(f"Invalid cube header in {cube_path}")

        n_atom = int(third[0])
        origin = np.array(third[1:4], dtype=float)

        size = []
        axis = []
        for _ in range(3):
            parts = f.readline().split()
            if len(parts) < 4:
                raise ValueError(f"Invalid grid header in {cube_path}")
            size.append(abs(int(parts[0])))
            axis.append([float(parts[1]), float(parts[2]), float(parts[3])])

        size = np.array(size, dtype=int)
        axis = np.array(axis, dtype=float)

        for _ in range(abs(n_atom)):
            f.readline()

        values = np.fromstring(f.read(), dtype=float, sep=" ")

    n_values = int(np.prod(size))
    if values.size < n_values:
        raise ValueError(
            f"{cube_path}: found {values.size} values, expected {n_values}"
        )
    if values.size > n_values:
        values = values[-n_values:]

    ijk = np.indices(tuple(size)).reshape(3, -1).T
    coords = ijk @ axis + origin

    meta = {"origin": origin, "size": size, "axis": axis}
    return coords, values, meta


def same_grid(a, b):
    return (
        np.array_equal(a["size"], b["size"])
        and np.allclose(a["origin"], b["origin"])
        and np.allclose(a["axis"], b["axis"])
    )


def coarse_grid(coords, density, esp):
    rho = np.asarray(density)

    electrostatic = np.asarray(esp) * np.where(
        rho < ESP_DENSITY_MASK_THRESHOLD,
        ESP_DENSITY_MASK_THRESHOLD - rho,
        0.0,
    )

    electronic = np.where(
        rho < ELECTRONIC_DENSITY_CLIP,
        rho,
        ELECTRONIC_DENSITY_CLIP,
    )

    xyz = coords / GRID_STEP
    xyz = np.where(xyz > 0, np.ceil(xyz), np.floor(xyz)).astype(int)

    df = pd.DataFrame(
        {
            "x": xyz[:, 0],
            "y": xyz[:, 1],
            "z": xyz[:, 2],
            "electronic": electronic,
            "electrostatic": electrostatic,
        }
    )

    return (
        df.groupby(["x", "y", "z"], as_index=False)[
            ["electronic", "electrostatic"]
        ]
        .sum()
        .sort_values(["x", "y", "z"])
        .reset_index(drop=True)
    )


def series_from_grid(grid, folded=False):
    g = grid.copy()
    suffix = "fold" if folded else "unfold"

    if folded:
        g["y"] = g["y"].abs()
        g = g.groupby(["x", "y", "z"], as_index=False)[
            ["electronic", "electrostatic"]
        ].sum()

    values = {}
    for row in g.itertuples(index=False):
        values[
            f"electronic_{suffix} {int(row.x)} {int(row.y)} {int(row.z)}"
        ] = row.electronic
        values[
            f"electrostatic_{suffix} {int(row.x)} {int(row.y)} {int(row.z)}"
        ] = row.electrostatic

    return pd.Series(values, dtype=float)


def process_sp_dir(sp_dir: Path, temperature: float):
    sp_dir = sp_dir.resolve()
    out_dir = sp_dir / "grid_features"
    out_dir.mkdir(exist_ok=True)

    logs = sorted(sp_dir.glob("conf_*_accepted_aligned.log"))
    if not logs:
        raise FileNotFoundError(
            f"No conf_*_accepted_aligned.log files found in {sp_dir}"
        )

    conformers = []
    failures = []

    for log in logs:
        stem = log.stem
        density_cube = sp_dir / f"{stem}_density.cube"
        esp_cube = sp_dir / f"{stem}_ESP.cube"

        try:
            if not density_cube.exists():
                raise FileNotFoundError(density_cube)
            if not esp_cube.exists():
                raise FileNotFoundError(esp_cube)

            energy = read_sp_energy_hartree(log)

            coords_d, density, meta_d = read_gaussian_cube(density_cube)
            coords_e, esp, meta_e = read_gaussian_cube(esp_cube)

            if not same_grid(meta_d, meta_e):
                raise ValueError("Density and ESP cube grids do not match")
            if not np.allclose(coords_d, coords_e):
                raise ValueError("Density and ESP coordinates do not match")

            grid = coarse_grid(coords_d, density, esp)

            conformers.append(
                {
                    "name": stem,
                    "energy": energy,
                    "grid": grid,
                }
            )

            print(f"SUCCESS  {stem}  E = {energy:.12f} Eh")

        except Exception as exc:
            failures.append((stem, str(exc)))
            print(f"FAILURE  {stem}: {exc}")

    if failures:
        details = "\n".join(f"{name}: {msg}" for name, msg in failures)
        raise RuntimeError(
            "At least one conformer failed. "
            "Refusing to renormalize only successful conformers.\n"
            + details
        )

    energies = np.array([c["energy"] for c in conformers], dtype=float)
    delta_e = energies - energies.min()

    exponent = -delta_e / (KB_HARTREE_PER_K * temperature)
    exponent = np.clip(exponent, -700.0, 0.0)
    boltzmann = np.exp(exponent)
    boltzmann /= boltzmann.sum()

    weights_df = pd.DataFrame(
        {
            "conformer": [c["name"] for c in conformers],
            "energy_hartree": energies,
            "delta_E_hartree": delta_e,
            "delta_E_kcal_mol": delta_e * 627.509474,
            "boltzmann_weight": boltzmann,
        }
    ).sort_values("energy_hartree")

    weights_df.to_csv(out_dir / "boltzmann_weights.csv", index=False)

    long_tables = []
    weighted_grids = []
    wide_rows = {}

    for conformer, weight in zip(conformers, boltzmann):
        grid = conformer["grid"].copy()

        long_grid = grid.copy()
        long_grid.insert(0, "conformer", conformer["name"])
        long_grid["energy_hartree"] = conformer["energy"]
        long_grid["boltzmann_weight"] = weight
        long_tables.append(long_grid)

        weighted = grid.copy()
        weighted[["electronic", "electrostatic"]] *= weight
        weighted_grids.append(weighted)

        wide_rows[conformer["name"]] = pd.concat(
            [
                series_from_grid(grid, folded=False),
                series_from_grid(grid, folded=True),
            ]
        )

    pd.concat(long_tables, ignore_index=True).to_csv(
        out_dir / "conformer_features_long.csv",
        index=False,
    )

    weighted_grid = (
        pd.concat(weighted_grids, ignore_index=True)
        .groupby(["x", "y", "z"], as_index=False)[
            ["electronic", "electrostatic"]
        ]
        .sum()
        .sort_values(["x", "y", "z"])
        .reset_index(drop=True)
    )

    weighted_grid.to_csv(
        out_dir / "boltzmann_weighted_grid.csv",
        index=False,
    )

    wide_rows["BOLTZMANN_WEIGHTED"] = pd.concat(
        [
            series_from_grid(weighted_grid, folded=False),
            series_from_grid(weighted_grid, folded=True),
        ]
    )

    wide = pd.DataFrame(wide_rows).T.fillna(0.0)
    wide.index.name = "conformer"
    wide.to_csv(out_dir / "features_wide.csv")

    print()
    print(f"Temperature: {temperature:.2f} K")
    print(f"Conformers: {len(conformers)}")
    print(f"Output directory: {out_dir}")
    print("Created:")
    print("  boltzmann_weights.csv")
    print("  conformer_features_long.csv")
    print("  boltzmann_weighted_grid.csv")
    print("  features_wide.csv")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Calculate grid features using only the current directory to decide scope. "
            "Run inside an sp directory for one substrate, inside a substrate directory "
            "for its sp directory, or inside Reisman_JACS_2014 to process all substrates."
        )
    )
    parser.add_argument(
        "-T",
        "--temperature",
        type=float,
        default=DEFAULT_T,
        help="Temperature in K for Boltzmann weighting",
    )
    args = parser.parse_args()

    cwd = Path.cwd().resolve()

    # Case 1: current directory itself is an SP directory.
    if list(cwd.glob("conf_*_accepted_aligned.log")):
        print(f"Mode: single SP directory")
        print(f"Target: {cwd}")
        process_sp_dir(cwd, args.temperature)
        return

    # Case 2: current directory is one substrate directory and contains sp/.
    local_sp = cwd / "sp"
    if local_sp.is_dir() and list(local_sp.glob("conf_*_accepted_aligned.log")):
        print(f"Mode: single substrate")
        print(f"Target: {local_sp}")
        process_sp_dir(local_sp, args.temperature)
        return

    # Case 3: current directory is the project root containing TSRE_R/S_* directories.
    sp_dirs = sorted(
        {
            p
            for pattern in ("TSRE_R_3*/sp", "TSRE_S_3*/sp")
            for p in cwd.glob(pattern)
            if p.is_dir() and list(p.glob("conf_*_accepted_aligned.log"))
        }
    )

    if not sp_dirs:
        raise FileNotFoundError(
            "Could not determine calculation scope from the current directory.\n"
            "Run this script from one of these locations:\n"
            "  1. TSRE_*_3*/sp        -> calculate that SP directory only\n"
            "  2. TSRE_*_3*           -> calculate its sp/ directory only\n"
            "  3. Reisman_JACS_2014   -> calculate all available R/S substrate sp/ directories"
        )

    print("Mode: all substrates found under current directory")
    print(f"Root: {cwd}")
    print(f"Found {len(sp_dirs)} SP directories.")

    failures = []

    for i, sp_dir in enumerate(sp_dirs, 1):
        print()
        print("=" * 72)
        print(f"[{i}/{len(sp_dirs)}] {sp_dir.parent.name}")
        print("=" * 72)

        try:
            process_sp_dir(sp_dir, args.temperature)
        except Exception as exc:
            failures.append((sp_dir, str(exc)))
            print(f"FAILED: {exc}")

    print()
    print("=" * 72)
    print(f"Finished: {len(sp_dirs) - len(failures)}/{len(sp_dirs)} successful")

    if failures:
        print("Failed directories:")
        for sp_dir, message in failures:
            print(f"  {sp_dir}: {message}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
