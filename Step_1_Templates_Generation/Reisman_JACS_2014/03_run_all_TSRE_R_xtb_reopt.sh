#!/usr/bin/env bash
set -u

# Re-optimize all CREST conformers for all Reisman_JACS_2014 substrates
# using constrained GFN2-xTB.
#
# Expected folder layout:
#   TSRE_R_3a_H_Me/
#     crest_conformers.xyz
#     TSRE_R_const.inp
#   TSRE_R_3b_4Me_Me/
#     crest_conformers.xyz
#     TSRE_R_const.inp
#   ...
#
# Run this script from:
# ~/HT_TSs_Opt/Step_1_Templates_Generation/Reisman_JACS_2014
#
# Recommended:
#   conda activate crest_env
#   chmod +x run_all_xtb_reopt.sh
#   ./run_all_xtb_reopt.sh

command -v xtb >/dev/null 2>&1 || {
    echo "ERROR: xtb not found in PATH."
    echo "Activate crest_env first: conda activate crest_env"
    exit 1
}

for d in TSRE_R_3{a..n}_*/
do
    [[ -d "$d" ]] || continue

    echo
    echo "============================================================"
    echo "Substrate: ${d%/}"
    echo "Started:   $(date)"
    echo "============================================================"

    if [[ ! -f "$d/crest_conformers.xyz" ]]; then
        echo "WARNING: $d/crest_conformers.xyz not found. Skipping."
        continue
    fi

    if [[ ! -f "$d/TSRE_R_const.inp" ]]; then
        echo "WARNING: $d/TSRE_R_const.inp not found. Skipping."
        continue
    fi

    cd "$d" || exit 1

    # Split multi-XYZ file only if conformer directory does not yet exist.
    if [[ ! -d conformers ]]; then
        echo "Splitting crest_conformers.xyz ..."

        mkdir -p conformers

        natoms=$(awk 'NR==1 {print $1}' crest_conformers.xyz)
        block=$((natoms + 2))

        awk -v n="$block" '
        {
            file=sprintf("conformers/conf_%04d.xyz", int((NR-1)/n)+1)
            print > file
        }
        ' crest_conformers.xyz
    fi

    nconf=$(find conformers -maxdepth 1 -name 'conf_*.xyz' -type f | wc -l)
    echo "Conformers found: $nconf"

    mkdir -p reopt

    for f in conformers/conf_*.xyz
    do
        [[ -f "$f" ]] || continue

        name=$(basename "$f" .xyz)

        # Restart-friendly: skip if optimized geometry already exists.
        if [[ -f "reopt/$name/xtbopt.xyz" ]]; then
            echo "SKIP: $name (xtbopt.xyz already exists)"
            continue
        fi

        mkdir -p "reopt/$name"
        cp "$f" "reopt/$name/input.xyz"
        cp TSRE_R_const.inp "reopt/$name/"

        echo "RUN:  $name"

        (
            cd "reopt/$name" || exit 1

            xtb input.xyz --gfn 2 --chrg 0 --uhf 1 \
                --input TSRE_R_const.inp --opt \
                > xtb.out 2>&1
        )

        status=$?

        if [[ $status -eq 0 && -f "reopt/$name/xtbopt.xyz" ]]; then
            echo "DONE: $name"
        else
            echo "WARNING: optimization may have failed: ${d%/}/$name"
            echo "         Check ${d%/}/reopt/$name/xtb.out"
        fi
    done

    cd ..

    echo "Finished substrate: ${d%/}"
    echo "Time: $(date)"
done

echo
echo "============================================================"
echo "All available substrates processed."
echo "Finished: $(date)"
echo "============================================================"
