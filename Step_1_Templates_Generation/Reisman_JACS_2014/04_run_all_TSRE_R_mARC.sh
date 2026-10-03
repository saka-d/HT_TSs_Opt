#!/usr/bin/env bash
set -u

# Run preliminary mARC conformer selection for all TSRE_R substrates (3a-3n)
# using the current constrained GFN2-xTB optimized structures.
#
# This is a preliminary reduction step, NOT the final Sigman workflow.
#
# Expected structure:
#   TSRE_R_3a_H_Me/
#     reopt/conf_0001/xtbopt.xyz
#     reopt/conf_0001/xtb.out
#     reopt/conf_0002/...
#   ...
#
# Run from:
# ~/HT_TSs_Opt/Step_1_Templates_Generation/Reisman_JACS_2014
#
# Requires:
#   conda activate crest_env
#   navicat-marc installed in the active environment

command -v python >/dev/null 2>&1 || {
    echo "ERROR: python not found in PATH."
    exit 1
}

python -m navicat_marc --help >/dev/null 2>&1 || {
    echo "ERROR: navicat_marc is not available in the current Python environment."
    echo "Activate crest_env first."
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

    if [[ ! -d "$d/reopt" ]]; then
        echo "WARNING: $d/reopt not found. Skipping."
        continue
    fi

    cd "$d" || exit 1

    mkdir -p mARC_input
    rm -f mARC_input/conf_*.xyz

    count=0

    for c in reopt/conf_*/
    do
        [[ -d "$c" ]] || continue

        name=$(basename "$c")

        if [[ ! -f "$c/xtbopt.xyz" ]]; then
            echo "SKIP: $name (xtbopt.xyz missing)"
            continue
        fi

        if [[ ! -f "$c/xtb.out" ]]; then
            echo "SKIP: $name (xtb.out missing)"
            continue
        fi

        energy=$(grep "TOTAL ENERGY" "$c/xtb.out" | tail -1 | awk '{print $4}')

        if [[ -z "$energy" ]]; then
            echo "SKIP: $name (TOTAL ENERGY not found)"
            continue
        fi

        {
            head -n 1 "$c/xtbopt.xyz"
            echo "$energy"
            tail -n +3 "$c/xtbopt.xyz"
        } > "mARC_input/${name}.xyz"

        count=$((count + 1))
    done

    echo "mARC input structures: $count"

    if [[ $count -eq 0 ]]; then
        echo "WARNING: no valid structures found. Skipping mARC."
        cd ..
        continue
    fi

    cd mARC_input || exit 1

    echo "Running mARC ..."
    python -m navicat_marc         -i conf_*.xyz         -m rmsd         -ewin 10         -mine         -v 2         > marc.out 2>&1

    status=$?

    if [[ $status -eq 0 ]]; then
        echo "mARC completed successfully for ${d%/}"
    else
        echo "WARNING: mARC failed for ${d%/}"
        echo "Check: ${d%/}/mARC_input/marc.out"
    fi

    cd ../..

    echo "Finished: $(date)"
done

echo
echo "============================================================"
echo "All available TSRE_R substrates processed."
echo "Finished: $(date)"
echo "============================================================"
