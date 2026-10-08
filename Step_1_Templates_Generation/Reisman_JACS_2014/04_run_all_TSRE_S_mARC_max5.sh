#!/usr/bin/env bash
set -u

# Run mARC for all TSRE_S substrates (3a-3n)
# using constrained GFN2-xTB optimized conformers.
#
# Selection settings:
#   RMSD metric
#   maximum target cluster count: 5
#   10 kcal/mol energy window
#   choose minimum-energy conformer from each cluster
#   verbosity level 2
#
# Run from:
# ~/HT_TSs_Opt/Step_1_Templates_Generation/Reisman_JACS_2014

command -v python >/dev/null 2>&1 || {
    echo "ERROR: python not found in PATH."
    exit 1
}

python -m navicat_marc --help >/dev/null 2>&1 || {
    echo "ERROR: navicat_marc is not available in the current Python environment."
    echo "Activate crest_env first."
    exit 1
}

for d in TSRE_S_3{a..n}_*/
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

    mkdir -p "$d/mARC_input"

    # Remove old mARC XYZ files so old accepted/rejected files cannot be re-read.
    rm -f "$d"/mARC_input/conf_*.xyz

    count=0

    for c in "$d"/reopt/conf_*/
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
        } > "$d/mARC_input/${name}.xyz"

        count=$((count + 1))
    done

    echo "mARC input structures: $count"

    if [[ $count -eq 0 ]]; then
        echo "WARNING: no valid structures found. Skipping mARC."
        continue
    fi

    (
        cd "$d/mARC_input" || exit 1

        echo "Running mARC ..."
        python -m navicat_marc \
          -i conf_[0-9][0-9][0-9][0-9].xyz \
          -m rmsd \
          -n 5 \
          -ewin 10 \
          -mine \
          -v 2 \
          > marc.out 2>&1
    )

    status=$?

    if [[ $status -eq 0 ]]; then
        accepted=$(find "$d/mARC_input" -maxdepth 1 -type f -name '*_accepted.xyz' | wc -l)
        echo "mARC completed successfully for ${d%/}"
        echo "Accepted conformers: $accepted"
    else
        echo "WARNING: mARC failed for ${d%/}"
        echo "Check: $d/mARC_input/marc.out"
    fi

    echo "Finished: $(date)"
done

echo
echo "============================================================"
echo "All available TSRE_S substrates processed."
echo "Finished: $(date)"
echo "============================================================"
