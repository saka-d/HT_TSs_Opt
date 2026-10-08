#!/bin/bash
set -e

BASE_DIR="$HOME/HT_TSs_Opt/Step_1_Templates_Generation/Reisman_JACS_2014"
cd "$BASE_DIR"

for d in TSRE_R_3{a..n}_*
do
    [ -d "$d" ] || continue

    echo "========================================"
    echo "Processing: $d"
    echo "========================================"

    reopt_dir="$d/reopt"
    marc_dir="$d/mARC_input"

    if [ ! -d "$reopt_dir" ]; then
        echo "SKIP: $reopt_dir not found"
        continue
    fi

    mkdir -p "$marc_dir"

    rm -f "$marc_dir"/conf_[0-9][0-9][0-9][0-9].xyz
    rm -f "$marc_dir"/conf_[0-9][0-9][0-9][0-9]_accepted.xyz
    rm -f "$marc_dir"/marc.out

    n=0

    for c in "$reopt_dir"/conf_[0-9][0-9][0-9][0-9]
    do
        [ -d "$c" ] || continue
        [ -f "$c/xtbopt.xyz" ] || continue
        [ -f "$c/xtb.out" ] || continue

        energy=$(grep "TOTAL ENERGY" "$c/xtb.out" | tail -1 | awk '{print $4}')

        if [ -z "$energy" ]; then
            echo "WARNING: no TOTAL ENERGY found in $c/xtb.out"
            continue
        fi

        name=$(basename "$c")

        {
            head -n 1 "$c/xtbopt.xyz"
            echo "$energy"
            tail -n +3 "$c/xtbopt.xyz"
        } > "$marc_dir/${name}.xyz"

        n=$((n+1))
    done

    echo "Prepared $n structures for mARC"

    if [ "$n" -eq 0 ]; then
        echo "SKIP: no valid xTB-reoptimized structures"
        continue
    fi

    (
        cd "$marc_dir"
        python -m navicat_marc \
            -i conf_[0-9][0-9][0-9][0-9].xyz \
            -m rmsd \
            -n 5 \
            -ewin 10 \
            -mine \
            -v 2 \
            > marc.out 2>&1
    )

    accepted=$(find "$marc_dir" -maxdepth 1 -type f -name 'conf_[0-9][0-9][0-9][0-9]_accepted.xyz' | wc -l)
    echo "Accepted representatives: $accepted"
done

echo "All R-series mARC jobs finished."
