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

    input_dir="$d/mARC_input"
    output_dir="$d/gjf"

    if [ ! -d "$input_dir" ]; then
        echo "SKIP: $input_dir not found"
        continue
    fi

    mkdir -p "$output_dir"

    count=0

    for xyz in "$input_dir"/*_accepted_aligned.xyz
    do
        [ -f "$xyz" ] || continue

        base=$(basename "$xyz" .xyz)
        gjf="$output_dir/${base}.gjf"

        {
            echo "%nprocshared=32"
            echo "%mem=120GB"
            echo "%chk=${base}.chk"
            echo "# wb97xd/def2tzvp scrf=(smd,solvent=n,n-DiMethylAcetamide) nosymm"
            echo
            echo "$base"
            echo
            echo "0 2"
            tail -n +3 "$xyz"
            echo
            echo
        } > "$gjf"

        echo "Created: $gjf"
        count=$((count+1))
    done

    echo "Generated $count Gaussian input file(s) in $output_dir"
done

echo "All R-series Gaussian .gjf files generated."
