#!/usr/bin/env bash
set -u

# Generate Gaussian .gjf files for all TSRE_S substrates (3a-3n)
# from mARC-aligned accepted conformers.
#
# Run from:
# ~/HT_TSs_Opt/Step_1_Templates_Generation/Reisman_JACS_2014
#
# Expected per-substrate structure:
#   TSRE_S_3a_*/
#     mARC_input/
#       conf_XXXX_accepted_aligned.xyz
#
# Output:
#   TSRE_S_3a_*/
#     gjf/
#       conf_XXXX_accepted_aligned.gjf
#
# Gaussian settings matched to sub16 resources:
#   %nprocshared=32
#   %mem=120GB
#   WB97XD/def2TZVP
#   SMD(DMA)
#   charge = 0
#   multiplicity = 2
#   nosymm
#
# Two blank lines are added after the final atom.

for d in TSRE_S_3{a..n}_*/
do
    [[ -d "$d" ]] || continue

    echo
    echo "============================================================"
    echo "Substrate: ${d%/}"
    echo "============================================================"

    if [[ ! -d "$d/mARC_input" ]]; then
        echo "WARNING: $d/mARC_input not found. Skipping."
        continue
    fi

    mkdir -p "$d/gjf"

    count=0

    for xyz in "$d"/mARC_input/*_accepted_aligned.xyz
    do
        [[ -f "$xyz" ]] || continue

        base=$(basename "$xyz" .xyz)

        {
            echo "%nprocshared=32"
            echo "%mem=120GB"
            echo "%chk=${base}.chk"
            echo "# wb97xd/def2tzvp scrf=(smd,solvent=n,n-DiMethylAcetamide) nosymm"
            echo
            echo "${base}"
            echo
            echo "0 2"

            tail -n +3 "$xyz"

            echo
            echo

        } > "$d/gjf/${base}.gjf"

        echo "Created: $d/gjf/${base}.gjf"
        count=$((count + 1))
    done

    if [[ $count -eq 0 ]]; then
        echo "WARNING: no *_accepted_aligned.xyz files found."
    else
        echo "Generated $count Gaussian input file(s)."
    fi
done

echo
echo "============================================================"
echo "All available TSRE_S substrates processed."
echo "============================================================"
