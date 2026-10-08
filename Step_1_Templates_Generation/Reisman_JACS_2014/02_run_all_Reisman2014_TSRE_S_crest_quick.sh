#!/usr/bin/env bash
set -u

# Run CREST quick GFN-FF conformer searches for all TSRE_S substrates (3a-3n)
#
# Run from:
# ~/HT_TSs_Opt/Step_1_Templates_Generation/Reisman_JACS_2014
#
# Each substrate folder is expected to contain:
#   TSRE_S.xyz
#   TSRE_S_const.inp
#
# CREST settings:
#   GFN-FF
#   charge = 0
#   UHF = 1
#   constraints from TSRE_S_const.inp
#   no reference-topology check
#   quick search mode
#
# Restart behavior:
#   If crest_quick.out already contains "CREST terminated normally",
#   that substrate is skipped.

command -v crest >/dev/null 2>&1 || {
    echo "ERROR: crest not found in PATH."
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

    if [[ ! -f "$d/TSRE_S.xyz" ]]; then
        echo "WARNING: $d/TSRE_S.xyz not found. Skipping."
        continue
    fi

    if [[ ! -f "$d/TSRE_S_const.inp" ]]; then
        echo "WARNING: $d/TSRE_S_const.inp not found. Skipping."
        continue
    fi

    if [[ -f "$d/crest_quick.out" ]] && grep -q "CREST terminated normally" "$d/crest_quick.out"; then
        echo "Already completed. Skipping."
        continue
    fi

    (
        cd "$d" || exit 1

        crest TSRE_S.xyz \
          --gfnff \
          --chrg 0 \
          --uhf 1 \
          --cinp TSRE_S_const.inp \
          --noreftopo \
          -quick \
          > crest_quick.out 2>&1
    )

    status=$?

    if [[ $status -eq 0 ]] && grep -q "CREST terminated normally" "$d/crest_quick.out"; then
        echo "CREST completed successfully for ${d%/}"
    else
        echo "WARNING: CREST may have failed for ${d%/}"
        echo "Check: $d/crest_quick.out"
    fi

    echo "Finished:  $(date)"
done

echo
echo "============================================================"
echo "All available TSRE_S quick CREST searches processed."
echo "Finished: $(date)"
echo "============================================================"
