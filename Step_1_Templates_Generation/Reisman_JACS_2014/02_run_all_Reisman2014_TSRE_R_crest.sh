#!/usr/bin/env bash
set -u

# Run CREST/GFN-FF sequentially for all Reisman_JACS_2014 TSRE_R substrate folders.
#
# Expected directory structure:
#   TSRE_R_3a_H_Me/
#     TSRE_R.xyz
#     TSRE_R_const.inp
#   TSRE_R_3b_4Me_Me/
#     TSRE_R.xyz
#     TSRE_R_const.inp
#   ...
#
# Run this script from:
# ~/HT_TSs_Opt/Step_1_Templates_Generation/Reisman_JACS_2014

XYZ="TSRE_R.xyz"
CINP="TSRE_R_const.inp"

command -v crest >/dev/null 2>&1 || {
    echo "ERROR: crest not found in PATH."
    echo "Activate crest_env first: conda activate crest_env"
    exit 1
}

for d in TSRE_R_3{a..n}_*/; do
    # If a particular folder does not exist, skip it.
    [[ -d "$d" ]] || continue

    echo
    echo "============================================================"
    echo "Running: ${d%/}"
    echo "Started: $(date)"
    echo "============================================================"

    if [[ ! -f "$d/$XYZ" ]]; then
        echo "WARNING: $d/$XYZ not found. Skipping."
        continue
    fi

    if [[ ! -f "$d/$CINP" ]]; then
        echo "WARNING: $d/$CINP not found. Skipping."
        continue
    fi

    # Resume-friendly behavior:
    # Skip folders whose previous crest.out indicates normal termination.
    if [[ -f "$d/crest.out" ]] && grep -q "CREST terminated normally" "$d/crest.out"; then
        echo "Already completed successfully. Skipping ${d%/}."
        continue
    fi

    (
        cd "$d" || exit 1

        crest "$XYZ" --gfnff --chrg 0 --uhf 1 \
          --cinp "$CINP" --noreftopo \
          > crest.out 2>&1
    )

    status=$?

    if [[ $status -eq 0 ]] && grep -q "CREST terminated normally" "$d/crest.out"; then
        echo "Completed successfully: ${d%/}"
    else
        echo "WARNING: CREST may have failed in ${d%/}"
        echo "Check: $d/crest.out"
    fi

    echo "Finished: $(date)"
done

echo
echo "============================================================"
echo "All available substrate folders processed."
echo "Finished: $(date)"
echo "============================================================"
