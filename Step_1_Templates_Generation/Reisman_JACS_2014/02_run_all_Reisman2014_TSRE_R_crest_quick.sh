#!/usr/bin/env bash
set -u

command -v crest >/dev/null 2>&1 || {
    echo "ERROR: crest not found in PATH."
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

    if [[ ! -f "$d/TSRE_R.xyz" ]]; then
        echo "WARNING: $d/TSRE_R.xyz not found. Skipping."
        continue
    fi

    if [[ ! -f "$d/TSRE_R_const.inp" ]]; then
        echo "WARNING: $d/TSRE_R_const.inp not found. Skipping."
        continue
    fi

    if [[ -f "$d/crest_quick.out" ]] && grep -q "CREST terminated normally" "$d/crest_quick.out"; then
        echo "Already completed. Skipping."
        continue
    fi

    (
        cd "$d" || exit 1
        crest TSRE_R.xyz \
          --gfnff \
          --chrg 0 \
          --uhf 1 \
          --cinp TSRE_R_const.inp \
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
echo "All available TSRE_R quick CREST searches processed."
echo "Finished: $(date)"
echo "============================================================"
