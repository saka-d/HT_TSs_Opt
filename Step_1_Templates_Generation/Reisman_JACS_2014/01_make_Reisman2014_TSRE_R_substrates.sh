#!/usr/bin/env bash
set -euo pipefail

# Generate Reisman_JACS_2014 TSRE_R substrate variants
# Each substrate folder uses the common structure filename: TSRE_R.xyz
# Requires:
#   - AaronTools substitute.py available in PATH
#   - TSRE_R.xyz in the current directory
#   - TSRE_R_const.inp in the current directory
#
# Run from:
# ~/HT_TSs_Opt/Step_1_Templates_Generation/Reisman_JACS_2014

BASE_XYZ="TSRE_R.xyz"
CONST_INP="TSRE_R_const.inp"

if [[ ! -f "$BASE_XYZ" ]]; then
    echo "ERROR: $BASE_XYZ not found in $(pwd)"
    exit 1
fi

if [[ ! -f "$CONST_INP" ]]; then
    echo "ERROR: $CONST_INP not found in $(pwd)"
    exit 1
fi

command -v substitute.py >/dev/null 2>&1 || {
    echo "ERROR: substitute.py not found in PATH"
    exit 1
}

make_parent() {
    local dir="$1"
    mkdir -p "$dir"
    cp "$BASE_XYZ" "$dir/TSRE_R.xyz"
    cp "$CONST_INP" "$dir/"
}

make_sub() {
    local dir="$1"
    local substitution="$2"

    mkdir -p "$dir"
    substitute.py "$BASE_XYZ" -s "$substitution" -o "$dir/TSRE_R.xyz"
    cp "$CONST_INP" "$dir/"
}

# 3a: R1 = H, R2 = Me (parent)
make_parent \
    "TSRE_R_3a_H_Me"

# Aromatic R1 substitutions
make_sub \
    "TSRE_R_3b_4Me_Me" \
    "34=Me"

make_sub \
    "TSRE_R_3c_3Me_Me" \
    "33=Me"

make_sub \
    "TSRE_R_3d_2Me_Me" \
    "32=Me"

make_sub \
    "TSRE_R_3e_4OMe_Me" \
    "34=OMe"

make_sub \
    "TSRE_R_3f_4F_Me" \
    "34=F"

make_sub \
    "TSRE_R_3g_4Cl_Me" \
    "34=Cl"

make_sub \
    "TSRE_R_3h_4Br_Me" \
    "34=Br"

make_sub \
    "TSRE_R_3i_4OCF3_Me" \
    "34=OCF3"

# R2 substitutions
make_sub \
    "TSRE_R_3j_H_Et" \
    "22=Et"

# Explicit radical/attachment atom is required for SMILES substituents in AaronTools.
make_sub \
    "TSRE_R_3k_H_Bn" \
    "22=smiles:[C.]c1ccccc1"

make_sub \
    "TSRE_R_3l_H_4pentenyl" \
    "22=smiles:[C.]CCC=C"

make_sub \
    "TSRE_R_3m_H_2hydroxyethyl" \
    "22=smiles:[C.]CO"

make_sub \
    "TSRE_R_3n_H_2chloroethyl" \
    "22=smiles:[C.]CCl"

echo
echo "Generation complete."
echo

# Summary: show each generated XYZ and atom count.
for d in TSRE_R_3*/; do
    xyz="$d/TSRE_R.xyz"
    if [[ -n "$xyz" ]]; then
        atoms=$(head -n 1 "$xyz")
        printf "%-35s atoms=%s\n" "${d%/}" "$atoms"
    fi
done

echo
echo "Each folder also contains: $CONST_INP"
