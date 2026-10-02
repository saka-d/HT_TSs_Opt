#!/usr/bin/env bash
set -euo pipefail

# Generate Reisman_JACS_2014 TSRE_S substrate variants
#
# Requires:
#   - AaronTools substitute.py available in PATH
#   - TSRE_S.xyz in the current directory
#   - TSRE_S_const.inp in the current directory
#
# Run from:
# ~/HT_TSs_Opt/Step_1_Templates_Generation/Reisman_JACS_2014
#
# Output:
#   TSRE_S_3a_H_Me/TSRE_S.xyz
#   TSRE_S_3b_4Me_Me/TSRE_S.xyz
#   ...
#   TSRE_S_3n_H_2chloroethyl/TSRE_S.xyz
#
# Each folder also receives a copy of TSRE_S_const.inp.

BASE_XYZ="TSRE_S.xyz"
CONST_INP="TSRE_S_const.inp"
COMMON_XYZ="TSRE_S.xyz"

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
    echo "Activate the environment containing AaronTools first."
    exit 1
}

make_parent() {
    local dir="$1"

    mkdir -p "$dir"
    cp "$BASE_XYZ" "$dir/$COMMON_XYZ"
    cp "$CONST_INP" "$dir/"
}

make_sub() {
    local dir="$1"
    local substitution="$2"

    mkdir -p "$dir"
    substitute.py "$BASE_XYZ" -s "$substitution" -o "$dir/$COMMON_XYZ"
    cp "$CONST_INP" "$dir/"
}

# 3a: R1 = H, R2 = Me (parent)
make_parent "TSRE_S_3a_H_Me"

# Aromatic R1 substitutions
make_sub "TSRE_S_3b_4Me_Me"   "34=Me"
make_sub "TSRE_S_3c_3Me_Me"   "33=Me"
make_sub "TSRE_S_3d_2Me_Me"   "32=Me"
make_sub "TSRE_S_3e_4OMe_Me"  "34=OMe"
make_sub "TSRE_S_3f_4F_Me"    "34=F"
make_sub "TSRE_S_3g_4Cl_Me"   "34=Cl"
make_sub "TSRE_S_3h_4Br_Me"   "34=Br"
make_sub "TSRE_S_3i_4OCF3_Me" "34=OCF3"

# R2 substitutions
make_sub "TSRE_S_3j_H_Et" "22=Et"

# Explicit radical/attachment atom required for SMILES substituents in AaronTools
make_sub "TSRE_S_3k_H_Bn"              "22=smiles:[C.]c1ccccc1"
make_sub "TSRE_S_3l_H_4pentenyl"       "22=smiles:[C.]CCC=C"
make_sub "TSRE_S_3m_H_2hydroxyethyl"   "22=smiles:[C.]CO"
make_sub "TSRE_S_3n_H_2chloroethyl"    "22=smiles:[C.]CCl"

echo
echo "Generation complete."
echo

# Summary
for d in TSRE_S_3*/; do
    [[ -d "$d" ]] || continue

    if [[ -f "$d/$COMMON_XYZ" ]]; then
        atoms=$(head -n 1 "$d/$COMMON_XYZ")
        printf "%-35s atoms=%s\n" "${d%/}" "$atoms"
    else
        printf "%-35s %s\n" "${d%/}" "ERROR: TSRE_S.xyz missing"
    fi
done

echo
echo "Each folder also contains: $CONST_INP"
