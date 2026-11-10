#!/usr/bin/env bash
set -e; export LC_ALL=C

cd "$(dirname "$0")"

ZINQ="../zig-out/bin/zinq"; ORCA="orca"

TARGETS=(
    "rhf_water" "mp2_water" "mp3_water" "cisd_water" "casscf_water" "dft_b3lyp_water" "dft_pbe_water"
    "rhf_h2"    "mp2_h2"    "mp3_h2"    "cisd_h2"    "casscf_h2"    "dft_b3lyp_h2"    "dft_pbe_h2"
    "rhf_hf"    "mp2_hf"    "mp3_hf"    "cisd_hf"    "casscf_hf"    "dft_b3lyp_hf"    "dft_pbe_hf"
)

if [ $# -gt 0 ]; then
    TARGETS=("$@")
fi

TABLE=()

for t in "${TARGETS[@]}"; do
    method="${t%_*}"; system="${t##*_}"

    echo "Method: $method, System: $system"

    "$ZINQ" "zinq_${t}.json" > "zinq_${t}.out"
    "$ORCA" "orca_${t}.inp"  > "orca_${t}.out"

    ez=$(awk '/(FINAL|MP0).*Eh/ {print $(NF-1)}' "zinq_${t}.out" | tail -n 1)

    eo=$(grep "FINAL SINGLE POINT ENERGY" "orca_${t}.out" | tail -n 1 | awk '{print $5}')

    de=$(awk -v z="$ez" -v o="$eo" 'BEGIN {printf "%.8f", z - o}')

    printf "  ZINQ: %14.8f Eh\n" "$ez"
    printf "  ORCA: %14.8f Eh\n" "$eo"
    printf "  DE:   %14.8f Eh\n" "$de"

    echo ""

    TABLE+=("$(printf "%-12s %-8s %16.8f %16.8f %16.8f" "$method" "$system" "$ez" "$eo" "$de")")
done

printf "%-12s %-8s %16s %16s %16s\n" "Method" "System" "ZINQ (Eh)" "ORCA (Eh)" "DE (Eh)"
printf '%s\n' "------------------------------------------------------------------------"

printf '%s\n' "${TABLE[@]}"
