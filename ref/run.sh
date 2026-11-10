#!/usr/bin/env bash
set -e; export LC_ALL=C

cd "$(dirname "$0")"

ZINQ="../zig-out/bin/zinq"; ORCA="orca"

METHODS=("rhf" "mp2" "mp3" "cisd" "casscf" "dft_b3lyp" "dft_pbe")

MOLECULES=($(grep -v '^[[:space:]]*#' molecules.txt)); BASES=($(grep -v '^[[:space:]]*#' bases.txt))

for mol in "${MOLECULES[@]}"; do
    xyz=$(find ../example/molecule -maxdepth 1 -iname "${mol}.xyz" | head -n 1)

    m="${mol,,}"

    for bas in "${BASES[@]}"; do
        b="${bas,,}"

        for t in template/orca_*.inp template/zinq_*.json; do
            out="${t#template/}"; out="${out%.*}_${m}_${b}.${out##*.}"

            sed -e "s|{{XYZ}}|$xyz|g" -e "s|{{BASIS}}|$bas|g" "$t" > "$out"
        done
    done
done

TARGETS=()

for mol in "${MOLECULES[@]}"; do
    for bas in "${BASES[@]}"; do
        for m in "${METHODS[@]}"; do
            TARGETS+=("${m}_${mol,,}_${bas,,}")
        done
    done
done

if [ $# -gt 0 ]; then
    TARGETS=("$@")
fi

TABLE=()

for t in "${TARGETS[@]}"; do
    basis="${t##*_}"; rest="${t%_*}"; method="${rest%_*}"; system="${rest##*_}"

    echo "Method: $method, System: $system, Basis: $basis"

    "$ZINQ" "zinq_${t}.json" > "zinq_${t}.out"
    "$ORCA" "orca_${t}.inp"  > "orca_${t}.out"

    ez=$(awk '/(FINAL|MP0).*Eh/ {print $(NF-1)}' "zinq_${t}.out" | tail -n 1)

    tz=$(awk '/TOTAL EXECUTION TIME:/ {print $NF}' "zinq_${t}.out")

    eo=$(grep "FINAL SINGLE POINT ENERGY" "orca_${t}.out" | tail -n 1 | awk '{print $5}')

    to=$(awk '/TOTAL RUN TIME:/ { s = $4 * 86400 + $6 * 3600 + $8 * 60 + $10; if (s == 0) printf "%dms", $12; else printf "%.3fs", s + $12 / 1000 }' "orca_${t}.out")

    de=$(awk -v z="$ez" -v o="$eo" 'BEGIN {printf "%.8f", z - o}')

    printf "  ZINQ: %14.8f Eh (%s)\n" "$ez" "$tz"
    printf "  ORCA: %14.8f Eh (%s)\n" "$eo" "$to"

    printf "  DE:   %14.8f Eh\n\n" "$de"

    row=$(printf "%-12s %-8s %-10s %16.8f %-12s %16.8f %-10s %16.8f" "$method" "$system" "$basis" "$ez" "($tz)" "$eo" "($to)" "$de")

    TABLE+=("$row")
done

printf "%-12s %-8s %-10s %16s %-12s %16s %-10s %16s\n" "Method" "System" "Basis" "ZINQ (Eh)" "" "ORCA (Eh)" "" "DE (Eh)"

printf '%s\n' "-----------------------------------------------------------------------------------------------------------"

printf '%s\n' "${TABLE[@]}"
