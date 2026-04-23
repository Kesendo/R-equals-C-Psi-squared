#!/usr/bin/env bash
# EQ-024 Brecher scan: N=5, 7, 9 with 3 receivers and 4 J-profiles each.
# Runs C# brecher mode and collects RESULT lines.
# Usage: ./run_brecher_scan.sh [N_MAX]   (default N_MAX=7)
set -e

N_MAX=${1:-7}
OUT_FILE="../../simulations/results/eq024_refinement/brecher_scan_csharp.txt"
mkdir -p "$(dirname "$OUT_FILE")"

{
    echo "EQ-024 Brecher scan (C#)"
    echo "========================"
    echo "Date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "N_MAX=$N_MAX"
    echo
} > "$OUT_FILE"

run_case() {
    local n="$1"
    local j="$2"
    local init="$3"
    local label="$4"
    echo "  [N=$n] $label J=$j Initial=$init ..."
    # Measurements every 0.1 up to t=2, then 0.5 up to t=15
    dotnet run -c Release --no-build -- brecher "$n" "$j" "$init" --gamma 0.05 --tmax 15 --dt 0.05 \
      | tee -a "$OUT_FILE"
}

# N=5
if [ "$N_MAX" -ge 5 ]; then
    echo "=== N=5 ==="
    {
        echo "=== N=5 ==="
    } >> "$OUT_FILE"
    J_UNI="1,1,1,1"
    J_CUT="1,1,0.01,1"
    J_SW="5,0.2,5,0.2"

    for init in "plus" "xpattern:+-+-+" "bits:01010" "xpattern:+0+0+"; do
        for label_j in "uniform:$J_UNI" "cut-center:$J_CUT" "strong-weak:$J_SW"; do
            label="${label_j%%:*}"
            j="${label_j##*:}"
            run_case 5 "$j" "$init" "$label"
        done
    done
fi

# N=7
if [ "$N_MAX" -ge 7 ]; then
    echo "=== N=7 ==="
    {
        echo
        echo "=== N=7 ==="
    } >> "$OUT_FILE"
    J_UNI="1,1,1,1,1,1"
    J_CUT="1,1,1,0.01,1,1"  # cut center bond 3-4
    J_SW="5,0.2,5,0.2,5,0.2"

    # N=7 F71-symmetric analogs:
    #   |+-+-+-+> alternating (sites 0=6=+, 1=5=-, 2=4=+, 3=center=-)
    #   |0101010> alternating (site 3 = center = |1>)
    #   |+0+0+0+> alternating (site 3 = center = |0>)
    for init in "plus" "xpattern:+-+-+-+" "bits:0101010" "xpattern:+0+0+0+"; do
        for label_j in "uniform:$J_UNI" "cut-center:$J_CUT" "strong-weak:$J_SW"; do
            label="${label_j%%:*}"
            j="${label_j##*:}"
            run_case 7 "$j" "$init" "$label"
        done
    done
fi

# N=9
if [ "$N_MAX" -ge 9 ]; then
    echo "=== N=9 ==="
    {
        echo
        echo "=== N=9 ==="
    } >> "$OUT_FILE"
    J_UNI="1,1,1,1,1,1,1,1"
    J_CUT="1,1,1,1,0.01,1,1,1"  # cut center bond 4-5
    J_SW="5,0.2,5,0.2,5,0.2,5,0.2"

    for init in "plus" "xpattern:+-+-+-+-+" "bits:010101010" "xpattern:+0+0+0+0+"; do
        for label_j in "uniform:$J_UNI" "cut-center:$J_CUT" "strong-weak:$J_SW"; do
            label="${label_j%%:*}"
            j="${label_j##*:}"
            run_case 9 "$j" "$init" "$label"
        done
    done
fi

echo
echo "Done. Results in $OUT_FILE"
