#!/bin/sh
# Run from the isolated comparison bench directory, after the GPU smoke check.
set -eu
: "${CUDA_VISIBLE_DEVICES:?Select the intended GPU explicitly}"
comparison_python=${FRC_COMPARISON_PYTHON:-../.venv/bin/python}
comparison_runs=${FRC_COMPARISON_RUN_ROOT:-runs}
mkdir -p "$comparison_runs"
run_case() {
    comparison_name=$1
    shift
    if test -f "$comparison_runs/$comparison_name/measurement.json"; then
        # Final archive audit checks completion and exact current source hashes.
        return
    fi
    "$comparison_python" -m examples.warpx_comparison "$@" \
        --output "$comparison_runs/$comparison_name" \
        > "$comparison_runs/$comparison_name.log" 2>&1
    printf 'Completed %s\n' "$comparison_name"
}
for comparison_cells in 32 64 128; do
    for comparison_helicity in -1 1; do
        run_case "hall_h${comparison_helicity}_n${comparison_cells}" \
            --kind hall --helicity "$comparison_helicity" --cells "$comparison_cells"
    done
    run_case "vacuum_n${comparison_cells}" --kind vacuum --cells "$comparison_cells"
done
for comparison_helicity in -1 1; do
    run_case "hall_h${comparison_helicity}_n128_half_dt" \
        --kind hall --helicity "$comparison_helicity" --cells 128 --dt-omega .0025
    run_case "hall_h${comparison_helicity}_n128_quarter_dt" \
        --kind hall --helicity "$comparison_helicity" --cells 128 --dt-omega .00125
done
run_case hall_h1_n128_ppc64 --kind hall --cells 128 --ppc 64
run_case hall_h1_n128_small_amplitude --kind hall --cells 128 --amplitude 1e-8
run_case hall_h1_n128_floor_control --kind hall --cells 128 --floor-fraction .01
run_case vacuum_n128_half_dt --kind vacuum --cells 128 --cfl .25
