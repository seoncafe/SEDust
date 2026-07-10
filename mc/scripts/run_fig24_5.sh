#!/bin/bash
# Generate namelists and run MC for the PIIM Figure 24.5 panels:
# (a, U) = {10, 20, 50, 100, 200 A} x {1, 100}, t_max = 1e5 s each.
set -e
cd "$(dirname "$0")/.."

mkdir -p output/fig24_5
SIZES=(0.001 0.002 0.005 0.010 0.020)        # um
SIZE_LBL=(10A  20A  50A  100A 200A)
U_VALS=(1.0 100.0)
U_LBL=(U1 U100)

for i in 0 1 2 3 4; do
  a=${SIZES[$i]}
  alabel=${SIZE_LBL[$i]}
  for j in 0 1; do
    U=${U_VALS[$j]}
    Ulabel=${U_LBL[$j]}
    tag="${alabel}_${Ulabel}"
    nml="output/fig24_5/${tag}.nml"
    cat > "$nml" <<EOF
&mc_input
  a_um       = $a
  U_isrf     = $U
  comp       = 'gra_dl01'
  N_events   = 200000
  t_max      = 1.0e5
  lam_c_um   = 1000.0
  seed       = $((42 + 10*i + j))
  out_prefix = 'output/fig24_5/${tag}'
  record_trajectory = .true.
/
EOF
    echo "=== $tag (a=$a um, U=$U) ==="
    ./main_mc.x "$nml" 2>&1 | grep -E '(rate_event|emit/abs|t_total|n_rec)'
  done
done

echo
echo "All runs done. Outputs in output/fig24_5/"
