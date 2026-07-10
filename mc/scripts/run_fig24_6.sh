#!/bin/bash
# Run MC for P(T) distributions (PIIM Fig 24.6).
# Sizes: 10, 20, 50, 100, 200 A.  ISRF: U=1 (Mathis 1983).
# N_events = 20000 per grain for histogram statistics.
set -e
cd "$(dirname "$0")/.."

mkdir -p output/fig24_6
SIZES=(0.001 0.002 0.005 0.010 0.020)
SIZE_LBL=(10A  20A  50A  100A 200A)
NEVT=(500000 500000 500000 500000 500000)

for i in 0 1 2 3 4; do
  a=${SIZES[$i]}
  alabel=${SIZE_LBL[$i]}
  N=${NEVT[$i]}
  tag="${alabel}_U1"
  nml="output/fig24_6/${tag}.nml"
  cat > "$nml" <<EOF
&mc_input
  a_um       = $a
  U_isrf     = 1.0
  comp       = 'gra_dl01'
  N_events   = $N
  t_max      = 0.0
  lam_c_um   = 1000.0
  seed       = $((101 + i))
  out_prefix = 'output/fig24_6/${tag}'
/
EOF
  echo "=== $tag (a=$a um, N=$N) ==="
  ./main_mc.x "$nml" 2>&1 | grep -E '(rate_event|emit/abs|t_total|n_rec)'
done

echo
echo "All P(T) runs done."
