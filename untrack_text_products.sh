#!/bin/bash
# Stop shipping the optics text products that make_qtable.x and calc_kext.x
# regenerate in minutes.  They stay on disk; only git forgets them.
#
#   bash untrack_text_products.sh        show what would happen
#   bash untrack_text_products.sh --go   do it (still leaves the commit to you)
#
# The T-matrix sweep goes too: it is hours of solver time rather than minutes,
# but the tree regenerates it (cd tmatrix && make && ./run_tmatrix.x, then
# make lyman_cut) and the same numbers are in the model's .h5, which ships.
# Its two files under tmatrix/output/ are ALREADY gone from the working tree --
# a model reads the installed copy under data/astrodust/, so those were a second
# copy of the same bytes.  They are still in the index, which is why they are
# listed: this is what takes them out of it.
#
# NOT removed: data/dielectric/q_D16graphite.dat, a material table and not a
# product of this tree.
set -e
cd "$(dirname "$0")"

FILES=(
  data/astrodust/kext_astrodust_MW.dat
  data/astrodust/kext_astrodust_MW_euv.dat
  data/astrodust/q_astrodust_P0.20_Fe0.00_1.400.dat
  data/astrodust/q_astrodust_P0.20_Fe0.00_1.400_euv.dat
  data/astrodust/q_astrodust_pah_ion.dat
  data/astrodust/q_astrodust_pah_neu.dat
  data/dl07/kext_dl07_MW.dat
  data/dl07/kext_dl07_MW_euv.dat
  data/dl07/q_dl07_gra.dat
  data/dl07/q_dl07_gra_euv.dat
  data/dl07/q_dl07_pah_ion.dat
  data/dl07/q_dl07_pah_ion_euv.dat
  data/dl07/q_dl07_pah_neu.dat
  data/dl07/q_dl07_pah_neu_euv.dat
  data/dl07/q_dl07_sil.dat
  data/dl07/q_dl07_sil_euv.dat
  data/zubko/kext_zubko_BARE_GR_S.dat
  data/zubko/kext_zubko_BARE_GR_S_euv.dat
  data/zubko/q_zubko_gra.dat
  data/zubko/q_zubko_gra_euv.dat
  data/zubko/q_zubko_pah.dat
  data/zubko/q_zubko_pah_euv.dat
  data/zubko/q_zubko_sil.dat
  data/zubko/q_zubko_sil_euv.dat
  tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat
  tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400_euv.dat
)

if [ "$1" != "--go" ]; then
  echo "Would untrack ${#FILES[@]} files (they stay on disk):"
  printf '  %s\n' "${FILES[@]}"
  du -ch "${FILES[@]}" 2>/dev/null | tail -1
  echo
  echo "Run with --go to stage the removal, then commit and push yourself."
  exit 0
fi

git rm --cached -q "${FILES[@]}"
#git add .gitignore
#echo "Staged.  Now, when you are ready:"
#echo "    git commit -m 'stop shipping the regenerable optics text products'"
#echo "    git push"
#echo
#echo "Note: this drops them from the tip only.  The blobs already pushed stay"
#echo "in history, so the clone size does not shrink without a history rewrite."
