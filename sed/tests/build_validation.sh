#!/bin/sh
# Build the validation driver against an already built SEDust library.
set -eu
if [ "$#" -ne 3 ]; then
   echo 'Usage: build_validation.sh SED_DIRECTORY exact|corrected EXECUTABLE' >&2
   exit 2
fi
sed_directory=$1
mode=$2
executable=$3
test_compiler=${FC:-ifort}
case "$(basename "$test_compiler")" in
   gfortran*)
      test_flags='-O3 -ffree-line-length-none -cpp -fopenmp -frecursive'
      hdf5_prefix=${HDF5_PREFIX:-/opt/hdf5_gfortran}
      ;;
   *)
      test_flags='-O3 -qopenmp -recursive -fp-model precise -fpp -w'
      hdf5_prefix=${HDF5_PREFIX:-/opt/hdf5_intel}
      ;;
esac
case "$mode" in
   exact) definitions= ;;
   corrected) definitions='-DFIELD_SUPPORT -DBIN_POWER' ;;
   *) echo 'Mode must be exact or corrected' >&2; exit 2 ;;
esac
script_directory=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
"$test_compiler" $test_flags \
   $definitions -I"$sed_directory" -I"$hdf5_prefix/include" \
   "$script_directory/validation_driver.f90" "$sed_directory/libsedust.a" \
   "$hdf5_prefix/lib/libhdf5_fortran.a" "$hdf5_prefix/lib/libhdf5.a" \
   -lsz -lz -ldl -o "$executable"
