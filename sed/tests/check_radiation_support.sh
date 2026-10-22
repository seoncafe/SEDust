#!/usr/bin/env bash
# Build the focused absorption/transition tests without touching library products.
set -euo pipefail
source_dir=$(cd "$(dirname "$0")/../src" && pwd)
test_dir=$(cd "$(dirname "$0")" && pwd)
test_build=$(mktemp -d "${TMPDIR:-/tmp}/sedust-support.XXXXXX")
test_compiler=${FC:-ifort}
case "$(basename "$test_compiler")" in
    gfortran*) test_flags='-O2 -fopenmp -fcheck=all -ffree-line-length-none' ;;
    *) test_flags='-O2 -qopenmp -recursive -fp-model precise -check bounds' ;;
esac
cd "$test_build"
for source in constants sed_mathlib radfield p_sub; do
    "$test_compiler" $test_flags \
        -c "$source_dir/$source.f90" -o "$source.o"
done
"$test_compiler" $test_flags \
    "$test_dir/test_radiation_support.f90" constants.o sed_mathlib.o radfield.o p_sub.o \
    -o test_radiation_support.x
./test_radiation_support.x
"$test_compiler" $test_flags \
    "$test_dir/test_absorption_balance.f90" constants.o sed_mathlib.o radfield.o p_sub.o \
    -o test_absorption_balance.x
./test_absorption_balance.x | tail -n 2
printf 'Test products: %s\n' "$test_build"
