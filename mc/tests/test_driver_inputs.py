#!/usr/bin/env python3
"""Exercise MC driver namelists without writing products into the source tree."""

from pathlib import Path
import subprocess
import sys
import tempfile


def run(executable, workdir, scratch, name, group, entries, expected, success=False):
    namelist = scratch / f"{name}.nml"
    namelist.write_text(f"&{group}\n{entries}\n/\n")
    result = subprocess.run(
        [str(executable), str(namelist)], cwd=workdir, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60,
        check=False,
    )
    assert (result.returncode == 0) == success, (name, result.stdout)
    assert expected in result.stdout, (name, result.stdout)
    assert "*****" not in result.stdout, (name, result.stdout)
    return result.stdout


def main():
    if len(sys.argv) != 3:
        raise SystemExit("Usage: test_driver_inputs.py MAIN_MC_X MAIN_PT_COMPARE_X")
    mc = Path(sys.argv[1]).resolve()
    compare = Path(sys.argv[2]).resolve()
    workdir = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory(prefix="sedust_mc_driver_") as temporary:
        scratch = Path(temporary)
        for name, entry, error in (
            ("radius", "a_um=0", "a_um must be finite and positive"),
            ("field", "U_isrf=-1", "U_isrf must be finite and nonnegative"),
            ("cutoff", "lam_c_um=0", "lam_c_um must be finite and positive"),
            ("events", "N_events=0", "N_events must be positive"),
        ):
            run(mc, workdir, scratch, name, "mc_input", entry, error)
        output = scratch / "single"
        run(mc, workdir, scratch, "single", "mc_input",
            f"T_init=20, N_events=2, t_max=0, out_prefix='{output}'",
            "emit/event_abs =", success=True)
        assert (scratch / "single_PT.dat").is_file()

        qtable = workdir.parent / "data/astrodust/q_astrodust_P0.20_Fe0.00_1.400.dat"
        assert qtable.is_file(), qtable
        output = scratch / "two_grains"
        run(compare, workdir, scratch, "two_grains", "pT_compare_input",
            f"qtable_path='{qtable}', a_um_list=0.01,0.03, "
            f"N_events=10, NT_init=40, out_prefix='{output}'",
            "n_grains=2", success=True)
        grains = sorted(scratch.glob("two_grains_??_a*A.dat"))
        assert len(grains) == 2, grains
        run(compare, workdir, scratch, "bad_radius", "pT_compare_input",
            "a_um_list=-0.01", "a_um_list must contain finite, nonnegative radii")
        run(compare, workdir, scratch, "bad_temperature", "pT_compare_input",
            "T_hi=0.5", "T_hi must be finite and greater than T_lo")
    print("MC driver input tests passed")


if __name__ == "__main__":
    main()
