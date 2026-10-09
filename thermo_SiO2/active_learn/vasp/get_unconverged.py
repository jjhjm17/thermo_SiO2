#!/usr/bin/env python3
"""This script obtains unconverged folder names after vasp calculations.

We read parameters from parameters.py file like ./makeInputForVasp.py .
We check if vasp*out file finished and if there is no mention that EDIFF was not reached, in both pre-convergence and the next step, if there is preconvergence.
Otherwise, print out the folders in unconverged.txt

Keep the above comment as it is hand-written. Additional comments can be added below.
"""

import argparse
from collections import Counter
import os
from pathlib import Path
import re
import tempfile

import parameters as param


SCF_ITERATION = re.compile(r"\bIteration\s+\d+\(\s*\d+\)")
CONVERGED = "aborting loop because EDIFF is reached"
FORCES = "POSITION                                       TOTAL-FORCE"
ENERGY = "free  energy   TOTEN"
FINISHED = "General timing and accounting informations for this job:"
PRECONVERGE_START = "We pre-converge WAVECAR at folder preconverge."
MAIN_START = "We calculate by reading pre-converged WAVECAR."
SCF_WARNING = "The electronic self-consistency was not achieved"
FINAL_F = re.compile(r"^\s*\d+\s+F=")
ELAPSED = "Elapsed time (sec):"
JOB_FINISHED = "Job Wall-clock time:"


def check_outcar(path):
    """Return a short reason if the final single-point result is unusable."""
    if not path.is_file():
        return "missing_OUTCAR"
    if path.stat().st_size == 0:
        return "empty_OUTCAR"

    last_iteration = last_converged = last_forces = last_energy = last_finish = 0
    try:
        with path.open(errors="replace") as stream:
            for number, line in enumerate(stream, 1):
                if SCF_ITERATION.search(line):
                    last_iteration = number
                if CONVERGED in line:
                    last_converged = number
                if FORCES in line:
                    last_forces = number
                if ENERGY in line:
                    last_energy = number
                if FINISHED in line:
                    last_finish = number
    except OSError as error:
        return f"unreadable_OUTCAR:{error.strerror or type(error).__name__}"

    if not last_iteration:
        return "no_electronic_steps"
    if last_converged < last_iteration:
        return "electronic_not_converged"
    if last_forces < last_converged or last_energy < last_converged:
        return "missing_final_energy_or_forces"
    if last_finish < max(last_forces, last_energy):
        return "incomplete_OUTCAR"
    return None


def check_stdout(folder, preconverge):
    """Check both VASP runs in the submitted job, when preconvergence is used.

    The analysis script b_analyze_mag_vs_no_mag.py associates DAV and F lines
    within each run. Here the explicit VASP SCF warning distinguishes a run
    that reached NELM but still printed an F line.
    """
    files = sorted(folder.glob("vasp*.out"))
    if not files:
        return "missing_vasp_stdout"
    if len(files) != 1:
        return "ambiguous_vasp_stdout"

    stages = {"preconvergence": {"f": 0, "elapsed": 0, "warning": False},
              "main": {"f": 0, "elapsed": 0, "warning": False}}
    stage = "preconvergence" if preconverge else "main"
    saw_preconvergence = saw_main = finished = False
    try:
        with files[0].open(errors="replace") as stream:
            for line in stream:
                if PRECONVERGE_START in line:
                    saw_preconvergence = True
                if MAIN_START in line:
                    stage = "main"
                    saw_main = True
                if FINAL_F.search(line):
                    stages[stage]["f"] += 1
                if ELAPSED in line:
                    stages[stage]["elapsed"] += 1
                if SCF_WARNING in line or ("EDIFF" in line and "not reached" in line):
                    stages[stage]["warning"] = True
                if JOB_FINISHED in line:
                    finished = True
    except OSError as error:
        return f"unreadable_vasp_stdout:{error.strerror or type(error).__name__}"

    if preconverge and (not saw_preconvergence or not saw_main):
        return "incomplete_preconvergence_stdout"
    for name in ("preconvergence", "main") if preconverge else ("main",):
        result = stages[name]
        if result["warning"]:
            return f"{name}_not_converged"
        if not result["f"] or not result["elapsed"]:
            return f"incomplete_{name}_stdout"
    if not finished:
        return "incomplete_vasp_stdout"
    return None


def write_atomic(path, lines):
    """Replace a report only after its complete contents have been written."""
    path = Path(path)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", dir=path.parent, prefix=f".{path.name}.",
        delete=False,
    ) as stream:
        temporary = Path(stream.name)
        try:
            stream.writelines(lines)
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    os.replace(temporary, path)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Find unfinished or unconverged VASP calculations")
    parser.add_argument("--output", default="unconverged.txt")
    parser.add_argument("--report", default="unconverged_report.tsv")
    args = parser.parse_args(argv)

    root = Path(param.calc_folder).resolve()
    if not root.is_dir():
        parser.error(f"calculation folder does not exist: {root}")
    try:
        first = int(param.start_config_number)
        last = int(param.end_config_number)
    except AttributeError:
        parser.error("parameters.py must define start_config_number and end_config_number")
    if first < 0 or last < first:
        parser.error("invalid configuration number range")

    selected = []
    reasons = Counter()
    for index in range(first, last + 1):
        folder = root / f"{index:04d}"
        source = folder
        if not (folder / "OUTCAR").exists() and (folder / "unconverged").is_dir():
            source = folder / "unconverged"
        reason = check_outcar(source / "OUTCAR")
        if reason is None:
            reason = check_stdout(source, getattr(param, "preconverge", False))
        if reason is not None:
            name = os.path.relpath(folder, Path.cwd())
            selected.append((name, reason))
            reasons[reason] += 1

    write_atomic(args.output, (f"{name}\n" for name, _ in selected))
    write_atomic(args.report, (f"path\treason\n", *(f"{name}\t{reason}\n" for name, reason in selected)))
    print(f"Checked {last - first + 1} calculations in {root}")
    print(f"Selected {len(selected)} for rerun: {args.output}")
    print(f"Reasons: {args.report}")
    for reason, count in sorted(reasons.items()):
        print(f"  {reason}: {count}")


if __name__ == "__main__":
    main()
