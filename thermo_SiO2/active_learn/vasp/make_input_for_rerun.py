"""
This script prepares input files for rerunning vasp.
After first vasp calculations, some vasp calculations are not converged.
Some stops in the middle, and some has text 'EDIFF is not reached', which is checked in
get_unconverged.py

We read parameters.py file.
We get unconverged folders from (for example) argument unconverged_dir = unconverged.txt
Then move the content of the folder to subfolder './unconverged'.
Do similar to makeInputForVasp.py, but first do not make new folder.
Copy input files, change blanks, make joblist, and prepare for the next submission.

Keep the above comment as it is hand-written. Additional comments can be added below.
"""

import argparse
import math
import os
from pathlib import Path
import re
import shutil
import tempfile
from textwrap import wrap

from ase.io import write
import numpy as np

import parameters as param
from ...util.SiO2_parameter import POTCAR_setup
from .get_unconverged import write_atomic
from .makeInputForVasp import get_atoms_and_forces, unique_ordered_list


TEMPLATE = Path("template_e_algo_all")
ZVAL = re.compile(r"\bZVAL\s*=\s*([+-]?\d+(?:\.\d+)?)")


def read_targets(list_file, calc_root):
    """Accept only numbered calculation directories in the selected root."""
    targets = []
    seen = set()
    with list_file.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            name = line.strip()
            if not name or name.startswith("#"):
                continue
            folder = Path(name).resolve()
            try:
                index = int(folder.name)
            except ValueError as error:
                raise ValueError(f"{list_file}:{line_number}: expected a numbered folder") from error
            if folder.parent != calc_root or folder.name != f"{index:04d}":
                raise ValueError(f"{list_file}:{line_number}: folder is outside {calc_root}: {name}")
            if not param.start_config_number <= index <= param.end_config_number:
                raise ValueError(f"{list_file}:{line_number}: index is outside the selected range: {name}")
            if folder in seen:
                raise ValueError(f"{list_file}:{line_number}: duplicate folder: {name}")
            if not folder.is_dir():
                raise ValueError(f"{list_file}:{line_number}: folder does not exist: {name}")
            archive = folder / "unconverged"
            prepared = archive.is_dir() and not archive.is_symlink() and not (folder / "OUTCAR").exists()
            if archive.exists() or archive.is_symlink():
                required = ("POSCAR", "POTCAR", "INCAR")
                if getattr(param, "preconverge", False):
                    required += ("INCAR.preconverge",)
                if not prepared or not all((folder / item).is_file() for item in required):
                    raise ValueError(f"{folder}: existing unconverged/ archive needs review")
            seen.add(folder)
            targets.append((index, folder, prepared))
    return targets


def potcars_for(atoms, potential_root):
    """Use the same POTCAR species order and NBANDS formula as makeInputForVasp."""
    symbols = atoms.get_chemical_symbols()
    unique = unique_ordered_list(symbols)
    sources = []
    electron_count = 0.0
    for symbol in unique:
        try:
            source = potential_root / POTCAR_setup[symbol] / "POTCAR"
        except KeyError as error:
            raise ValueError(f"no POTCAR mapping for {symbol}") from error
        data = source.read_bytes()
        match = ZVAL.search(data.decode("latin-1"))
        if match is None:
            raise ValueError(f"ZVAL missing from {source}")
        electron_count += symbols.count(symbol) * float(match.group(1))
        sources.append(data)
    nbands = math.ceil(electron_count / 2 + len(atoms) / 4)
    return b"".join(sources), nbands


def magmom_for(index, atoms):
    value = getattr(param, "set_magmom_value", None)
    if value is None:
        raise ValueError("set_magmom_value is required by the rerun INCAR template")
    if not getattr(param, "magmom_random_sign", False):
        return f"{len(atoms)}*{value}"
    seed = param.magmom_seed * index
    signs = np.random.default_rng(seed).choice((-1, 1), size=len(atoms))
    values = " ".join(str(value * int(sign)) for sign in signs)
    separator = " " + chr(92) + "\n         "
    return separator.join(wrap(values, width=80, break_long_words=False, break_on_hyphens=False))


def render_incar(source, destination, nbands, magmom):
    text = source.read_text(encoding="utf-8")
    substitutions = {
        "xxx__NBANDS__xxx": f"{nbands}  # NELECT / 2 + NIONS / 4",
        "xxx__MAGMOM__xxx": magmom,
    }
    for placeholder, replacement in substitutions.items():
        text = text.replace(placeholder, replacement)
    if "xxx__" in text:
        raise ValueError(f"unfilled placeholder in {source}")
    destination.write_text(text, encoding="utf-8")


def prepare_inputs(stage, destination, atoms, index, potential_root, template):
    """Build complete new inputs away from the existing calculation."""
    stage.mkdir()
    write(stage / "POSCAR", atoms, format="vasp", vasp5=True)
    potcar, nbands = potcars_for(atoms, potential_root)
    (stage / "POTCAR").write_bytes(potcar)
    magmom = magmom_for(index, atoms)
    render_incar(template / "INCAR", stage / "INCAR", nbands, magmom)

    if getattr(param, "preconverge", False):
        change = template / "INCAR.preconverge.change"
        full = template / "INCAR.preconverge"
        if change.is_file():
            render_incar(change, stage / change.name, nbands, magmom)
        elif full.is_file():
            render_incar(full, stage / full.name, nbands, magmom)
        else:
            raise ValueError(f"preconvergence INCAR missing from {template}")

    for name in ("KPOINTS",):
        source = template / name
        if source.is_file():
            shutil.copy2(source, stage / name)
    kernel = template / "vdw_kernel.bindat"
    if kernel.is_file():
        (stage / kernel.name).symlink_to(os.path.relpath(kernel, destination))


def archive_and_install(folder, stage):
    """Move originals into unconverged/ and roll back a failed folder move."""
    archive = folder / "unconverged"
    originals = list(folder.iterdir())
    moved = []
    installed = []
    external_links = []
    archive.mkdir()
    try:
        for source in originals:
            original_link = os.readlink(source) if source.is_symlink() else None
            resolved = source.resolve() if original_link is not None else None
            destination = archive / source.name
            source.rename(destination)
            moved.append((source, destination, original_link))
            if resolved is not None and not resolved.is_relative_to(folder):
                destination.unlink()
                destination.symlink_to(os.path.relpath(resolved, archive))
                external_links.append((destination, original_link))
        for source in list(stage.iterdir()):
            destination = folder / source.name
            source.rename(destination)
            installed.append((source, destination))
    except BaseException:
        for source, destination in reversed(installed):
            destination.rename(source)
        for destination, original_link in external_links:
            destination.unlink()
            destination.symlink_to(original_link)
        for source, destination, _ in reversed(moved):
            destination.rename(source)
        archive.rmdir()
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description="Prepare VASP rerun inputs from unconverged.txt")
    parser.add_argument("--list", default="unconverged.txt", type=Path)
    parser.add_argument("--job-list", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    root = Path(param.calc_folder).resolve()
    if not root.is_dir():
        parser.error(f"calculation folder does not exist: {root}")
    if not args.list.is_file():
        parser.error(f"unconverged list does not exist: {args.list}")
    job_list = args.job_list or Path(f"jobList_rerun_{root.name}")
    template = TEMPLATE.resolve()
    if not (template / "INCAR").is_file():
        parser.error(f"INCAR template does not exist: {template / 'INCAR'}")
    potential_root = Path(os.environ.get("pbepot", ""))
    if not os.environ.get("pbepot") or not potential_root.is_dir():
        parser.error("set $pbepot to the PBE POTCAR elements directory")

    targets = read_targets(args.list, root)
    if not targets:
        print(f"No calculations listed in {args.list}; no inputs changed")
        return
    already_prepared = {folder for _, folder, prepared in targets if prepared}
    if job_list.exists():
        old_jobs = [Path(line).resolve() for line in job_list.read_text().splitlines() if line.strip()]
        if len(old_jobs) != len(set(old_jobs)) or not set(old_jobs) <= already_prepared:
            parser.error(f"existing job list does not match prepared reruns: {job_list}")
    new_targets = [(index, folder) for index, folder, prepared in targets if not prepared]
    print(f"Selected {len(targets)} calculations in {root}; "
          f"{len(already_prepared)} already prepared, {len(new_targets)} new")
    if new_targets:
        atoms_and_forces = get_atoms_and_forces(param)
        if max(index for index, _ in new_targets) >= len(atoms_and_forces):
            raise ValueError("unconverged list refers to an index beyond the configuration file")

    if args.dry_run:
        for folder in sorted(already_prepared):
            print(f"Already prepared {folder}")
        for index, folder in new_targets:
            atoms = atoms_and_forces[index]["atoms"]
            potcars_for(atoms, potential_root)
            magmom_for(index, atoms)
            print(f"Would prepare {folder}")
        print("Dry run complete; no inputs changed")
        return

    completed = set(already_prepared)
    try:
        if new_targets:
            with tempfile.TemporaryDirectory(dir=root, prefix=".rerun_stage_") as staging:
                staging = Path(staging)
                for index, folder in new_targets:
                    prepare_inputs(
                        staging / folder.name, folder, atoms_and_forces[index]["atoms"],
                        index, potential_root, template,
                    )
                for _, folder in new_targets:
                    archive_and_install(folder, staging / folder.name)
                    completed.add(folder)
                    print(f"Prepared {folder}")
    finally:
        if completed:
            write_atomic(job_list, (f"{folder}\n" for _, folder, _ in targets if folder in completed))
            print(f"Wrote {len(completed)} job directories to {job_list}")


if __name__ == "__main__":
    main()
