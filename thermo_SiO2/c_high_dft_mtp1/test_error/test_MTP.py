#!/usr/bin/env python3
"""This script test MTP in ./original_train_folder/ to configurations in allcfgs.cfg.
"""

import os.path as path
from a_parameters import mtp_levels, train_folder
from ...util.util import shell, read_output, if_exists_delete_file

VERBOSE = False
# VERBOSE = True  # debug

RMSE_output = (
    "# Level of MTP   RMSE of energy per atom (eV/atom)    "
    "RMSE of force (eV/Ang) \n"
)

for mtp_level in mtp_levels:
    print("\n MTP level:", mtp_level)
    last_mtp_folder = read_output(
        f"ls {train_folder}/" f"mtp_level_{mtp_level}/last_mtp_in_mlp_*"
    )
    last_mtp_folder = (
        last_mtp_folder.split("/")[-1].replace("last_mtp_in_", "").replace(" ", "")
    )
    if VERBOSE:
        print("last_mtp_folder =", last_mtp_folder)
    output_file = f"output_error_level_{mtp_level}.out"

    # Justus 2
    # command_str = 'module load compiler/intel && ' \
    #     'module load numlib/mkl && ' \
    #     f'mlp calc-errors ./original_train_folder/mtp_level_{mtp_level}' \
    #     f'/{last_mtp_folder}/Trained.mtp_  allcfgs.cfg '

    # Fritz server
    # command_str = 'module load    intel  mkl  intelmpi && ' \
    #     f'mlp.2020-07.serial calc-errors ./original_train_folder/mtp_level_{mtp_level}' \
    #     f'/{last_mtp_folder}/Trained.mtp_  allcfgs.cfg '

    source_file_dir = path.dirname(path.abspath(__file__))
    command_str = (
        f"{source_file_dir}/run_mlp_calc-errors.sh {train_folder}/mtp_level_{mtp_level}"
        f"/{last_mtp_folder}/Trained.mtp_  allcfgs.cfg "
    )
    if VERBOSE:
        command_str = command_str + f"| tee {output_file}"
    else:
        command_str = command_str + f"> {output_file}"
    if VERBOSE:
        print("command_str =", command_str)
    shell(command_str)

    # read RMS error
    with open(output_file, "r") as file:
        lines = file.readlines()
        lines = [line.rstrip() for line in lines]
    all_lines = " ".join(lines)

    # example output
    # Energy per atom:
    #     Errors checked for 42 configurations
    #     Maximal absolute difference = 0.0202603
    #     Average absolute difference = 0.00668585
    #     RMS     absolute difference = 0.00815697
    lines = all_lines.split(sep="Energy per atom:", maxsplit=1)[1]
    # [1]: after the separator
    lines = lines.split(sep="RMS     absolute difference =", maxsplit=1)[1]
    RMS_energy = lines.split()[0]

    lines = all_lines.split(sep="Forces:", maxsplit=1)[1]
    # [1]: after the separator
    lines = lines.split(sep="RMS     absolute difference =", maxsplit=1)[1]
    RMS_force = lines.split()[0]

    # RMS_energy = float(RMS_energy) * 1000
    print("RMS_energy =", RMS_energy, "eV/atom")
    print("RMS_force =", RMS_force, "eV/Ang")

    # RMSE_output += f'{mtp_level}  {RMS_energy:.6g}  {RMS_force}\n'
    RMSE_output += f"{mtp_level}  {RMS_energy}  {RMS_force}\n"

    # sys.exit()  # debug

print()
if_exists_delete_file("RMSE.txt")
with open("RMSE.txt", "w") as file:
    file.write(RMSE_output)

print("Output file RMSE.txt is written.")
