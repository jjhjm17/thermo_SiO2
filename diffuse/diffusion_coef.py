"""This script obtains the average diffusion coefficient D from lammps stdout."""
import sys
import numpy as np
from pathlib import Path
from thermo_SiO2.util.util import read_output

# def diffusion_coef(calc_dirs='calc_*'):
def diffusion_coef_lmp_out(calc_dirs='calc_*'):
    txt = ''
    calc_dirs_glob = sorted(Path('./').glob(calc_dirs + '/'))
    for calc_dir in calc_dirs_glob:
        coef_s = []
        step_s = []
        txt += f'calc_dir: {calc_dir}\n'
        for folder in sorted(Path(calc_dir).glob('./*/')): # ex) 000, 001, ...
            # folders.append(folder)
            last_line_num = 75
            line = read_output(f'tail -n{last_line_num} {folder}/lammps*.out | head -n1').split()
            T_last = float(line[1])
            if T_last < 200 or T_last > 400:
                print('Error: T < 200 or T > 400, please check.')
                breakpoint()
                sys.exit()
            step_s.append(int(line[0]))
            coef_s.append(float(line[-1]))

        txt += f'avg: {np.average(coef_s):5f}, '
        txt += f'std: {np.std(coef_s, ddof=1):5f} Ang^2/ps, '
        txt += f'avg_step: {np.average(step_s)}\n\n'
        #  breakpoint()
    with open('out_diff_coef.txt', 'w') as fout:
        print(txt)
        fout.write(txt)





