"""This script calculates MTP error on each subset."""
# import numpy as np
# import matplotlib as mpl
import matplotlib.pyplot as plt


def get_between_words(word1, word2, line):
    return line.split(word1)[1].split(word2)[0]


def error_on_subset(error_each_files: list[str]):
    """error_each_files: output of mlp3 check_errors --log
    ex) run
    mlp3ser check_errors --log=error_relax_each.txt pot.almtp train_subset.cfg  | tee error_relax.txt
    ... --log=error_water_each.txt ...

    then
    error_each_files = ['error_rela_each.txt', 'error_water_each.txt']
    """
    for error_each in error_each_files:
        error_per_atom = []
        rms_force =[ ]
        with open(error_each, 'r') as fin:
            for line in fin:
                error_per_atom.append(float(get_between_words('Diff(epa):', 'RMS(force)', line)))
                rms_force.append(float(get_between_words('RMS(force):', 'MAX(force)', line)))

        # np.savetxt(f'{error_each}_Diff_epa', error_per_atom, fmt='%.6g')

        fig, axs = plt.subplots(nrows=2, ncols=1)
        ax=axs[0]
        ax.plot(error_per_atom)
        ax.set_ylabel('Diff(epa)')
        ax=axs[1]
        ax.plot(rms_force)
        ax.set_ylabel('RMS(force)')
                                      
        fig.savefig(f'fig_{error_each}.pdf')
        plt.show()

