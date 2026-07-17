import glob
import numpy as np
from scipy import stats
from ase.io import read, write
from thermo_SiO2.util.util import read_output
from thermo_SiO2.diffuse.diffuse_coef_vasp import get_diffusion_coeff_traj
import a_parameters as param

def get_outcar_length(outcar_path):
    """Get the number of ionic steps in an OUTCAR file."""
    #  len = read_output(f"awk ' /Iteration/ {iter=$3}  /EDIFF is reached/ {last_iter=iter} END {print last_iter}'  {outcar_path}")
    # len = int(len[:-1])  # ex) 2185( to 2185
    with open(outcar_path, 'r') as f:
        lines = f.readlines()
    count = sum(1 for line in lines if 'POSITION' in line)
    return count

def read_concat_outcars(base_dir, folder, start, step, cont_base_dirs=None):
    """Read and concatenate OUTCAR trajectories for a single folder (e.g., 0000).
      base_dir: a folder containing 1st outcars"""
    all_traj = []

    # --- First OUTCAR (main run) ---
    main_outcar = f"{base_dir}/{folder}/OUTCAR"
    print(f"Reading main OUTCAR: {main_outcar}")
    traj_main = read(main_outcar, index=f"{start}::" + str(step))
    # traj_main = read(main_outcar, index=f"{start}::")
    print(f"  frames read 1st: {len(traj_main)}")
    all_traj.extend(traj_main)
    len_last = get_outcar_length(main_outcar)
    start_last = param.start
    start_next = int(start_last + np.ceil((len_last - start_last) / step) * step - len_last)
    if (start_next < 0 or start_next >= step):
        assert ValueError("Error in calculating start_next for continuation runs.")
    # Keep a constant interval 'step' for continuation runs with varying lengths.

    # --- Continued runs ---
    if cont_base_dirs is None:
        return all_traj
    for cont_dir in cont_base_dirs:
        cont_outcar = f"{cont_dir}/{folder}/OUTCAR"
        try:
            traj_cont = read(cont_outcar, index=f"{start_next}::" + str(step))
            # traj_cont = read(cont_outcar, index=':')
            all_traj.extend(traj_cont)
            print(f"  Added continuation: {cont_outcar} ({len(traj_cont)} frames)")
            len_last = get_outcar_length(cont_outcar)
            start_last = start_next
            start_next = int(start_last + np.ceil((len_last - start_last) / step) * step - len_last)
            if (start_next < 0 or start_next >= step):
                assert ValueError("Error in calculating start_next for continuation runs.")
        except Exception as e:
            print(f"  Skipping {cont_outcar}: {e}")

    # all_traj = all_traj[::step]  # Resample to maintain consistent step size
    print(f" → Total frames combined: {len(all_traj)}")
    return all_traj

def main():
    # base_dir = "../b.low_AIMD/calc"
    base_dir = param.base_dir
    # cont_base_dirs = sorted(glob.glob("../b.low_AIMD/calc_cont_*"))
    if hasattr(param, 'cont_base_dirs'):
        cont_base_dirs = param.cont_base_dirs
    else:
        cont_base_dirs = None
    calc_dirs = sorted(glob.glob(f"{base_dir}/*/"))
    if hasattr(param, 'n_max_dir'):
        calc_dirs = calc_dirs[:param.n_max_dir]
    run_idxs = [d.split("/")[-2] for d in calc_dirs]

    results = []

    for run_idx in run_idxs:
        print(f"\nProcessing folder {run_idx} ...")
        traj = read_concat_outcars(base_dir, run_idx, param.start, param.step, cont_base_dirs)
        if hasattr(param, 'write_traj') and param.write_traj:
            # write(f'traj_{folder}.traj', traj)
            write(f'traj_{run_idx}.xyz', traj)
            print(f"Trajectory written to traj_{run_idx}.xyz")
            
        if hasattr(param, 'show_msd_plot'):
            show_msd_plot = param.show_msd_plot
        else:
            show_msd_plot = True
        slope, std = get_diffusion_coeff_traj(traj, param.step, show_msd_plot=show_msd_plot, run_idx=run_idx)
        results.append((run_idx, slope, std))

    # --- Collect average & std across all runs ---
    slopes = np.array([r[1] for r in results])
    stds   = np.array([r[2] for r in results])
    total_std  = np.std(slopes)
    n = len(slopes)
    mean_slope = np.mean(slopes)
    t_value = stats.t.ppf(0.975, df=n-1)  # 95% CI, two-tailed
    sem = total_std / np.sqrt(n)          # standard error of the mean
    ci_half = t_value * sem
    print(f'{t_value=}')
    print(f'{sem=}')
    print(f'{n=}')

    print("\n====== Summary ======")
    for run_idx, slope, std in results:
        print(f"{run_idx}: D = {slope:.4e} ± {std:.4e} (std) Å²/ps")
    print(f"\nAverage D = {mean_slope:.4e} ± {total_std:.4e} (std)  Å²/ps")
    print(f"\n  ± {ci_half:.4e} (95% CI of the mean)  Å²/ps")

    # Optional: save to file
    with open("diffusion_summary.txt", "w") as f:
        f.write("Folder\tD(Å²/ps)\tstd(Å²/ps)\n")
        for run_idx, slope, std in results:
            f.write(f"{run_idx}\t{slope:.6e}\t{std:.6e}\n")
        f.write(f"\nAverage\t{mean_slope:.6e}\t{total_std:.6e}\n")
        f.write(f"\n +-\t{ci_half:.6e} (95% CI of the mean)\n")

if __name__ == "__main__":
    main()

