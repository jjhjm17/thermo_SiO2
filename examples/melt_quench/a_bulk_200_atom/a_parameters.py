from glob import glob


# calc_folder = 'calc'
# calc_folder = 'calc_Al_0'
# template_folder = './template/a_Al_0/c_melt/'
# # num_seeds = 2
# num_seeds = 30
# num_mtps = 1

# calc_folder = 'calc_Al_0.2'
# template_folder = './template/b_Al_0.2/c_melt/'
# # num_seeds = 2
# num_seeds = 30
# num_mtps = 1

calc_folder = 'calc_Al_0.05'
template_folder = './template/c_Al_0.05/c_melt/'
# num_seeds = 2
num_seeds = 30
num_mtps = 1


server = 'fritz'
# server = 'justus'

use_initial_config = False
initial_config = None

symbols = 'Si O H Al'
# symbols = 'Si O Al'
different_seeds = True
seed_of_seeds = 12345

# new seed_of_seeds when the iteration proceeds!

# almtp = '../../02_r2SCAN-D4/c.train/a.16g_iter_no_cont/trained.almtp'
# mtp_s = sorted(glob('../../a.merge_MTP/d.train/a.cut_5/calc/*/trained.almtp'))
mtp_s = sorted(glob('/home/atuin/a102cb/a102cb12/SiO2_mesoporous_final_MTPs/c.r2SCAN_D4_syn_MTP/trained.almtp'))
# train_cfg = '../../02_r2SCAN-D4/a.add_D4/d.collect/struct_new.cfg'

# select_start_seed = 0
# # select_end_seed = 15
# # select_end_seed = 31
# # select_end_seed = 63
# select_end_seed = 127  # out of memory in fritz
# # select_end_seed = 256

variable_file = 'in.header'  # substitute the variables
