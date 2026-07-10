from collections import Counter

def get_E_form(cfg, ref_E):
    # ref_E: 
    # export PYTHONPATH=${PYTHONPATH}:/home/atuin/a102cb/a102cb12/SiO2_mesoporous/b.MTP/a.1.3_ENMAX/e.struct/b.ref_E_cx
    #
    # from ref_E import ref_E
    #
    # copied from 
    # /home/atuin/a102cb/a102cb12/SiO2_mesoporous/b.MTP/e.SiOHAl_1000K/a.pore-eq-MTP_vdW/04_pore_OH_v2_relax/e.high_DFT_relax_train/a.subsample_train/c.train/b.error_on_subset/plot_error_and_E_form.py
    # 
    energy = cfg.energy  # eV/cell
    symbols_counter = Counter(cfg.get_chemical_symbols())
    num_SiO2 = symbols_counter['Si']
    num_Al2O3 = symbols_counter['Al'] / 2
    num_H2O = symbols_counter['H'] / 2
    num_O2 = (symbols_counter['O']
              - num_SiO2 * 2
              - num_Al2O3 * 3
              - num_H2O ) / 2

    # formation energy
    # /home/atuin/a102cb/a102cb12/SiO2_mesoporous/b.MTP/a.1.3_ENMAX/e.struct/b.ref_E_cx
    E_f_ref = (num_SiO2 * ref_E['SiO2'] 
               + num_Al2O3 * ref_E['Al2O3'] 
               + num_H2O * ref_E['H2O'] 
               + num_O2 * ref_E['O2'])
    # formation E ref
    return (energy - E_f_ref) / len(cfg)  # eV/atom


# Claude
def net_formal_charge(cfg, formal_charge_dic=None):
    if formal_charge_dic is None:
        formal_charge_dic = {'Si': 4, 'O': -2, 'H': 1, 'Al': 3}
    symbols_counter = Counter(cfg.get_chemical_symbols())
    return sum(symbols_counter[el] * charge
               for el, charge in formal_charge_dic.items())
