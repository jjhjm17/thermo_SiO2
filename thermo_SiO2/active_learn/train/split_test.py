from thermo_SiO2.io import read_sil
from thermo_SiO2.mlip.write_config import write_cfg_SiO2
import random

def split_test(cfg_file, symbols, val_fraction):
    """This function splits train and test sets.
    Example
    cfg_file = '/home/atuin/a102cb/a102cb12/SiO2_mesoporous/b.MTP/o.final_r2SCAN_stoichiometry/a.train/a.train_set_r2SCAN-D4/stoi.cfg'
    symbols = 'Si O H Al'
    val_fraction = 0.05  # validation fraction
    """

    # Read configurations
    cfgs = read_sil(cfg_file, atom_symbols=symbols)
    # cfgs = cfgs[:300]

    # --- Split parameters ---
    RANDOM_SEED = 42
    VAL_FRACTION = val_fraction

    # Shuffle indices with fixed seed for reproducibility
    indices = list(range(len(cfgs)))
    random.seed(RANDOM_SEED)
    random.shuffle(indices)

    # Split into validation and training sets
    n_val = max(1, int(len(cfgs) * VAL_FRACTION))
    val_indices   = indices[:n_val]
    train_indices = indices[n_val:]

    val_cfgs   = [cfgs[i] for i in sorted(val_indices)]
    train_cfgs = [cfgs[i] for i in sorted(train_indices)]

    print(f"Total configurations : {len(cfgs)}")
    print(f"Training set         : {len(train_cfgs)}  -> train_95pro.cfg")
    print(f"Validation set       : {len(val_cfgs)}   -> test_5pro.cfg")

    # Save to files
    write_cfg_SiO2(file='train_95pro.cfg', configs=train_cfgs, atom_symbols=symbols)
    write_cfg_SiO2(file='test_5pro.cfg',    configs=val_cfgs,   atom_symbols=symbols)

if __name__ == '__main__':
    main()
