from a_parameters import server

if server == 'Fritz':
    # species_count = 2  # SiO2

    # Justus 2 server
    # MLP=mlp_par_fits_to_energy_based_on_ismear
    # module load compiler/intel
    # module load numlib/mkl

    # Fritz server at the University of Erlangen-Nuernberg.
    # A parallel version does not run at the login node.
    mlp = 'mlp.2020-07.serial'

    print("\nPlease run the following line at the terminal by hand.")
    print("module load    intel  mkl  intelmpi\n")


    # Justus 2 server
    # init_mtp_folder = '/home/st/st_us-031400/st_ac138035/program/mlip-2/untrained_mtps/'
    # Fritz server at the University of Erlangen-Nuernberg
    init_mtp_folder = '/home/hpc/a102cb/a102cb12/program/mlip-2/untrained_mtps/'
elif server == 'Justus2':
    mlp = 'mlp3'
    print("\nPlease run the following line at the terminal by hand.")
    print("module load compiler/intel/19.0   numlib/mkl/2019   mpi/impi/2019.5\n")
    init_mtp_folder = '/home/st/st_us-031400/st_ac138035/program/mlip-3/MTP_templates/'

