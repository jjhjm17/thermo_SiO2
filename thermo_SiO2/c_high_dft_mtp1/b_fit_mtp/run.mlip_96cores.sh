#!/usr/bin/env bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=48
# #SBATCH --ntasks-per-node=2  # debug
#SBATCH --time=24:00:00
#SBATCH --mem=16gb  # large memory
#SBATCH -J MTP_fitting
#SBATCH --output=MTP_fitting.out
#SBATCH --error=MTP_fitting.err

MLP=mlp_par_fits_to_energy_based_on_ismear
# MLP=mlp
module load compiler/intel
module load numlib/mkl
module load mpi/impi
SRUN_COMMAND=srun
# SRUN_COMMAND=""  # debug

set -u  # stop at undefined variable

# srun ~/executables/mlp_cpar train init.mtp training_set.cfg
# srun $MLP train init.mtp training_set.cfg 
ln -s ../training_set.cfg .
ln -s ../allcfgs.cfg .

mkdir mlp_a_every_10
cd mlp_a_every_10
mv ../init.mtp .
ln -s ../training_set.cfg .

$SRUN_COMMAND $MLP train --weighting=structures init.mtp training_set.cfg | tee mlp.train.out

get_result () {
    rmsEnergy=`grep --after-context=4 "Energy per atom" mlp.train.out | tail -n1 | awk '{ print $5 }'`
    rmsForce=`grep --after-context=4 "Forces" mlp.train.out | tail -n1 | awk '{ print $5 }'`
    mtpLevel=`cat ../mtp_level`
    echo "# Level of MTP   RMS of energy per atom (eV/atom)    RMS of force (eV/Ang)" > result.txt
    echo "$mtpLevel  $rmsEnergy  $rmsForce" >> result.txt
    cp result.txt ..
}

get_result 
cd ..


prevFolder=mlp_a_every_10
newFolder=mlp_b_select_add

select_and_train () {
    mkdir $newFolder
    cd $newFolder
    # $SRUN_COMMAND $MLP select-add --weighting=structures ../$prevFolder/Trained.mtp_ ../$prevFolder/training_set.cfg ../allcfgs.cfg to_be_added.cfg | tee mlp.select-add.out
    # select-add is for single core. Otherwise it is repeated with the number
    # of cores.
    srun -n 1 $MLP select-add --weighting=structures ../$prevFolder/Trained.mtp_ ../$prevFolder/training_set.cfg ../allcfgs.cfg to_be_added.cfg | tee mlp.select-add.out
    numConfig=`grep Size to_be_added.cfg | wc -l`
    if [ $numConfig = "0" ]; then
        echo "The number of configurations in to_be_added.cfg is 0."
        echo "The training is finished."
        exit
    fi
    cat ../$prevFolder/training_set.cfg to_be_added.cfg > training_set.cfg
    rm to_be_added.cfg
    $SRUN_COMMAND $MLP train --weighting=structures ../$prevFolder/Trained.mtp_ training_set.cfg | tee mlp.train.out
    rm ../last_mtp_in_$prevFolder
    touch  ../last_mtp_in_$newFolder

    get_result 
    cd ..
    sleep 5  # so that Trained.mtp_ is written.
}

select_and_train


prevFolder=mlp_b_select_add
newFolder=mlp_c_select_add

select_and_train


prevFolder=mlp_c_select_add
newFolder=mlp_d_select_add

select_and_train

