#!/bin/bash

# Fritz server

# submit.sh ~/Thermodynamics/run_scripts/run.lammps.mlp3.fast.twoSteps.par.cluster.fritz.2tb.n1t1  -j jobList  | tee c.submitJob.out
submit.sh ~/Thermodynamics/run_scripts/run.lammps.mlp3.fast.twoSteps.par.cluster.fritz.1tb.n1t1  -j jobList  | tee c.submitJob.out

# justus
# submit.sh ~/Thermodynamics/run_scripts/run.lammps.mlp3.par.cluster.justus2.c48t4  | tee c.submitJob.out

exit
