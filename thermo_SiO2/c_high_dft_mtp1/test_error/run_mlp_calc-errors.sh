#!/bin/bash
# In python script, mlp cannot be directly called using subprocess.
# ex) run_mlp_calc-errors.sh 

MLP='mlp.2020-07.serial'  # Fritz server
$MLP calc-errors "$1"  "$2"
