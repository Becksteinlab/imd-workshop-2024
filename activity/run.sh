#!/bin/bash

cd sample_simulation

## RUN GROMPP
if [ ! -f sample_simulation/topol.tpr ]; then
echo "Running grompp"
gmx grompp -f imd.mdp -c start.gro -p topol.top -o >& grompp.out
fi
echo "Starting mdrun"

## START MDRUN
gmx mdrun -v -nt 1 -imdwait -imdport 9999
cd ..
