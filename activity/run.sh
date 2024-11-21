#!/bin/bash

mkdir -p sample_simulation
cd sample_simulation
# if [ ! -f sample_simulation/topol.tpr ]; then
# echo "Running grompp"
# gmx grompp -f imd.mdp -c start.gro -p topol.top -o >& grompp.out
# fi
echo "Starting mdrun"
# modify gromacs execution string according to available resources
# here we use 2 parallel threads
gmx mdrun -v -nt 2 -imdwait -imdport 8889
cd ..
