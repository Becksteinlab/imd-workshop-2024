#!/bin/bash

source /usr/local/gromacs-imd-v3/bin/GMXRC

mkdir -p imd
cd imd

gmx grompp -f ../mdp-files/run-NPT_imd-v3.mdp -c ../equi/confout.gro -p ../sys/topol.top -o 
gmx mdrun -v -nt 4 -imdwait -imdport 8889

cd ..
