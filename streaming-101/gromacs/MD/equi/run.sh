#!/bin/bash

gmx grompp -f ../mdp-files/equi-NPT.mdp -c ../em/confout.gro -p ../sys/topol.top -o -maxwarn 2
gmx mdrun -v -nt 16
