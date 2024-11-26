#!/bin/bash

gmx grompp -f ../mdp-files/em.mdp -c ../sys/conf.gro -p ../sys/topol.top -o 
gmx mdrun -v -nt 1
