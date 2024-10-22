#!/bin/bash

echo "setting GROMACS version"
source /usr/local/gromacs-imd-v3/bin/GMXRC

mkdir -p imd
cd imd
if [ ! -f topol.tpr ]; then
echo "running grompp"
gmx grompp -f ../imd.mdp -c ../start.gro -p ../topol.top -o >& grompp.out
fi
echo "starting mdrun"
echo ""
echo "**************************************************************"
echo "wait a few seconds, then:"
echo "run (in separate terminal): python demo.py"
echo "**************************************************************"
echo ""
# modify gromacs execution string according to available resources
# here we use 4 parallel threads
gmx mdrun -v -nt 4 -imdwait -imdport 8888
cd ..
