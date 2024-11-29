#!/bin/bash

## TMP
source /usr/local/gromacs/bin/GMXRC

## WAIT UNTIL PORT IS OPEN
SECONDS=0
while lsof -n -i :8889; do
    echo "Waiting for port 8889 to open"
    sleep 1

    if [ $SECONDS -ge 30 ]; then
        echo "Timeout reached. Port 8889 did not open."
        exit 1
    fi
done

mkdir -p sample_simulation
cd sample_simulation

## RUN GROMPP
if [ ! -f sample_simulation/topol.tpr ]; then
echo "Running grompp"
gmx grompp -f imd.mdp -c start.gro -p topol.top -o >& grompp.out
fi
echo "Starting mdrun"

## START MDRUN
gmx mdrun -v -nt 1 -imdwait -imdport 8889
cd ..
