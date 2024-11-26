#!/bin/bash

files=(
traj.xyz
log.lammps
dump_initial_config.dump
topology.data
logfile.txt
)

for f in ${files[@]}
do
if [ -f $f ]; then
rm $f
fi
done
