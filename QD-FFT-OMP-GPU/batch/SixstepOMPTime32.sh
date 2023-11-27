#!/bin/bash
#PBS -A NUMLIB
#PBS -q gpu
#PBS -l elapstim_req=02:00:00

cd $PBS_O_WORKDIR
../build/bench/SixstepOMPTime 28 4 16 32 48 96 --range -loops 10 -start 26