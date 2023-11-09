#!/bin/bash
#PBS -A NUMLIB
#PBS -q gpu
#PBS -l elapstim_req=02:00:00

cd $PBS_O_WORKDIR
../build/bench/SixstepTime 24 1 2 4 8 16 32 48 96