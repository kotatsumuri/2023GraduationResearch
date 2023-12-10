#!/bin/bash
#PBS -A NUMLIB
#PBS -q gpu
#PBS -l elapstim_req=02:00:00

cd $PBS_O_WORKDIR
../build/bench/SixstepOMPTime 29 5 16 32 48 64 96 -loops 10