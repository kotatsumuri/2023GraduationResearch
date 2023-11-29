#!/bin/bash
#PBS -A NUMLIB
#PBS -q gpu
#PBS -l elapstim_req=02:00:00

cd $PBS_O_WORKDIR
../build/bench/SixstepGPUTime 30 --range -loops 100