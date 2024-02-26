#!/bin/bash
#PBS -A NUMLIB
#PBS -q gpu
#PBS -l elapstim_req=02:00:00

module load intel

cd $PBS_O_WORKDIR
../build/bench/StockhamOMPBenchFFTSpeed 28 -start 10