#!/bin/bash
#PBS -A NUMLIB
#PBS -q gpu
#PBS -l elapstim_req=02:00:00

cd $PBS_O_WORKDIR
../build/bench/StockhamGPUBenchFFTSpeed 28 -start 10