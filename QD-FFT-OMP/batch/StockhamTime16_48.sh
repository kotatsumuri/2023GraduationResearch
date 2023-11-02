#!/bin/bash
#PBS -A NUMLIB
#PBS -q gpu
#PBS -l elapstim_req=00:10:00

cd $PBS_O_WORKDIR
../build/bench/StockhamTime 30 48