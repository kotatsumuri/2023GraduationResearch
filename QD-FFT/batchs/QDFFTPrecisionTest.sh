#!/bin/bash
#PBS -A NUMLIB
#PBS -q pmem-b
#PBS -l elapstim_req=00:10:00

cd $PBS_O_WORKDIR
../build/test/QDFFTPrecisionTest