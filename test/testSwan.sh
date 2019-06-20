#!/bin/bash
#PBS -l select=1:ncpus=44
#PBS -l place=scatter
#PBS -l walltime=04:00:00
#PBS -V

cd $PBS_O_WORKDIR
start_test 
