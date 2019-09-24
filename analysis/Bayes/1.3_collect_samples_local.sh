#!/bin/bash

# cd cmdstan

echo $1
echo $2

MODEL_ARRAY=(hyperbolic exponential hyperbolic_gm 
  prop_diff tradeoff ITCH const_sens mazur1987 loewenstein1992
  mcclure2007 killeen2009)

if [ "$2" -lt 12 ]
then
	MODEL_NUM="$2"-1
	EXP=1
fi

if [ "$2" -gt 11 ]
then
	MODEL_NUM="$2"-12
	EXP=2
fi

echo $EXP
echo $MODEL_NUM
echo ${MODEL_ARRAY[MODEL_NUM]}

THREADS=-1

export STAN_NUM_THREADS=$THREADS

~/Nextcloud/DISCRECOV-Q0992/models/"${MODEL_ARRAY[MODEL_NUM]}"_mpi sample num_samples=1000 num_warmup=1000 thin=1 id=$1 data file=~/Nextcloud/DISCRECOV-Q0992/data/clean/data_list_e"${EXP}"_rdump.R random seed=12345 output file=~/Nextcloud/DISCRECOV-Q0992/analysis/Bayes/samples_e"${EXP}"_"${MODEL_ARRAY[MODEL_NUM]}"_$1.csv > ~/Nextcloud/DISCRECOV-Q0992/analysis/Bayes/samples_e"${EXP}"_"${MODEL_ARRAY[MODEL_NUM]}"_$1.txt

