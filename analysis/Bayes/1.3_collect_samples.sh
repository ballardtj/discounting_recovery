#!/bin/bash

# cd cmdstan

echo $1
echo $2

MODEL_ARRAY=(hyperbolic_nontruncated exponential_nontruncated hyperbolic_gm 
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

/QRISdata/Q0992/models/"${MODEL_ARRAY[MODEL_NUM]}"_mpi sample algorithm=hmc engine=nuts max_depth=20 num_samples=5000 num_warmup=10000 thin=1 adapt delta=.975 data file=/QRISdata/Q0992/data/clean/data_list_e"${EXP}"_rdump.R random seed=12345 id=$1 output file=/QRISdata/Q0992/analysis/Bayes/samples_e"${EXP}"_"${MODEL_ARRAY[MODEL_NUM]}"_$1.csv > /QRISdata/Q0992/analysis/Bayes/samples_e"${EXP}"_"${MODEL_ARRAY[MODEL_NUM]}"_$1.txt

