#!/bin/bash

echo 1.5
echo ${1}

#Runs the R script that compiles the samples
Rscript /QRISdata/Q0992/analysis/Bayes/1.6_compile_samples.R --${1}
