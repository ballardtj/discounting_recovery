#!/bin/bash

echo 5.6
echo ${1}

#Runs the R script that compiles the samples
Rscript /30days/uqtballa/Q0992/analysis/5.7_compile_samples.R --${1} > /30days/uqtballa/Q0992/analysis/output/5.7_compile_samples_${1}.Rout 2>&1