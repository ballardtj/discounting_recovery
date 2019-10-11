#!/bin/bash

echo 1.6
echo ${1}

#Runs the R script that compiles the samples
Rscript /30days/uqtballa/Q0992/analysis/1.7_compile_samples.R --${1} > /30days/uqtballa/Q0992/analysis/output/1.7_compile_samples_${1}.Rout 2>&1

#Rscript ~/Nextcloud/DISCRECOV-Q0992/analysis/1.7_compile_samples.R --1 > ~/Nextcloud/DISCRECOV-Q0992/analysis/output/1.7_compile_samples_1.Rout 2>&1