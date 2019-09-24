#!/bin/bash
pwd

#source ~/.bash_profile;

#22 runs. First 11 are for Exp 1, second 11 are for Exp 2

for R in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do

#submit sampling job
job1=$(qsub -v RUN=$R -o /QRISdata/Q0992/analysis/Bayes/ -e /QRISdata/Q0992/analysis/Bayes/ /QRISdata/Q0992/analysis/Bayes/1.2_collect_samples.pbspro)

#wait for sampling job to complete and submit compile job
qsub -W depend=afterok:$job1 -v RUN=$R -o /QRISdata/Q0992/analysis/Bayes/ -e /QRISdata/Q0992/analysis/Bayes/ /QRISdata/Q0992/analysis/Bayes/1.4_compile_samples.pbspro

#qsub -v RUN=$R -o /QRISdata/Q0992/analysis/Bayes/ -e /QRISdata/Q0992/analysis/#Bayes/ /QRISdata/Q0992/analysis/Bayes/1.4_compile_samples.pbspro

done
