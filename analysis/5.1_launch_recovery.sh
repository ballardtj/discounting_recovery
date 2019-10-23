#!/bin/bash
pwd

#source ~/.bash_profile;

#submit sampling job
job1=$(qsub /30days/uqtballa/analysis/5.2_TD_collect_samples.pbspro)

#wait for sampling job to complete and submit compile job
qsub -W depend=afterok:$job1 -v /30days/uqtballa/Q0992/analysis/5.5_compile_samples.pbspro

done
