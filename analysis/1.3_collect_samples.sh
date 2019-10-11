echo 1.3
#echo ${1}
#echo ${2}

#Runs the R script that compiles the samples
Rscript /30days/uqtballa/Q0992/analysis/1.4_collect_samples.R --${1} --${2} > /30days/uqtballa/Q0992/analysis/output/1.4_collect_samples_${1}_${2}.Rout 2>&1

#Rscript ~/Nextcloud/DISCRECOV-Q0992/analysis/1.4_collect_samples.R --1 --1 > ~/Nextcloud/DISCRECOV-Q0992/analysis/output/1.4_collect_samples_1_2.Rout 2>&1