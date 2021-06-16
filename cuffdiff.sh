#!/bin/bash 
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe orte 10
#$ -V

cufflinks-1.3.0.Linux_x86_64/cuffdiff -L sample1_OutputName sample2_OutputName -
p 10 -o outputFileName path/to/UCSC_Genes.gtf path/to/sample1/tophatOutput/accep
ted_hits.bam path/to/sample2/tophatOutput/accepted_hits.bam
