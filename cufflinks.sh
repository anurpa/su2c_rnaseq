#!/bin/bash 
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe orte 10
#$ -V

cufflinks-1.3.0.Linux_x86_64/cufflinks -p 10 -o outputFileName -G path/to/UCSC_Genes.gtf  path/to/tophatOutput/accepted_hits.bam
