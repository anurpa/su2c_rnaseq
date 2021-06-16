#!/bin/bash 
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe orte 24
#$ -V

tophat-1.3.0.Linux_x86_64/tophat -p 24 -r 100 -o outputFileName -G path/to/UCSC_Genes.gtf path/to/hg19/bowtieFiles re
ads1.fastq reads2.fastq
