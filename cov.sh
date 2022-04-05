#! /bin/bash

seqfile=/root/cov/test/data/2019-ncov-test.txt
hlafile=/root/cov/test/data/hla_allele-test.csv
outdir=/root/cov/test/output

python /root/cov/scripts/main.py -bl -t1 -t2 -s $seqfile -hlafile $hlafile -o $outdir