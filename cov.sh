#! /bin/bash

curdir=$(cd $(dirname $0) && pwd)

seqfile=${curdir}/test/data/2019-ncov-test.txt
hlafile=${curdir}/test/data/hla_allele-test.csv
outdir=${curdir}/test/output

python ${curdir}/scripts/main.py -t1 all -t2 all -s $seqfile -hlafile $hlafile -o $outdir