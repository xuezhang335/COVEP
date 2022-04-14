import os, argparse


def command_line():
    parser = argparse.ArgumentParser()
    parser.add_argument('-bl', dest='bl', action='store_true',
                        help='Prediction for B cell linear epitopes.')
    parser.add_argument('-t1', dest='t1', action='store_true',
                        help='Prediction for T cell epitopes presented by MHC1 molecular.')
    parser.add_argument('-t2', dest='t2', action='store_true',
                        help='Prediction for T cell epitopes presented by MHC2 molecular.')
    parser.add_argument('-s', '--seqfile', dest='seqfile',
                        help='txt file contains one or more fasta sequences.')
    parser.add_argument('-hla', dest='hlastr',
                        help='This is necessary if t1 or t2 is selected.')
    parser.add_argument('-hlafile', dest='hlafile',
                        help='This is necessary if t1 or t2 is selected.')
    parser.add_argument('-o', '--outdir', dest='outdir', 
                        help='Directory to store predicted results. User must have write privilege. If omitted, the current directory will be applied.')

    return parser.parse_args()



