import numpy as np
import pandas as pd
import sys
import time
import os
import io
import csv
import math
import datetime
import csv
import re

from sklearn.utils import shuffle
from sklearn.metrics import roc_auc_score, roc_curve, auc
from scipy import stats

# from keras import initializers
from keras.models import load_model
from tensorflow.keras.utils import CustomObjectScope
import tensorflow.python.util.deprecation as deprecation
from attention import Attention
import multiprocessing as mp
import tensorflow as tf

tf.get_logger().setLevel('ERROR')   # warning
deprecation._PRINT_DEPRECATION_WARNINGS = False  # future discard
tf.compat.v1.disable_v2_behavior()  # ValueError: GRU(reset_after=False) is not compatible with GRU(reset_after=True)

os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
curDir = os.path.dirname(os.path.realpath(__file__))+'/'
HLA_seq = pd.read_csv(curDir + 'model/MHC_pseudo.dat', sep='\t')
seqs = {}
for i in range(len(HLA_seq)):
    seqs[HLA_seq.HLA[i]] = HLA_seq.sequence[i]

predScores = []
predScores1 = []
aa_idx = {'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 9, 'L': 10, 'M': 11,
          'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17, 'V': 18, 'W': 19, 'Y': 20, 'X': 21}


def run_model(i, X_test):
    score = np.zeros((5, len(X_test)))
    with CustomObjectScope({'Attention': Attention}):
        model = load_model(
            curDir + 'model/binding_model' + str(i+1) + '.hdf5')
        score[i, :] = np.squeeze(model.predict(X_test))
    return score[i, :]


def run_model1(i, X_test):
    score1 = np.zeros((5, len(X_test)))
    with CustomObjectScope({'Attention': Attention}):
        model1 = load_model(
            curDir + 'model/immunogenicity_model' + str(i+1) + '.hdf5')
        score1[i, :] = np.squeeze(model1.predict(X_test))
    return score1[i, :]


def collect_result(result):
    global predScores
    predScores.append(result)


def collect_result1(result1):
    global predScores1
    predScores1.append(result1)


def transform(HLA, peptide):
    data = HLA+peptide
    seq = data+'X'*(49-len(data))
    seq = [aa_idx[x] for x in seq]
    return seq


def read_and_prepare(file):
    data = pd.read_csv(file)
    data['cost_cents'] = data.apply(
        lambda row: transform(
            HLA=seqs[row['HLA']],
            peptide=row['peptide']),
        axis=1)
    return np.vstack(data.cost_cents)


def read_and_prepare_single(peptide, hla):
    complex = np.full((1, 49), 21, int)
    seq = [aa_idx[x] for x in list(seqs[hla] + peptide)]
    for i in range(len(seq)):
        complex[0, i] = seq[i]
    return complex


def deephlapan_main(opt):
    i = datetime.datetime.now()
    print(str(i) + ' Prediction starting.....\n')
    # print("="*30,"GPU USED","="*30,"tf version:",tf.__version__,"use GPU:",tf.test.is_gpu_available())

    peptide = opt.sequence
    hla = opt.hla
    dir = opt.WD
    if len(dir) == 0:
        dir = '.'

    fname = peptide+'_'+hla
    if (opt.file):
        fname = opt.file.split('/')[-1]
        fname = fname.split('.')[0]
        df = pd.read_csv(opt.file)
        X_test = read_and_prepare(opt.file)
    else:
        X_test = read_and_prepare_single(peptide, hla)
    
    pool = mp.Pool(mp.cpu_count())
    for i in range(5):
        pool.apply_async(run_model, args=(i, X_test), callback=collect_result)
        pool.apply_async(run_model1, args=(i, X_test),
                         callback=collect_result1)
    pool.close()
    pool.join()

    result = np.average(predScores, axis=0)
    result1 = np.average(predScores1, axis=0)
    with open(dir + '/' + fname.replace(':', '') + '_predicted_result.csv', 'w') as f:
        f.write('Seq_id,Allele,Position,Peptide,Bindingscore,Immunogenicscore\n')
        if (opt.file):
            for i in range(len(result)):
                result[i] = ("%.4f" % result[i])
                result1[i] = ("%.4f" % result1[i])
                f.write(str(df.seq_id[i]) + ',' + str(df.HLA[i]) + ',' + str(df.position[i]) + ',' + str(
                    df.peptide[i]) + ',' + str(format(result[i],'.4f')) + ',' + str(format(result1[i],'.4f')) + '\n')
        else:
            f.write('single peptide,' + str(hla) + ',' + str(peptide) +
                    ',' + str(result[0]) + ',' + str(result1[0]) + '\n')

    if (opt.file):
        command = 'perl ' + curDir + 'model/rank.pl ' + \
            dir + '/' + fname + '_predicted_result.csv'
        os.system(command)
    j = datetime.datetime.now()
    print(str(j) + ' Prediction end\n')
