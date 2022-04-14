import os
import re
import time
import numpy as np
import pandas as pd

from Bio import SeqIO

# import bepipred2 as bp2


cov_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def iedb(task_name, seq_file, outdir):
    seq_list = []
    methods = ['Chou-Fasman', 'Emini', 'Karplus-Schulz',
               'Kolaskar-Tongaonkar', 'Parker', 'Bepipred', 'Bepipred-2.0']  # 'Bepipred-2.0'
    id_list = []

    def cla_threshold(method, score):
        if method == 'Bepipred':
            threshold = 0.35
        elif method == 'Bepipred-2.0':
            threshold = 0.50
        elif method == 'Emini':
            threshold = 1.00
        else:
            threshold = score.astype(float).mean(0)
        return threshold

    records = SeqIO.parse(seq_file, 'fasta')
    for record in records:
        seq_id = str(record.id)
        seq = str(record.seq)
        id_list.append(seq_id)
        seq_list.append(seq)

    process_iedb = {}
    for method in methods:
        process_iedb[method] = True
        try:
            outfile = f'{outdir}/cov_tools_predResult/bl_{task_name}_{method}.csv'
            print('='*50, f' RUNNING IEDB {method}... ', '='*50)

            with open(outfile, 'w')as f:
                f.write('Seq_id,Position,Residue,Score,pred_label\n')

            for seq_id, seq in zip(id_list, seq_list):
                with open(outfile, 'a')as f:
                    outfile_tmp = f'{outdir}/tmp/bl_{task_name}_{method}.csv'
                    os.system(
                        f'curl --data "method={method}&sequence_text={seq}" http://tools-cluster-interface.iedb.org/tools_api/bcell/ > {outfile_tmp}')

                    df = pd.read_csv(outfile_tmp)
                    column_name = df.columns.values[0].split()
                    result = df.iloc[:, 0].str.split(expand=True)
                    result.columns = column_name

                    start_pos = int(result.loc[0, 'Position'])
                    end_pos = int(result.loc[len(result)-1, 'Position'])

                    if start_pos != 1:
                        for i, res in enumerate(seq[:start_pos-1]):
                            f.write(seq_id + ',' + str(i+1) + ',' +
                                    res + ',' + str(0) + ',' + str(0) + '\n')

                    result['Seq_id'] = [seq_id]*len(result['Score'])
                    threshold = cla_threshold(method, result['Score'])
                    result['pred_label'] = result['Score'].apply(
                        lambda x: 1 if float(x) >= threshold else 0)

                    for index in result.index.values:
                        f.write(seq_id + ',' + result.loc[index, 'Position'] + ',' + result.loc[index, 'Residue'] + ',' + str(
                            result.loc[index, 'Score']) + ',' + str(result.loc[index, 'pred_label']) + '\n')

                    if end_pos != len(seq):
                        for i, res in enumerate(seq):
                            if i >= end_pos:
                                f.write(seq_id + ',' + str(i+1) + ',' +
                                        res + ',' + str(0) + ',' + str(0) + '\n')
        except:
            process_iedb[method] = False
            print(f'''Warning: Due to network reasons, the prediction of B cell epitope by {method} is not completed,
                        and this method prediction result analysis will not be carried out later.''')

    return process_iedb


# def bepipred2(task_name, seq_file, outdir):
#     seq_file= cov_dir + '/data/2019-ncov-test.txt'
#     os.system(f'conda run -n base BepiPred-2.0 {seq_file}')
#     return


def runEpitopevec(task_name, seq_file, outdir):
    print('='*50, ' RUNNING EpitopeVec... ', '='*50)
    print('In progress.....')
    outfile_tmp = f'{outdir}/tmp/bl_{task_name}_epitopevec.csv'
    os.system(
        f'conda run -n cov python {cov_dir}/tools/epitopevec/main.py -i {seq_file} -o {outfile_tmp}')

    outfile = f'{outdir}/bl_{task_name}_epitopevec.csv'
    df = pd.read_csv(outfile_tmp)
    df['pred_label'] = df['Score'].apply(lambda x: 1 if x >= 0.5 else 0)
    df.to_csv(outfile, index=False)
    return
