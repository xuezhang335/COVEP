import os
import pandas as pd
from Bio import SeqIO


def creatPepForHla1(seq_file, outdir):
    records = SeqIO.parse(seq_file, 'fasta')
    pep_list = []  # seq-->8-11mer
    id_list = []   # id
    pos_list = []  # position
    for record in records:
        seq_id = str(record.id)
        seq = str(record.seq)
        for i in range(len(seq)):
            for j in range(8, 12):
                if i+j <= len(seq):
                    pep_list.append(seq[i:i+j])
                    pos_list.append(str(i+1) + '_' + str(i+j+1-1))
                    id_list.extend([seq_id])

    pep_for_hla1 = pd.DataFrame(
        {'seq_id': id_list, 'peptide': pep_list, 'pos': pos_list})
    # outfile_pep_for_hla1 = seq_file.split('.')[0]+'_pep_for_hla1.csv'
    outfile_pep_for_hla1 = outdir + '/tmp/' + os.path.basename(seq_file).split('.')[0] + '_pep_for_hla1.csv'
    pep_for_hla1.to_csv(outfile_pep_for_hla1, index=False)
    return outfile_pep_for_hla1


def creatPepForHla2(seq_file, outdir):
    records = SeqIO.parse(seq_file, 'fasta')
    pep_list = []  # seq-->15-25mer
    id_list = []   # id
    pos_list = []  # position
    for record in records:
        seq_id = str(record.id)
        seq = str(record.seq)
        for i in range(len(seq)):
            for j in range(12, 22):
                if i+j <= len(seq):
                    pep_list.append(seq[i:i+j])
                    pos_list.append(str(i+1) + '_' + str(i+j+1-1))
                    id_list.extend([seq_id])

    pep_for_hla2 = pd.DataFrame(
        {'seq_id': id_list, 'peptide': pep_list, 'pos': pos_list})
    # outfile_pep_for_hla2 = seq_file.split('.')[0]+'_pep_for_hla2.csv'
    outfile_pep_for_hla2 = outdir + '/tmp/' + os.path.basename(seq_file).split('.')[0] + '_pep_for_hla2.csv'
    pep_for_hla2.to_csv(outfile_pep_for_hla2, index=False)
    return outfile_pep_for_hla2
