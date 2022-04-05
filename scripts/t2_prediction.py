import os
import numpy as np
import pandas as pd


cov_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def runMixmhc2pred(task_name, pep2_file, input_hla, outdir):
    print('='*50, ' RUNNING MixMHC2pred... ', '='*50)
    df = pd.read_csv(pep2_file)
    pep_list = list(df.peptide)
    id_list = list(df.seq_id)
    pos_list = list(df.pos)

    infile = os.path.dirname(
        pep2_file) + '/' + '_'.join(os.path.basename(pep2_file).split('_')[:-1]) + '_mixmhc2pred.txt'
    with open(infile, 'w')as f:
        for i in pep_list:
            f.write(i + '\n')

    if type(input_hla)==type(''):
        df = pd.read_csv(input_hla)
        hla2 = list(df.hla2.dropna())
    else:
        hla2 = input_hla

    str_hla2 = ' '.join(hla2)
    outfile_tmp = outdir + '/tmp/mixmhc2pred_orgin.txt'

    os.system(
        f'{cov_dir}/tools/MixMHC2pred-1.2/MixMHC2pred_unix -i {infile} -o {outfile_tmp} -a {str_hla2}')

    lines = (i for i in open(outdir + '/tmp/mixmhc2pred_orgin.txt')
             if str(i)[0] != '#')
    f = open(outdir + '/tmp/mixmhc2pred_orgin(1).txt', 'w', encoding='utf-8')
    f.writelines(lines)
    f.close()

    result_txt = np.loadtxt(
        outdir + '/tmp/mixmhc2pred_orgin(1).txt', dtype=str)
    result_txtdf = pd.DataFrame(result_txt[1:], columns=result_txt[0])
    result_txtdf.insert(0, 'Position', pos_list)
    result_txtdf.insert(0, 'Seq_id', id_list)

    result = {'Seq_id': [], 'Allele': [],
              'Position': [], 'Peptide': [], '%Rank': []}
    for key in ['Seq_id', 'Position', 'Peptide']:
        for i in result_txtdf[key]:
            result[key].extend([i]*len(hla2))

    result['Allele'] = hla2*len(result_txtdf['Peptide'])
    rank_clomun = ['%Rank_' + hla for hla in hla2]
    for i in list(result_txtdf.index):
        for c in rank_clomun:
            result['%Rank'].append(result_txtdf.loc[i, c])
    result = pd.DataFrame(result)
    result['pred_label'] = result['%Rank'].map(lambda x:1 if float(x)<2 else 0)
    result.to_csv(
        outdir + f'/cov_tools_predResult/t2_{task_name}_mixmhc2pred.csv', index=False)

    return


def runNetmhc2pan(task_name, pep2_file, input_hla, outdir):
    print('='*50, ' RUNNING NetMHC2pan... ', '='*50)
    print('In progress.....')
    df = pd.read_csv(pep2_file)
    pep_list = list(df.peptide)
    id_list = list(df.seq_id)
    pos_list = list(df.pos)

    if type(input_hla)==type(''):
        df = pd.read_csv(input_hla)
        hla2 = list(df.hla2.dropna())
    else:
        hla2 = input_hla

    hla2_net = []
    for hla in hla2:
        if hla[0:2] == 'DR':
            hla = hla.split('_')[0] + '_' + \
                hla.split('_')[1] + hla.split('_')[2]
            hla2_net.append(hla)
        else:
            hla = 'HLA-' + hla
            hla = hla.replace('__', '-')
            hla = hla.replace('_', '')
            hla2_net.append(hla)
    str_hla2 = ','.join(hla2_net)

    infile = outdir + '/tmp/' + '_'.join(os.path.basename(pep2_file).split('_')[:-1]) + '_netmhc2pan.txt'
    outfile_tmp = outdir + '/tmp/netmhc2pan_out.xls'
    terminalout = outdir + '/tmp/netmhc2pan_out.print'

    with open(infile, 'w')as f:
        for p in pep_list:
            f.write(p + '\n')

    os.system(
        f'conda run -n cov {cov_dir}/tools/netMHCIIpan-4.0/netMHCIIpan -f {infile} -inptype 1 -BA -xls -a {str_hla2} -xlsfile {outfile_tmp} > {terminalout}')  # 

    processNetmhc2panOutput(
        task_name, pep_list, id_list, pos_list, hla2, outdir)
    
    return


def processNetmhc2panOutput(task_name, pep_list, id_list, pos_list, hla2, outdir):
    hla2_sort = sorted(hla2)
    hla_pos_dic = {}
    for hla in hla2:
        pos = hla2_sort.index(hla)
        hla_pos_dic[hla] = pos
    outfile_tmp = outdir + '/tmp/netmhc2pan_out.xls'
    predresult = pd.read_csv(
        outfile_tmp, sep='\t', names=list(range(4+len(hla2)*6)))  # low_memory=False
    outfile = outdir + '/cov_tools_predResult/t2_' + task_name + '_netmhc2pan.csv'
    result = {'Seq_id': [], 'Allele': [], 'Position': [], 'Peptide': [], 'Rank': []}
    for id, pep in enumerate(pep_list):
        for hla in hla_pos_dic:
            result['Seq_id'].append(id_list[id])
            result['Allele'].append(hla)
            result['Position'].append(pos_list[id])
            result['Peptide'].append(pep)
            # result['Score'].append(
            #     predresult.iloc[id+2, 4+hla_pos_dic[hla]*5])
            result['Rank'].append(
                predresult.iloc[id+2, 4+hla_pos_dic[hla]*5+1])
            # result['Score_BA'].append(
            #     predresult.iloc[id+2, 4+hla_pos_dic[hla]*5+2])
            # result['nM'].append(predresult.iloc[id+2, 4+hla_pos_dic[hla]*5+3])
            # result['Rank_BA'].append(
            #     predresult.iloc[id+2, 4+hla_pos_dic[hla]*5+4])
    result = pd.DataFrame(result)
    result.round({'Rank':4})
    result['pred_label'] = result['Rank'].map(lambda x: 1 if float(x)<10 else 0)
    result.to_csv(outfile, index=False)




