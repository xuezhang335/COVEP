from cProfile import label
import os
import numpy as np
import pandas as pd

import tensorflow as tf
import tensorflow.python.util.deprecation as deprecation

from mhcflurry import Class1PresentationPredictor


os.environ['CUDA_VISIBLE_DEVICES'] = '0'
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'   # AVX2
tf.get_logger().setLevel('ERROR')   # warning
deprecation._PRINT_DEPRECATION_WARNINGS = False  # future discard

# print("="*30,"GPU USED","="*30,"tf version:",tf.__version__,"use GPU:",tf.test.is_gpu_available())

cov_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def runDeephlapan(task_name, pep1_file, input_hla, outdir):
    print('='*50, ' RUNNING DeepHLApan... ', '='*50)
    df = pd.read_csv(pep1_file)
    pep_list = list(df.peptide)
    id_list = list(df.seq_id)
    pos_list = list(df.pos)

    if type(input_hla)==type(''):
        df = pd.read_csv(input_hla)
        hla1 = list(df.hla1.dropna())
    else:
        hla1 = input_hla

    infile = outdir + '/tmp/' + '_'.join(os.path.basename(pep1_file).split('_')[:-1]) + '_deephlapan.csv'

    pep_for_deephlapan = [i for i in pep_list for j in range(len(hla1))]
    id_list = [i for i in id_list for j in range(len(hla1))]
    pos_list = [i for i in pos_list for j in range(len(hla1))]
    hla1_for_deephlapan = hla1*len(pep_list)

    assert len(pep_for_deephlapan) == len(hla1_for_deephlapan)
    pd.DataFrame({'seq_id': id_list,
                  'HLA': hla1_for_deephlapan,
                  'position': pos_list,
                  'peptide': pep_for_deephlapan
                  }).to_csv(infile, index=False)
    outdir_tmp = outdir+'/tmp'
    # os.system('source activate /home/biopharm/anaconda3/envs/py27_1')
    os.system(f'python {cov_dir}/tools/deephlapan/deephlapan/deephlapan.py -F {infile} -O {outdir_tmp}') 
    # os.system('source activate /share/anaconda3/envs/cov')
    
    result = pd.read_csv(outdir+'/tmp/'+task_name+'_pep_for_deephlapan_predicted_result.csv')
    result['binding_label'] = result['Bindingscore'].apply(lambda x:1 if x>=0.5 else 0)
    result['immuno_label'] = result['Immunogenicscore'].apply(lambda x:1 if x>=0.5 else 0)

    def aff2label(affscore, immunoscore):
        if affscore>0.5 and immunoscore>0.5:
            pred_label = 1
        else:
            pred_label = 0
        return pred_label
    # result['pred_label'] = result.apply(lambda x:aff2label(x['Bindingscore'],x['Immunogenicscore']),axis=1)
    
    outfile = outdir + '/cov_tools_predResult/t1_' + task_name + '_deephlapan.csv'
    result.to_csv(outfile, sep=',', index=False)

    print('DeepHLApan Prediction end.')

    return


def runMhcflurry(task_name, pep1_file, input_hla, outdir):
    print('='*50, 'RUNNING MHCflurry...', '='*50)
    df = pd.read_csv(pep1_file)
    pep_list = list(df.peptide)
    id_list = list(df.seq_id)
    pos_list = list(df.pos)

    hla_dict = {}
    if type(input_hla)==type(''):
        df = pd.read_csv(input_hla)
        hla1 = list(df.hla1.dropna())
    else:
        hla1 = input_hla

    for id, hla in enumerate(hla1):
        key = 'sample' + str(id)
        hla_dict[key] = [hla.split("-")[-1].replace(':', '')]

    predictor = Class1PresentationPredictor.load()
    result = predictor.predict(peptides=pep_list, alleles=hla_dict, verbose=0)

    result['sample_name'] = result['sample_name'].apply(lambda x: int(x[6:]))
    result = result.sort_values(by=['peptide_num','sample_name'])

    id_list = [id for id in id_list for i in range(len(hla_dict.keys()))]
    pos_list = [pos for pos in pos_list for i in range(len(hla_dict.keys()))]
    assert len(result) == len(id_list) == len(pos_list)
    result.insert(0, 'Seq_id', id_list)
    result.insert(1, 'best_allele', result.pop('best_allele'))
    result.insert(2, 'Position', pos_list)
    result.rename(columns={'best_allele': 'Allele','peptide':'Peptide'}, inplace=True)
    result['Allele'] = result['Allele'].apply(lambda x: 'HLA-'+x[0:3]+':'+x[3:])
    result.round({'presentation_percentile':4})
    result['pred_label'] = result['presentation_percentile'].map(lambda x: 1 if x<2 else 0)
    outfile = outdir + '/cov_tools_predResult/t1_' + task_name + '_mhcflurry.csv'
    result.iloc[:,[0,1,2,3,4,-2,-1]].to_csv(outfile, sep=',', index=False)
    print('MHCflurry Prediction end.')

    return


def runNetmhcpan(task_name, pep1_file, input_hla, outdir):
    # import signal
    # signal.signal(signal.SIGPIPE, signal.SIG_IGN)  # 忽略SIGPIPE信号

    print('='*50, ' RUNNING NetMHCpan... ', '='*50)
    print('In progress.....')
    df = pd.read_csv(pep1_file)
    pep_list = list(df.peptide)
    id_list = list(df.seq_id)
    pos_list = list(df.pos)

    if type(input_hla)==type(''):
        df = pd.read_csv(input_hla)
        hla1 = list(df.hla1.dropna())
    else:
        hla1 = input_hla

    str_hla1 = ','.join(hla1)

    infile = outdir + '/tmp/' + '_'.join(os.path.basename(pep1_file).split('_')[:-1]) + '_netmhcpan.txt'
    outfile_tmp = outdir + '/tmp/netmhcpan_out.xls'
    terminalout = outdir + '/tmp/netmhcpan_out.print'

    with open(infile, 'w')as f:
        for p in pep_list:
            f.write(p + '\n')

    os.system(
        f'{cov_dir}/tools/netMHCpan-4.1/netMHCpan -p {infile} -BA -xls -a {str_hla1} -xlsfile {outfile_tmp} > {terminalout}')  # 

    processNetmhcpanOutput(
        task_name, pep_list, id_list, pos_list, hla1, outdir)

    print('NetMHCpan Prediction end.')
    return


def processNetmhcpanOutput(task_name, pep_list, id_list, pos_list, hla1, outdir):
    hla1_sort = sorted(hla1)
    hla_pos_dic = {}
    for hla in hla1:
        pos = hla1_sort.index(hla)
        hla_pos_dic[hla] = pos

    outfile_tmp = outdir + '/tmp/netmhcpan_out.xls'
    predresult = pd.read_csv(outfile_tmp, sep='\t')  # low_memory=False
    outfile = outdir + '/cov_tools_predResult/t1_' + task_name + '_netmhcpan.csv'
    result = {'Seq_id': [], 'Allele': [], 'Position': [], 'Peptide': [], 'BArank': []}

    for id, pep in enumerate(pep_list):
        for hla in hla_pos_dic:
            result['Seq_id'].append(id_list[id])
            result['Allele'].append(hla)
            result['Position'].append(pos_list[id])
            result['Peptide'].append(pep)
            # result['core'].append(predresult.iloc[id+1, 3+hla_pos_dic[hla]*6])
            # result['icore'].append(
            #     predresult.iloc[id+1, 3+hla_pos_dic[hla]*6+1])
            # result['EL-score'].append(predresult.iloc[id+1,
            #                           3+hla_pos_dic[hla]*6+2])
            # result['EL_Rank'].append(
            #     predresult.iloc[id+1, 3+hla_pos_dic[hla]*6+3])
            # result['BAscore'].append(predresult.iloc[id+1,
            #                           3+hla_pos_dic[hla]*6+4])
            result['BArank'].append(
                predresult.iloc[id+1, 3+hla_pos_dic[hla]*6+5])
    result = pd.DataFrame(result)
    result.round({'BArank':4})
    
    # result['pred_label']=result.apply(lambda x:function(x['BAscore'],x['BArank']),axis=1)
    result['pred_label'] = result['BArank'].apply(lambda x: 1 if float(x)<2 else 0)

    result.to_csv(outfile, index=False)



