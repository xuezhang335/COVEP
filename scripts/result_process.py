import os
import numpy as np
import pandas as pd


def bl_merge_rank(iedb_process,task_name, outdir):
    methods = []
    for method in iedb_process.keys():
        if iedb_process[method] == True:
            methods.append(method)
        else:
            os.system(f'rm {outdir}/cov_tools_predResult/bl_{task_name}_{method}.csv')

    if method != []:
        merge_result = pd.DataFrame()
        for i, method in enumerate(methods):
            outfile = f'{outdir}/cov_tools_predResult/bl_{task_name}_{method}.csv'
            df = pd.read_csv(outfile)
            if i==0:
                merge_result = df.iloc[:,0:3]
            merge_result[method] = df['pred_label']
        
        merge_result['Score'] = merge_result.iloc[:,3:].apply(np.sum, axis=1)
        merge_file = outdir+f'/bl_{task_name}_predResult.csv'
        merge_result.to_csv(merge_file, index=False)
    else:
        print('''Warning: Due to network reasons, the prediction result of B cell epitope is empty,
                 so that this method prediction result analysis will not be carried out.''')
    return


def t1_merge_rank(task_name, outdir):
    methods = ['mhcflurry', 'netmhcpan', 'deephlapan']
    for i, method in enumerate(methods):
        outfile = f'{outdir}/cov_tools_predResult/t1_{task_name}_{method}.csv'
        df = pd.read_csv(outfile)

        if i==0:
            merge_result = df.iloc[:,0:4]

        if method=='deephlapan':
            merge_result[method+'_BS'] = df['binding_label']
            merge_result[method+'_IS'] = df['immuno_label']
        else:
            merge_result[method] = df['pred_label']

    merge_result['Score'] = merge_result.iloc[:,4:].apply(np.sum, axis=1)
    merge_file = outdir+f'/t1_{task_name}_predResult.csv'
    merge_result.to_csv(merge_file, index=False)

    return


def t2_merge_rank(task_name, outdir):
    methods = ['mixmhc2pred', 'netmhc2pan']
    for i, method in enumerate(methods):
        outfile = f'{outdir}/cov_tools_predResult/t2_{task_name}_{method}.csv'
        df = pd.read_csv(outfile)
        if i==0:
            merge_result = df.iloc[:,0:4]
        merge_result[method] = df['pred_label']

    merge_result['Score'] = merge_result.iloc[:,4:].apply(np.sum, axis=1)
    merge_file = outdir+f'/t2_{task_name}_predResult.csv'
    merge_result.to_csv(merge_file, index=False)

    return