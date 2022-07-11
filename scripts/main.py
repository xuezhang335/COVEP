import os

from parse_args import *
from preprocess import *
from bl_prediction import *
from t1_prediction import *
from t2_prediction import *
from result_process import *


def argCheck(args):
    bl = args.bl
    t1 = args.t1
    t2 = args.t2
    pred_options = []
    if bl == True:
        pred_options.append('linear B cell epitopes')
    if t1 == True:
        pred_options.append('MHC class I epitopes')
    if t2 == True:
        pred_options.append('MHC class II epitopes')

    seq_file = args.seqfile
    task_name = seq_file.split('.')[0].split('/')[-1]

    if args.hlastr:
        input_hla = args.hlastr.split(',')
    if args.hlafile:
        input_hla = args.hlafile

    if t1 and t2:
        if not input_hla.endswith('.csv'):
            raise ValueError(
                'Class I and class II T cell epitopes prediction are selected, please input hla file in the specified format.')

    outdir = args.outdir
    if outdir == None:
        outdir = os.path.dirname(os.path.abspath(__file__))
        print(
            'Warning: Output directory not specified,the current directory will be applied.')

    print(f'''

    ----------------------------------------------------------------------------------------------------------------------
    Prediction for {','.join(pred_options)} will be performed.
    Input File: {seq_file}
    HLA : {input_hla}
    Output Directory: {outdir}
    ----------------------------------------------------------------------------------------------------------------------

    ''')
    
    return bl, t1, t2, task_name, seq_file, input_hla, outdir


def prediction(bl, t1, t2, task_name, seq_file, input_hla, outdir):
    if not os.path.exists(outdir + '/cov_tools_predResult'):
        os.mkdir(outdir + '/cov_tools_predResult')
    if not os.path.exists(outdir + '/tmp'):
        os.mkdir(outdir + '/tmp')

    if bl == True:
        process_iedb = iedb(task_name, seq_file, outdir)
        bl_merge_rank(process_iedb, task_name, outdir)

    if t1 == True:
        pepfile_for_hla1 = creatPepForHla1(seq_file, outdir)
        runMhcflurry(task_name, pepfile_for_hla1, input_hla, outdir)
        runNetmhcpan(task_name, pepfile_for_hla1, input_hla, outdir)
        runDeephlapan(task_name, pepfile_for_hla1, input_hla, outdir)
        t1_merge_rank(task_name, outdir)

    if t2 == True:
        pepfile_for_hla2 = creatPepForHla2(seq_file, outdir)
        runMixmhc2pred(task_name, pepfile_for_hla2, input_hla, outdir)
        runNetmhc2pan(task_name, pepfile_for_hla2, input_hla, outdir)
        t2_merge_rank(task_name, outdir)

    os.system(f'rm -r {outdir}/tmp')
    print('-'*50, ' Prediction end! ', '-'*50)

    return


if __name__ == '__main__':
    
    args = command_line()
    bl, t1, t2, task_name, seq_file, input_hla, outdir = argCheck(args)
    prediction(bl, t1, t2, task_name, seq_file, input_hla, outdir)
    
    
