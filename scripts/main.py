import os, argparse
from preprocess import *
from result_process import *


def command_line():
    parser = argparse.ArgumentParser()
    parser.add_argument('-bl', dest='bl', 
                        help='Prediction for B cell linear epitopes.')
    parser.add_argument('-t1', dest='t1',
                        help='Prediction for T cell epitopes presented by MHC1 molecular.')
    parser.add_argument('-t2', dest='t2',
                        help='Prediction for T cell epitopes presented by MHC2 molecular.')
    parser.add_argument('-seq', dest='seqfile',
                        help='txt file contains one or more fasta sequences.')
    parser.add_argument('-hla', dest='hlastr',
                        help='HLA allele name, such as `HLA-A11:01,HLA-A24:02,` for t1 or `DRB1_01_01,DRB1_01_02` for t2')
    parser.add_argument('-hlafile', dest='hlafile',
                        help='`-hlastr` will be covered if this option is specified. This `hlafile` is necessary if t1 and t2 are both selected.')
    parser.add_argument('-o', dest='outdir', 
                        help='Directory to store predicted results. User must have write privilege. If omitted, the current directory will be applied.')

    return parser.parse_args()


def argCheck(args): 
    bl, t1, t2 = args.bl, args.t1, args.t2

    if args.bl != None and args.bl !='iedb7':
        KeyError(args.bl,'is not a validaed argument.')

    if args.t1 != None:
        if args.t1=='all':
            t1 = ['netmhcpan', 'mhcflurry', 'deephlapan'] 
        else:
            t1_list = args.t1.split(',')
            t1 = [i for i in t1_list if i in ['netmhcpan', 'mhcflurry', 'deephlapan']]

    if args.t2 != None:
        if args.t2=='all':
            t2 = ['netmhc2pan', 'mixmhc2pred'] 
        else:
            t2_list = args.t2.split(',')
            t2 = [i for i in t2_list if i in ['netmhc2pan', 'mixmhc2pred']]
                
    pred_options = [i for i,j in (zip(['bl','t1','t2'],[bl, t1, t2])) if j!=None]

    seq_file = args.seqfile
    task_name = seq_file.split('.')[0].split('/')[-1]

    if args.t1 or args.t2:
        input_hla = args.hlafile if args.hlafile else args.hlastr.split(',')
    else: input_hla = None

    outdir = os.path.dirname(os.path.abspath(__file__)) if args.outdir==None else args.outdir

    if t1 and t2:
        if not input_hla.endswith('.csv'):
            raise ValueError(
                'Class I and class II T cell epitopes prediction are selected, please input hla file in the specified format.')

    print(f'''
                                                        COVEP
    ----------------------------------------------------------------------------------------------------------------------
    Prediction for {','.join(pred_options)} will be performed.
    Task name: {task_name}
    Linear B epitope prediction: {bl}
    MHC I binding prediction: {t1}
    MHC II binding prediction: {t2}
    Seq File: {seq_file}
    HLA : {input_hla}
    Output Directory: {outdir}
    ----------------------------------------------------------------------------------------------------------------------

    ''')

    if not os.path.exists(outdir + '/cov_tools_predResult'):
        os.mkdir(outdir + '/cov_tools_predResult')
    if not os.path.exists(outdir + '/tmp'):
        os.mkdir(outdir + '/tmp')
    
    return bl, t1, t2, task_name, seq_file, input_hla, outdir


def prediction(bl, t1, t2, task_name, seq_file, input_hla, outdir):

    if bl == 'iedb7':
        from bl_functions import iedb
        process_iedb = iedb(task_name, seq_file, outdir)
        bl_merge_rank(process_iedb, task_name, outdir)

    if t1!=[] and t1!=None:
        from t1_functions import runMhcflurry, runNetmhcpan, runDeephlapan
        pepfile_for_hla1 = creatPepForHla1(seq_file, outdir)

        for i in t1:
            if i == 'mhcflurry':
                runMhcflurry(task_name, pepfile_for_hla1, input_hla, outdir)
            if i == 'netmhcpan':
                runNetmhcpan(task_name, pepfile_for_hla1, input_hla, outdir)
            if i == 'deephlapan':
                runDeephlapan(task_name, pepfile_for_hla1, input_hla, outdir)

        t1_merge_rank(t1, task_name, outdir)

    if t2!=[] and t2!=None:
        from t2_functions import runMixmhc2pred, runNetmhc2pan
        pepfile_for_hla2 = creatPepForHla2(seq_file, outdir)

        for i in t2:
            if i == 'netmhc2pan':
                runNetmhc2pan(task_name, pepfile_for_hla2, input_hla, outdir)
            if i == 'mixmhc2pred':
                runMixmhc2pred(task_name, pepfile_for_hla2, input_hla, outdir)
        t2_merge_rank(t2, task_name, outdir)

    os.system(f'rm -r {outdir}/tmp')
    print('-'*50, ' Prediction end! ', '-'*50)

    return


if __name__ == '__main__':
    
    args = command_line()
    bl, t1, t2, task_name, seq_file, input_hla, outdir = argCheck(args)
    prediction(bl, t1, t2, task_name, seq_file, input_hla, outdir)
    
    
