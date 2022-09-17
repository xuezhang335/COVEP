import enum
import os, sys
import pandas as pd
import regex as re
from build_profile import build_profile
from compare import compare
from align2d import align2d
from model_single import model_single
from collections import defaultdict


def get_template():
    compiler = re.compile('[\S]+')
    lines = open('build_profile.prf').readlines()
    df = pd.DataFrame(columns=['index', 'pdb_chain', 2, 3, 4, 5, 6, 7, 8, 9, 'identity','evalue','sequence'])
    # 2：模板的PDB入口代码 。11：相似度百分比，大于25%即有可能成为有效模板。 12：错误率，越低越好，最好等于0
    for line in lines[6+1:]:
        line.strip()
        row = compiler.findall(line)
        df.loc[len(df.index)] = row
    df.sort_values(by=['evalue'], inplace=True)
    df_template = df[df['evalue']=='0.0']
    df_template['pdb'] = df_template['pdb_chain'].map(lambda x: x[:4])
    df_template['chain'] = df_template['pdb_chain'].map(lambda x: x[-1])
    template_dict = defaultdict(dict)

    for pdb_chain in df_template['pdb_chain'].values:
        if len(pdb_chain)==5:
            pdb, chain = pdb_chain[:4], pdb_chain[-1]
        elif len(pdb_chain)==4:
            pdb, chain = pdb_chain, ''
            
        template_dict[pdb_chain]['pdb'] = pdb
        template_dict[pdb_chain]['chain'] = chain
        template_dict[pdb_chain]['identity'] = df_template[df_template['pdb_chain']==pdb_chain]['identity'].values[0]
        template_dict[pdb_chain]['evalue'] = df_template[df_template['pdb_chain']==pdb_chain]['evalue'].values[0]

        os.system(f'wget -O {pdb}.pdb "http://www.rcsb.org/pdb/files/{pdb}.pdb"')
    
    return template_dict


# print('='*50, sys.argv[1], '='*50)
seq_file = 'seq.ali'
build_profile(seq_file=seq_file)
template_dict =get_template()
template = [(template_dict[pdb_chain]['pdb'], template_dict[pdb_chain]['chain']) for pdb_chain in template_dict]
compare(template=tuple(template))

# 根据compare.log中的信息，选择最好的模版


