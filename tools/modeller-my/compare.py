from modeller import *
import sys

def compare(template):
    sys.stdout = open('compare.log', 'w')
    env = Environ()
    aln = Alignment(env)
    for (pdb, chain) in template:
        m = Model(env, file=pdb, model_segment=('FIRST:'+chain, 'LAST:'+chain))
        aln.append_model(m, atom_files=pdb, align_codes=pdb+chain)
    aln.malign()
    aln.malign3d()
    aln.compare_structures()
    aln.id_table(matrix_file='family.mat')
    env.dendrogram(matrix_file='family.mat', cluster_cut=-1.0)
