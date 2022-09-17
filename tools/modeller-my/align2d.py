from modeller import *
import sys

def align2d(target_seq, best_model):
    sys.stdout = open('align2d.log', 'w')
    pdb, chain = best_model[0,4], best_model[-1]
    env = Environ()
    aln = Alignment(env)
    mdl = Model(env, file=pdb, model_segment=('FIRST:'+chain,'LAST'+chain))
    aln.append_model(mdl, align_codes=pdb+chain, atom_files=pdb+'.pdb')
    aln.append(file=target_seq+'.ali', align_codes=target_seq)
    aln.align2d(max_gap_length=50)
    aln.write(file=f'{target_seq}-{best_model}.ali', alignment_format='PIR')
    aln.write(file=f'{target_seq}-{best_model}.pap', alignment_format='PAP')
