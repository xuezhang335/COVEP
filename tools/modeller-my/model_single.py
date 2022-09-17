import sys
from modeller import *
from modeller.automodel import *
#from modeller import soap_protein_od

def model_single():
    sys.stdout = open('model_single.log', 'w')
    env = Environ()
    a = AutoModel(env, alnfile='TvLDH-1bdmA.ali',
                knowns='1bdmA', sequence='TvLDH',
                assess_methods=(assess.DOPE,
                                #soap_protein_od.Scorer(),
                                assess.GA341))
    a.starting_model = 1
    a.ending_model = 5
    a.make()
