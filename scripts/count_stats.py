import treeswift
import phylostate_libs as psl
from phylostate_libs.PhyloStateModel import PhyloStateModel
from phylostate_libs.ML_solver import ML_solver
from ast import literal_eval
import pandas as pd
import numpy as np
import timeit

metadata = "/Users/gc3045/scmail_v1/sc-mail-experiments/Real_biodata/TLSCL/evaluation/AM-DNA-097/AM-DNA-097_metadata.txt"
meta = pd.read_csv(metadata, sep="\t", index_col=0, header=0, converters={'cell_state': pd.eval})
cellstate_dict = meta['cell_state'].to_dict()



mapping = {'Pluripotent-progenitor':0,
           'PCGLC': 1,
           #'Germlayer-progenitor': 2,
           'NMPs': 2,
           'Endoderm': 3,
           'Endothelial': 3,
           'NeuralTube1': 4,
           'NeuralTube2': 4,
           'Somite': 5,
           'Somite-1': 5,
           'Somite0': 5,
           'SomiteDermo': 5,
           'SomiteSclero': 5,
           'aPSM': 5,
           'pPSM': 5,
           }


mapping_counts = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0}
# count how many of each label is observed 
unknown = 0
for node_label in cellstate_dict:

    if cellstate_dict[node_label][0] != 'Unknown':

        mapping_counts[mapping[cellstate_dict[node_label][0]]] += 1
    else:
        unknown += 1


print(mapping_counts)
print("unknown:", unknown)
