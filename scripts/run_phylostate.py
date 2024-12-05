import treeswift
import phylostate_libs as psl
from phylostate_libs.PhyloStateModel import PhyloStateModel
from phylostate_libs.ML_solver import ML_solver
from ast import literal_eval
import pandas as pd
import numpy as np
import timeit

class Logger(object):
    def __init__(self, output_prefix):
        self.terminal = sys.stdout
        self.log = open(output_prefix + ".log", "a")
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
    def flush(self):
        pass

import sys

output_prefix = "LAML_097_casshybrid_12h_celldivision_ancestral_cellstates"
sys.stdout = Logger(output_prefix)

#filename="/Users/gc3045/scmail_v1/sc-mail-experiments/Real_biodata/TLSCL/evaluation/LAML_097_casshybrid_extended_12h_celldivision.nwk"
filename="/Users/gc3045/scmail_v1/sc-mail-experiments/Real_biodata/TLSCL/evaluation/LAML_097_rep15_scaled_extended_12h_celldivision.nwk"
print("Running with...", filename)
t = treeswift.read_tree_newick(filename)

root_distr = np.array([1.0-5*psl.EPS, psl.EPS, psl.EPS, psl.EPS, psl.EPS, psl.EPS])
print("root_distr:", root_distr) 
state_set = {0,1,2,3,4,5}

metadata = "/Users/gc3045/scmail_v1/sc-mail-experiments/Real_biodata/TLSCL/evaluation/AM-DNA-097/AM-DNA-097_metadata.txt"
meta = pd.read_csv(metadata, sep="\t", index_col=0, header=0, converters={'cell_state': pd.eval})
cellstate_dict = meta['cell_state'].to_dict()
unknown_node_labels = set()
for node_label in cellstate_dict:
    if cellstate_dict[node_label][0] == 'Unknown':
        unknown_node_labels.add(node_label)
print("unknown node labels", unknown_node_labels) 

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

# get rid of the unknowns
print("Dropping", len(unknown_node_labels), "nodes with unknown cell state..." )
print(t.num_nodes(leaves=True, internal=False))

t_pruned = t.extract_tree(unknown_node_labels, True, suppress_unifurcations=False)
treeStr = t_pruned.newick()
print(t_pruned.num_nodes(leaves=True, internal=False))

lb2state_leaves = dict()
for node_label in cellstate_dict:
    if node_label not in unknown_node_labels:
        lb2state_leaves[node_label] = mapping[cellstate_dict[node_label][0]]




optimal_llh = -float("inf")
optimal_answer = None
for i in range(20): 
    start_time = timeit.default_timer()
    mySolver = ML_solver(treeStr,lb2state_leaves,state_set,root_distr)
    print("lb2idx:", mySolver.model.P_trans.lb2idx)
    print("idx2lbl:", mySolver.model.P_trans.idx2lb)
    print("root distr:", mySolver.model.root_distr)
    
    my_llh = mySolver.solve_EM()
    stop_time = timeit.default_timer()
    print("New llh",my_llh)
    print("Time:" + str(stop_time-start_time))
    print("P:", mySolver.model.P_trans.P_matrix)
    if my_llh > optimal_llh:
        optimal_llh = my_llh
        optimal_answer = mySolver

print('optimal',optimal_llh,np.round(optimal_answer.model.P_trans.P_matrix,5))      
optimal_tree = optimal_answer.model.tree 

outfile_mat = "LAML_097_casshybrid_12h_celldivision_ancestral_cellstates.mat"
with open(outfile_mat, "w+") as w:
    # m = np.round(optimal_answer.model.P_trans.P_matrix, 5)
    for row in optimal_answer.model.P_trans.P_matrix:
        row_str = '\t'.join([str(x) for x in row]) + "\n"
        w.write(row_str)

internal_node_idx = 0
for node in optimal_tree.traverse_preorder():
    if not node.label or node.label == "AutoLabel_1" or node.label == 'unseen': 
        node.label = f"InternalNode{internal_node_idx}"
        internal_node_idx += 1
    p = node.node_posterior
    #s = np.argmax(p)
    print(node.label, p)



outfilename = "LAML_097_casshybrid_12h_celldivision_ancestral_cellstates.labels"
internal_node_idx = 0
with open(outfilename, "w+") as w: 
    for node in optimal_tree.traverse_preorder():
        if not node.label or node.label == "AutoLabel_1" or node.label == 'unseen': 
            node.label = f"InternalNode{internal_node_idx}"
            internal_node_idx += 1
        p = node.node_posterior
        s = np.argmax(p)
        w.write('\t'.join([str(node.label), str(s)]) + "\n")

outfilename = "LAML_097_casshybrid_12h_celldivision_ancestral_cellstates.nwk"
#outfilename = "LAML_097_casshybrid_6h_celldivision_ancestral_cellstates_v3.nwk"
with open(outfilename, "w+") as w:
    w.write(optimal_tree.newick())
