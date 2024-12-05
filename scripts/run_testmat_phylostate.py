from laml_libs.sim_lib import get_balanced_tree
from phylostate_libs.PhyloStateModel import PhyloStateModel
from phylostate_libs.ML_solver import ML_solver
import numpy as np
import phylostate_libs as psl
import timeit
import random
from treeswift import *
import statistics 

# reasonable model
matrix = np.loadtxt('optimal_llh.mat', usecols=range(6))
P = dict()
with open('optimal_llh.mat', 'r') as f:
    for i, line in enumerate(f):
        tokens = line.split()
        for j, val in enumerate(tokens):
            if i not in P:
                P[i] = dict()
            P[i][j] = float(val)
print(P)
filename = "/Users/gc3045/scmail_v1/sc-mail-experiments/Real_biodata/TLSCL/evaluation/LAML_097_rep15_scaled_extended_12h_celldivision.nwk"
T = read_tree_newick(filename)

root_distr = np.array([1.0-5*psl.EPS, psl.EPS, psl.EPS, psl.EPS, psl.EPS, psl.EPS])

total_num = 0
state_set = {0, 1, 2, 3, 4, 5}
mapping_counts = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0}

error_proportions = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
mapping_proportions = {0: 0, 1: 0.05, 2: 0.68, 3: 0.06, 4: 0.28, 5: 0.56}
for j in range(10):
    trueModel = PhyloStateModel(T.newick(),P,root_distr)
    lb2state_all,lb2state_leaves = trueModel.simulate()
    true_llh_all = trueModel.compute_llh_all_observed(lb2state_all)
    true_llh_leaves = trueModel.compute_in_llh(lb2state_leaves)
    print('true_llh',true_llh_leaves)

    # compute the number of nodes at the leaves with the labels

    for node_label in lb2state_leaves:
        total_num += 1
        mapping_counts[lb2state_leaves[node_label]] += 1

    for val in mapping_counts:
        mapping_counts[val] = mapping_counts[val]/total_num
    print(mapping_counts)

    for val in mapping_counts:
        error_proportions[val].append(mapping_proportions[val] - mapping_counts[val])

print([statistics.mean(error_proportions[val]) for val in error_proportions])

"""
    optimal_llh = -float("inf")
    optimal_answer = None
    for i in range(10):    
        start_time = timeit.default_timer()
        mySolver = ML_solver(treeStr,lb2state_leaves,state_set,root_distr)
        my_llh = mySolver.solve_EM()
        stop_time = timeit.default_timer()
        print("New llh",my_llh)
        print("Time:" + str(stop_time-start_time))
        if my_llh > optimal_llh:
            optimal_llh = my_llh
            optimal_answer = mySolver
    print('optimal',optimal_llh,np.round(optimal_answer.model.P_trans.P_matrix,3))      
    optimal_tree = optimal_answer.model.tree 

"""
#outfilename = "test_tmp"
#with open(outfilename, "w+") as w:
#    w.write(optimal_tree.newick())

#for node in optimal_tree.traverse_preorder():
#    p = node.node_posterior
#    print(p)
#    s = np.argmax(p)
#    print(node.label,s,lb2state_all[node.label]) 
"""
outfile_mat = "simulation_ancestral_cellstates.mat"
with open(outfile_mat, "w+") as w:
    m = optimal_answer.model.P_trans.P_matrix
    for row in m:
        row_str = '\t'.join([str(x) for x in row]) + "\n"
        w.write(row_str)

outfilename = "simulation_ancestral_cellstates.labels"
internal_node_idx = 0
with open(outfilename, "w+") as w:
    for node in optimal_tree.traverse_preorder():
        if not node.label or node.label == "AutoLabel_1" or node.label == 'unseen':
            node.label = f"InternalNode{internal_node_idx}"
            internal_node_idx += 1
        p = node.node_posterior
        s = np.argmax(p)
        w.write('\t'.join([str(node.label), str(s)]) + "\n")

outfilename = "simulation_ancestral_cellstates.nwk"
with open(outfilename, "w+") as w:
    w.write(optimal_tree.newick())

"""
