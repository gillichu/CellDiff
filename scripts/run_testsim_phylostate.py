from laml_libs.sim_lib import get_balanced_tree
from phylostate_libs.PhyloStateModel import PhyloStateModel
from phylostate_libs.ML_solver import ML_solver
import numpy as np
import phylostate_libs as psl
import timeit
import random
from treeswift import *

# true model
P = {0:{0:0.73, 1:0.0127, 2:0.243, 3:0.0143, 4:0.0,   5:0.0}, # pluripotent progenitor
     1:{0:0.0, 1:1.0,  2:0.0,   3:0.0,   4:0.0,   5:0.0}, # sink state, pcglc
     2:{0:0.0, 1:0.0,  2:0.753,   3:0.0,   4:0.083, 5:0.17},
     3:{0:0.0, 1:0.0,  2:0.0,   3:1.0,   4:0.0,   5:0.0}, # sink state, endoderm / endothelial
     4:{0:0.0, 1:0.0,  2:0.0,   3:0.0,   4:1.0,   5:0.0}, # sink state
     5:{0:0.0, 1:0.0,  2:0.0,   3:0.0,   4:0.0,   5:1.0}}

#P1 = {0:{0:0.5,1:0.5},
#     1:{1:0.3,2:0.7},
#     2:{2:1}}

treeStr = get_balanced_tree(13,1)
#treeStr = get_balanced_tree(11,1)
T = read_tree_newick(treeStr)
L = [node.label for node in T.traverse_leaves()]
L_small = random.sample(L,1800)
#L_small = random.sample(L,18000)
T_pruned = T.extract_tree_with(L_small,suppress_unifurcations=False)
treeStr = T_pruned.newick()
root_distr = np.array([1.0-5*psl.EPS, psl.EPS, psl.EPS, psl.EPS, psl.EPS, psl.EPS])

#print(treeStr)
#treeStr = "[&R] ((((((((a:2)unseen:2)unseen:2)unseen:2)unseen:2,((((b:2)unseen:2)unseen:2)unseen:2)unseen:2)ab:1):1.5)unseen:1.5,(((((((c:2)unseen:2)unseen:2)unseen:2)unseen:2,((((d:2)unseen:2)unseen:2)unseen:2)unseen:2)cd:1):2)unseen:2):2;"
#root_distr = np.array([1-2*psl.EPS,psl.EPS,psl.EPS])
#state_set = {0,1,2}

total_num = 0
state_set = {0, 1, 2, 3, 4, 5}
mapping_counts = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0}
for j in range(1):
    trueModel = PhyloStateModel(treeStr,P,root_distr)
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


#outfilename = "test_tmp"
#with open(outfilename, "w+") as w:
#    w.write(optimal_tree.newick())

#for node in optimal_tree.traverse_preorder():
#    p = node.node_posterior
#    print(p)
#    s = np.argmax(p)
#    print(node.label,s,lb2state_all[node.label]) 

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


