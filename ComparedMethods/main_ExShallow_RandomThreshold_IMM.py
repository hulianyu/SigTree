from ShallowTree.ShallowTree import ShallowTree
from ShallowTree.RandomThresholdTree import RandomThresholdTree
from ExKMC.Tree import Tree
from sklearn.preprocessing import OneHotEncoder
import os
import pandas as pd
import numpy as np
import time
from scipy.io import savemat

# Set the path to the target folder (can be relative or absolute)
current_directory = os.path.dirname(__file__)
file_path = os.path.join(current_directory, 'ODS')
data_obj_list = ["lenses", "lung_cancer", "soybean_small", "zoo", "dna_promoter", "hayes_roth",
                 "lymphography", "heart_disease", "solar_flare", "primary_tumor", "dermatology", "house_votes",
                 "balance_scale", "credit_approval", "breast_cancer_wisconsin", "mammographic_mass", "tic_tac_toe", "car"]
data_K_list = [3, 3, 4, 7, 2, 3, 4, 5, 6, 21, 6, 2, 3, 2, 2, 2, 2, 4]
ODSs = ["ODS_" + item for item in data_obj_list]

SHA_pi = []
RDM_pi = []
IMM_pi = []
Execution_times = np.zeros((18, 3))
MaxDepth = np.zeros((18, 3))
AvgLeafDepth = np.zeros((18, 3))
NumLeaf = np.zeros((18, 3))

num_runs = 50

for i in range(0, 18):
    print(i)
    # Read the file content into a DataFrame
    txt_path = os.path.join(file_path, ODSs[i] + ".txt")
    df = pd.read_csv(txt_path, delimiter=',', header=None, index_col=False)
    # Initialize OneHotEncoder with dtype=int to ensure numerical output
    encoder = OneHotEncoder(dtype=int)
    df_encoded = encoder.fit_transform(df)
    data = df_encoded.toarray()  # This is a NumPy array
    K = data_K_list[i]

    for run in range(num_runs):
        # ExShallow
        start_SHA = time.time()
        SHA = ShallowTree(K)
        # Construct the tree, and return cluster labels
        SHA_pi.append(SHA.fit_predict(data))
        end_SHA = time.time()
        Execution_times[i, 0] = Execution_times[i, 0] + (end_SHA - start_SHA)
        # Max Depth, Average Leaf Depth
        MaxDepth[i, 0] = MaxDepth[i, 0] + SHA._max_depth()
        AvgLeafDepth[i, 0] = AvgLeafDepth[i, 0] + SHA.average_leaf_depth()
        NumLeaf[i, 0] = NumLeaf[i, 0] + SHA.count_leaves()
        # Tree plot saved to filename
        # SHA.plot('eg_SHA_pi')

        # RandomThreshold #
        start_RDM = time.time()
        RDM = RandomThresholdTree(K)
        # Construct the tree, and return cluster labels
        RDM_pi.append(RDM.fit_predict(data))
        end_RDM = time.time()
        Execution_times[i, 1] = Execution_times[i, 1] + (end_RDM - start_RDM)
        # Max Depth, Average Leaf Depth
        MaxDepth[i, 1] = MaxDepth[i, 1] + RDM._max_depth()
        AvgLeafDepth[i, 1] = AvgLeafDepth[i, 1] + RDM.average_leaf_depth()
        NumLeaf[i, 1] = NumLeaf[i, 1] + RDM.count_leaves()
        # Tree plot saved to filename
        # RDM.plot('eg_RDM_pi')

        # IMM #
        start_IMM = time.time()
        IMM = Tree(k=K)
        # Construct the tree, and return cluster labels
        IMM_pi.append(IMM.fit_predict(data))
        end_IMM = time.time()
        Execution_times[i, 2] = Execution_times[i, 2] + (end_IMM - start_IMM)
        # Max Depth, Average Leaf Depth
        MaxDepth[i, 2] = MaxDepth[i, 2] + IMM._max_depth()
        AvgLeafDepth[i, 2] = AvgLeafDepth[i, 2] + IMM.average_leaf_depth()
        NumLeaf[i, 2] = NumLeaf[i, 2] + IMM.count_leaves()
        # Tree plot saved to filename
        # IMM.plot('eg_IMM_pi')

MaxDepth = MaxDepth/num_runs
Execution_times = Execution_times/num_runs
AvgLeafDepth = AvgLeafDepth/num_runs
NumLeaf = NumLeaf/num_runs

savemat('SHA_pi.mat', {'SHA_pi': SHA_pi})
savemat('RDM_pi.mat', {'RDM_pi': RDM_pi})
savemat('IMM_pi.mat', {'IMM_pi': IMM_pi})
savemat('MaxDepth.mat', {'MaxDepth': MaxDepth})
savemat('Execution_times.mat', {'Execution_times': Execution_times})
savemat('AvgLeafDepth.mat', {'AvgLeafDepth': AvgLeafDepth})
savemat('NumLeaf.mat', {'NumLeaf': NumLeaf})