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

for i in range(2, 3):
    print(i)
    # Read the file content into a DataFrame
    txt_path = os.path.join(file_path, ODSs[i] + ".txt")
    df = pd.read_csv(txt_path, delimiter=',', header=None, index_col=False)
    # Initialize OneHotEncoder with dtype=int to ensure numerical output
    encoder = OneHotEncoder(dtype=int)
    df_encoded = encoder.fit_transform(df)
    data = df_encoded.toarray()  # This is a NumPy array
    K = data_K_list[i]

    # ExShallow
    SHA = ShallowTree(K)
    SHA.fit_predict(data)
    print(SHA._max_depth())
    print(SHA.average_leaf_depth())
    # Tree plot saved to filename
    SHA.plot('So_SHA_Tree')

    # RandomThreshold #
    RDM = RandomThresholdTree(K)
    RDM.fit_predict(data)
    print(RDM._max_depth())
    print(RDM.average_leaf_depth())
    # Tree plot saved to filename
    RDM.plot('So_RDM_Tree')

    # IMM #
    IMM = Tree(k=K)
    IMM.fit_predict(data)
    print(IMM._max_depth())
    print(IMM.average_leaf_depth())
    # Tree plot saved to filename
    IMM.plot('So_IMM_Tree')