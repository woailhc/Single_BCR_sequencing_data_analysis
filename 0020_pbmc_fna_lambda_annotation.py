#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 17:43:43 2025

@author: wuw5
"""
import os
import pandas as pd
import numpy as np
import glob
# Define folder path
folder_path = "/~/20251028_FNA_PBMC/Lambda/LeanFiles_Per_Donor"

# Find all CSV files in the folder
csv_files = glob.glob(os.path.join(folder_path, "*.csv"))

# Read and concatenate them
df = pd.concat((pd.read_csv(f) for f in csv_files), ignore_index=True)

df.columns = df.columns.str.replace('\n', '_', regex=False)

# Optionally, make names more uniform (e.g., lowercase)
df.columns = df.columns.str.lower()

# Check the cleaned-up list
df.columns.tolist()

df["first_real_subjectid"] = df["first_real_subjectid"].astype(int)
df["second_real_subjectid"] = df["second_real_subjectid"].astype("Int64")



os.chdir(r"/Users/wuw5/Desktop/notes/20251028_FNA_PBMC/")
from classify_sequences_labmda_99 import classify_all_sequences

# Apply to your DataFrame
df = classify_all_sequences(
    df,
    groupby_col='first_sequence_id',
    identity_col='aa_local_percentage_identity'
)
grouped = (
    df[df['FNA_shared_PBMC_cross_donor'] != ""]
    .groupby('first_sequence_id')['FNA_shared_PBMC_cross_donor']
    .nunique()
)
ids_with_multiple = grouped[grouped > 1].index.tolist()

df[['FNA_shared_PBMC_cross_donor','FNA_shared_PBMC_cross_epitope']].value_counts()

# FNA_shared_PBMC_cross_donor  FNA_shared_PBMC_cross_epitope
#                                                               1056632
# yes and no                   yes                                  406
# no                                                                184
# yes                          yes                                   29
# no                           yes and no                            28
#                              yes                                   10
#                              no                                     3
# yes and no                                                          2
filtered_df = df[
    (df["FNA_shared_PBMC_cross_donor"].isin(["yes", "no", "yes and no"])) |
    (df["FNA_shared_PBMC_cross_epitope"].isin(["yes", "no", "yes and no"]))
]
filtered_df['first_sequence_id'].nunique()
#Out[42]: 44
filtered_df.to_csv("/~/20251103_FNA_PBMC/Lambda/44_lambda_chains_matched_to_PBMC_flags.csv", index=False)





