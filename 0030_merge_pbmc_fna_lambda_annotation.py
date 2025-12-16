#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 18:15:13 2025

@author: wuw5
"""

import os
import pandas as pd
import numpy as np
import glob

####=============
###read within_fna_annotated file
within_fna_df = pd.read_csv(r"/~/20251103_FNA_PBMC/Lambda/1115Lambda_chain_combined_df_max.csv")
fna_matched_pbmc_df = pd.read_csv("/~/20251103_FNA_PBMC/Lambda/44_lambda_chains_matched_to_PBMC_flags.csv")
fna_matched_pbmc_df["PBMC_Match"] = "MATCHED TO PBMC"

## merge the pbmc flags to fna dataset
# Step 1: Take only the relevant columns
cols_to_keep = [
    'first_sequence_id',
    'FNA_shared_PBMC_cross_donor',
    'FNA_shared_PBMC_cross_epitope',
    'PBMC_Match'
]
pbmc_subset = fna_matched_pbmc_df[cols_to_keep].drop_duplicates(subset='first_sequence_id')

# Step 2: Merge with within_fna_df
within_fna_df = within_fna_df.merge(
    pbmc_subset,
    how='left',
    left_on='1st_sequence_id',
    right_on='first_sequence_id'
)

# Step 3: Optional — drop the redundant merge key
within_fna_df = within_fna_df.drop(columns='first_sequence_id')

# Step 4: Preview result
# print("✅ Merge complete!")
# print(f"Final shape: {within_fna_df.shape}")
# print("New columns added:", [c for c in ['FNA_shared_PBMC_cross_donor', 'FNA_shared_PBMC_cross_epitope', 'PBMC_Match'] if c in within_fna_df.columns])


within_fna_df["Keep/Discard"] = np.where(
    within_fna_df["cross_donor_fna_Bcell"].astype(str).str.contains("yes", case=False, na=False)
    | within_fna_df["cross_epitope_fna_Bcell"].astype(str).str.contains("yes", case=False, na=False),
    "Discard",
    "Keep"
)


# Identify groups (wells) that contain at least one "yes"
group_has_yes = within_fna_df.groupby("1st_sequence_or_well_id")["unique_within_fna_Bcell"].transform(
    lambda x: (x == "yes").any()
)

# Label non-"yes" rows within those groups as "Discard"
within_fna_df["Keep/Discard_multi_within_well"] = np.where(
    (group_has_yes) & (within_fna_df["unique_within_fna_Bcell"] != "yes"),
    "Discard",
    np.nan  # Leave the unique/yes rows blank
)


within_fna_df.to_csv("/~/20251103_FNA_PBMC/Lambda/0030_1115_lambda_chains_with_fna_pbmc_flags.csv", index=False)

