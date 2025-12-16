#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 16:55:21 2025

@author: wuw5
"""
import os
import pandas as pd
import numpy as np
# Correct usage
df = pd.read_excel("/~/20251023_nt_alignment_FNA/Lambda/000460_FNA_LAMBDA_Alignments_Permutations_with_MetaData_LeanFile_DELIVERY_9.xlsx")
# Remove newline characters
df.columns = df.columns.str.replace('\n', '_', regex=False)

# Optionally, make names more uniform (e.g., lowercase)
df.columns = df.columns.str.lower()

# Check the cleaned-up list
df.columns.tolist()


# Step 1: Shallow copy and swap "1st" and "2nd" in column names

df_swapped = df.copy(deep=False)

# Temporarily replace "1st" with a placeholder
df_swapped.columns = df_swapped.columns.str.replace('1st', 'temp', regex=False)

# Now swap "2nd" → "1st"
df_swapped.columns = df_swapped.columns.str.replace('2nd', '1st', regex=False)

# Finally replace "temp" → "2nd"
df_swapped.columns = df_swapped.columns.str.replace('temp', '2nd', regex=False)


#df_swapped_path_csv = os.path.join(out_dir, "swapped_000460_FNA_LAMBDA_Alignments_Permutations_with_MetaData_LeanFile.csv")
#df_swapped.to_csv(df_swapped_path_csv, index=False)


# Step 2: Concatenate (rbind equivalent)
combined_df = pd.concat([df, df_swapped], axis=0, ignore_index=True)


# Read non_compared chains with unique cdr3 length unique into a DataFrame
non_compared_df = pd.read_excel("/~/20251031_FNA_PBMC_version4/FNA_vs_FNA__UniqCDR3LenSingletons.xlsx")
non_compared_df.columns = non_compared_df.columns.str.replace('\n', '_', regex=False)
non_compared_df.columns = non_compared_df.columns.str.lower()
non_compared_df.columns = ['1st_' + col for col in non_compared_df.columns]
common_cols = non_compared_df.columns.intersection(combined_df.columns)
print(common_cols)
# Filter rows from non_compared_df where '1st_sequence_id' contains 'lambda'
filtered_non_compared = non_compared_df[
    non_compared_df['1st_sequence_id'].astype(str).str.contains('lambda', case=False, na=False)
]
filtered_non_compared = filtered_non_compared[
    [col for col in filtered_non_compared.columns if col in combined_df.columns]
]
concated_df = pd.concat(
    [combined_df, filtered_non_compared],
    ignore_index=True,
    sort=False
)

# Identify rows where 'nt_global_% of _identity' is missing or empty
missing_identity_mask = (
    concated_df['nt_global_% of _identity'].isna() |
    (concated_df['nt_global_% of _identity'].astype(str).str.strip() == "")
)

# Define (target, source) column pairs
col_pairs = [
    ('2nd_cell_type', '1st_cell_type'),
    ('2nd_flowindex_epitope', '1st_flowindex_epitope'),
    ('2nd_real_subjectid', '1st_real_subjectid')
]

# Apply only to rows where identity is missing
for target, source in col_pairs:
    concated_df.loc[missing_identity_mask, target] = concated_df.loc[missing_identity_mask, source]

# Count how many rows were affected
num_rows_affected = missing_identity_mask.sum()

print(f"✅ Filled 2nd_* columns for {num_rows_affected} row(s) where 'nt_global_% of _identity' was missing.")


combined_df = concated_df.copy()
# Check the first few rows
print(non_compared_df.head())

os.chdir(r"/~/20251029_FNA_PBMC/")
from module001_bt_cell_match import add_bt_cell_match
combined_df_original = combined_df.copy()
combined_df = add_bt_cell_match(
    combined_df,
    identity_col='nt_global_% of _identity',
    global_col='nt_global_alignment_score',
    cell_type_col_1='1st_cell_type',
    cell_type_col_2='2nd_cell_type',
    global_threshold=50
)

────────── B-T CELL MATCH SUMMARY ──────────
B-T CELL MATCH
NA                                 757857
WARNING: B-T CELL 100% NT MATCH     11588
Name: count, dtype: int64
────────────────────────────────────────────

# 2. add a "cross_cell_type", for combined_df['1st_cell_type'] == "Bcell", any of the mask_100 = group_to_update[identity_col] == 100 rows '2nd_cell_type' = 'Bcell', put "yes"

import module001_removed_well_duplicate_compare_with_Bcell
from module001_removed_well_duplicate_compare_with_Bcell import classify_all_sequences

# Apply your function only to the Bcell subset

combined_df = classify_all_sequences(
    combined_df,
    groupby_col='1st_sequence_id',
    identity_col='nt_local_% of _identity',
    global_col='nt_local_% of _identity',
    subset_cell_type="Bcell"      # → 1st_cell_type must be Bcell
)

🔬 665753 rows where 1st_cell_type='Bcell' and 2nd_cell_type='Bcell'

──────────── Summary (by 1st_sequence_id) ────────────
unique_within_fna_Bcell   yes        → 585
unique_within_fna_Bcell   no         → 454

cross_donor_fna_Bcell     yes        → 1
cross_donor_fna_Bcell     no         → 12
cross_donor_fna_Bcell     yes and no → 165

cross_epitope_fna_Bcell   no         → 1
cross_epitope_fna_Bcell   yes and no → 96

───────────────────────────────────────────────────────

────────── Step 1 Summary ──────────
Sequences removed (identity==100 & nt_local_% of _identity>50, same well, lower read): 0
Remaining rows: 769445
Processed subset rows: 665753
────────────────────────────────────

from module001_add_multi_unique_in_one_well import add_multi_unique_in_one_well


combined_df = add_multi_unique_in_one_well(combined_df)
────────── multi_unique_in_one_well Summary ──────────
multi_unique_in_one_well
no     707398
yes     62047
Name: count, dtype: int64
──────────────────────────────────────────────────────

──────── Unique Wells per Category ────────
no    → 820
yes   → 40
───────────────────────────────────────────

combined_df_max = (
    combined_df
    .sort_values('nt_global_% of _identity', ascending=False)
    .groupby('1st_sequence_id', as_index=False)
    .head(1)
)

# optional-----------------
# combined_unique_withfna_df = combined_df[combined_df['unique_within_fna_Bcell'] == 'yes']
# combined_unique_withfna_df_max = (
#     combined_unique_withfna_df
#     .sort_values('nt_global_% of _identity', ascending=False)
#     .groupby('1st_sequence_id', as_index=False)
#     .head(1)
# )
# sub_df = combined_df[
#     (combined_df["nt_local_% of _identity"] == 100)
# ].copy()
# combined_sub_df = pd.concat([sub_df, combined_unique_withfna_df_max], ignore_index=True)
# unique_ids_count = combined_sub_df['1st_sequence_id'].nunique()
-----------------------


#Number of unique 1st_sequence_id: 1115
#Number of unique 1st_sequence_id of Bcell is 1039

combined_df_max.to_csv("/~/20251103_FNA_PBMC/Lambda/1115Lambda_chain_combined_df_max.csv",index= False)

20251031_FNA_PBMC_version4



