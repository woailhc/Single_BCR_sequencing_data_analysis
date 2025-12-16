#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 17:39:13 2025

@author: wuw5
"""
import pandas as pd

def add_multi_unique_in_one_well(
    df,
    well_col='1st_sequence_or_well_id',
    seq_col='1st_sequence_id',
    unique_flag_col='unique_within_fna_Bcell'
):
    """
    Adds a column 'multi_unique_in_one_well' that flags wells having
    multiple unique 1st_sequence_id values with 'yes' in unique_within_fna_Bcell.
    
    Logic:
        For each well (1st_sequence_or_well_id):
            - Count distinct sequence IDs where unique_within_fna_Bcell == 'yes'
            - If count > 1 → mark all rows of that well as 'yes'
            - Else → 'no'
    
    Also prints summary counts of unique wells by 'yes' or 'no'.
    """

    df = df.copy()

    # Filter to only 'yes' rows for counting
    valid = df[df[unique_flag_col] == "yes"]

    # Count unique 1st_sequence_id per well
    well_counts = valid.groupby(well_col)[seq_col].nunique()

    # Wells that have more than one "yes"
    multi_wells = well_counts[well_counts > 1].index

    # Create output flag column
    df["multi_unique_in_one_well"] = df[well_col].isin(multi_wells).map({True: "yes", False: "no"})

    # ─────────── SUMMARY ───────────
    print("\n────────── multi_unique_in_one_well Summary ──────────")
    print(df["multi_unique_in_one_well"].value_counts())
    print("──────────────────────────────────────────────────────")

    # Count unique wells by label
    well_summary = (
        df.drop_duplicates(subset=well_col)
          .groupby("multi_unique_in_one_well")[well_col]
          .nunique()
          .rename("unique_well_count")
    )

    print("\n──────── Unique Wells per Category ────────")
    for label, count in well_summary.items():
        print(f"{label:<5} → {count}")
    print("───────────────────────────────────────────\n")

    return df
