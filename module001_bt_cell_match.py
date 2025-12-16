#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 17:03:45 2025

@author: wuw5
"""

import pandas as pd
import numpy as np

def add_bt_cell_match(
    df,
    identity_col='identity',
    global_col='global',
    cell_type_col_1='1st_cell_type',
    cell_type_col_2='2nd_cell_type',
    global_threshold=50
):
    """
    Adds a 'B-T CELL MATCH' column to the dataframe.

    Logic:
        - If identity == 100 AND global > global_threshold
          AND 1st_cell_type != 2nd_cell_type:
              → "WARNING: B-T CELL 100% NT MATCH"
        - Otherwise:
              → "NA"
    """

    df = df.copy()
    # Clean 'identity' column:
    df[identity_col] = pd.to_numeric(df[identity_col], errors="coerce")  # convert to numeric, NaN if invalid

    # Define condition: only valid if identity == 100 (not NaN, not empty)
  

    # Default to "NA"
    df["B-T CELL MATCH"] = "NA"

    # Define condition
    cond_bt_warning = (
        (df[identity_col] == 100)
        & (df[global_col] > global_threshold)
        & (df[cell_type_col_1] != df[cell_type_col_2])
    )

    # Apply results
    df.loc[cond_bt_warning, "B-T CELL MATCH"] = "WARNING: B-T CELL 100% NT MATCH"

    # Optional summary output
    print("────────── B-T CELL MATCH SUMMARY ──────────")
    print(df["B-T CELL MATCH"].value_counts())
    print("────────────────────────────────────────────")

    return df