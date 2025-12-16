#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 20:55:27 2025

@author: wuw5
"""
import numpy as np
import pandas as pd


def classify_group(group, identity_col='identity'):
    # Initialize new columns
    group["FNA_shared_PBMC_cross_donor"] = ""
    group["FNA_shared_PBMC_cross_epitope"] = ""

    # Conditions
    mask_99 = group[identity_col] >= 99
    same_donor = group["first_real_subjectid"] == group["second_real_subjectid"]

    both_have_epitope = (
        group["first_flowindex_epitope"].notna() &
        group["second_flowindex_epitope"].notna()
    )

    same_epitope = both_have_epitope & (
        group["first_flowindex_epitope"] == group["second_flowindex_epitope"]
    )
    not_same_epitope = both_have_epitope & (
        group["first_flowindex_epitope"] != group["second_flowindex_epitope"]
    )

    # Per-row assignments
    group.loc[mask_99 & same_donor, "FNA_shared_PBMC_cross_donor"] = "no"
    group.loc[mask_99 & ~same_donor, "FNA_shared_PBMC_cross_donor"] = "yes"

    group.loc[mask_99 & same_epitope, "FNA_shared_PBMC_cross_epitope"] = "no"
    group.loc[mask_99 & not_same_epitope, "FNA_shared_PBMC_cross_epitope"] = "yes"

    # -- Fix for group consistency --
    for col in ["FNA_shared_PBMC_cross_donor", "FNA_shared_PBMC_cross_epitope"]:
        mask = mask_99 & group[col].isin(["yes", "no"])
        yes_rows = (group.loc[mask, col] == "yes").any()
        no_rows = (group.loc[mask, col] == "no").any()
        if yes_rows and no_rows:
            group.loc[mask, col] = "yes and no"

    return group

def classify_all_sequences(df, groupby_col='first_sequence_id', identity_col='identity'):
    """
    Apply the new classification to all groups (by first_sequence_id).
    Adds two new columns: FNA_shared_PBMC_cross_donor and FNA_shared_PBMC_cross_epitope.
    """
    result = (
        df.groupby(groupby_col, group_keys=False)
          .apply(lambda g: classify_group(g, identity_col=identity_col))
    )

    # Preserve original order
    result = result.loc[df.index]
    return result