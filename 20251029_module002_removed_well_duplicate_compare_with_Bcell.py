#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 12:09:31 2025

@author: wuw5
"""

import numpy as np
import pandas as pd

def classify_group(group, identity_col='identity', global_col='global'):
    """
    Classify rows for a given 1st_sequence_id group (Bcell context).
    """

    for col in ["unique_within_fna", "cross_donor_fna", "cross_epitope_fna"]:
        if col not in group.columns:
            group[col] = pd.Series("", index=group.index, dtype="string")
        else:
            group[col] = group[col].astype("string")

    g = group.copy()

    # Base masks
    mask_100 = g[identity_col] == 100
    mask_global = g[global_col] > 50
    mask_100_global = mask_100 & mask_global

    same_well = g["1st_sequence_or_well_id"] == g["2nd_sequence_or_well_id"]
    same_donor = g["1st_real_subjectid"] == g["2nd_real_subjectid"]

    both_have_epitope = g["1st_flowindex_epitope"].notna() & g["2nd_flowindex_epitope"].notna()
    same_epitope = both_have_epitope & (g["1st_flowindex_epitope"] == g["2nd_flowindex_epitope"])
    not_same_epitope = both_have_epitope & (g["1st_flowindex_epitope"] != g["2nd_flowindex_epitope"])

    # Step 1️⃣ Remove lower-read duplicates within same well
    removed_rows = 0
    mask_within = mask_100_global & same_well
    if mask_within.any():
        to_drop = []
        for well, subset in g.loc[mask_within].groupby("1st_sequence_or_well_id", group_keys=False):
            if subset.shape[0] > 1:
                maxr = subset["1st_mixcr_read"].max()
                to_drop.extend(subset.index[subset["1st_mixcr_read"] < maxr])
        if to_drop:
            removed_rows = len(to_drop)
            g = g.drop(index=to_drop)

    # Step 2️⃣ Unique within FNA
    g["unique_within_fna"] = np.where(mask_100_global.any(), "no", "yes")

    # Step 3️⃣ Cross‑donor
    g.loc[mask_100_global &  same_donor, "cross_donor_fna"] = "no"
    g.loc[mask_100_global & ~same_donor, "cross_donor_fna"] = "yes"
    donor_vals = [v for v in g.loc[mask_100_global, "cross_donor_fna"].unique() if v]
    if set(donor_vals) == {"yes", "no"}:
        g.loc[mask_100_global, "cross_donor_fna"] = "yes and no"

    # Step 4️⃣ Cross‑epitope
    g.loc[mask_100_global &  same_epitope,     "cross_epitope_fna"] = "no"
    g.loc[mask_100_global &  not_same_epitope, "cross_epitope_fna"] = "yes"
    donor_vals = [
    v for v in g.loc[mask_100_global, "cross_donor_fna"].unique()
        if pd.notna(v) and v != ""
        ]

    epi_vals = [
    v for v in g.loc[mask_100_global, "cross_epitope_fna"].unique()
        if pd.notna(v) and v != ""
        ]
    if set(epi_vals) == {"yes", "no"}:
        g.loc[mask_100_global, "cross_epitope_fna"] = "yes and no"

    return g, removed_rows


def classify_all_sequences(
    df,
    groupby_col='1st_sequence_id',
    identity_col='identity',
    global_col='global',
    subset_cell_type=None,
    cell_type_col_1='1st_cell_type',
    cell_type_col_2='2nd_cell_type'
):
    """
    Run Bcell classification group‑by‑group, rename all new columns with '_Bcell' suffix.
    """

    df = df.copy()
    total_removed = 0
    new_cols = ["unique_within_fna", "cross_donor_fna", "cross_epitope_fna"]

    # Filter: 2nd_cell_type must be 'Bcell'
    base_mask = df[cell_type_col_2] == "Bcell"

    if subset_cell_type is not None:
        mask = base_mask & (df[cell_type_col_1] == subset_cell_type)
        print(f"🔬 {mask.sum()} rows where {cell_type_col_1}='{subset_cell_type}' and {cell_type_col_2}='Bcell'")
    else:
        mask = base_mask
        print(f"🔬 {mask.sum()} rows where {cell_type_col_2}='Bcell'")

    work = df.loc[mask].copy()

    def run(g):
        nonlocal total_removed
        up, rem = classify_group(g, identity_col, global_col)
        total_removed += rem
        return up

    classified = work.groupby(groupby_col, group_keys=False).apply(run)
    classified = classified.loc[work.index]

    # give columns Bcell suffix before merging and summarizing
    # give columns <celltype> suffix before merging and summarizing
    suffix = f"_{subset_cell_type}" if subset_cell_type else "_all"
    renamed_cols = {c: f"{c}{suffix}" for c in new_cols}
    classified = classified.rename(columns=renamed_cols)


    df.loc[mask, list(renamed_cols.values())] = classified[list(renamed_cols.values())]

    # --- Summary per 1st_sequence_id ---
    summary = (
        df.loc[mask]
        .groupby(groupby_col)[list(renamed_cols.values())]
        .first()
    )

    print("\n──────────── Summary (by 1st_sequence_id) ────────────")
    for col in summary.columns:
        vc = summary[col].value_counts()
        for label in ["yes", "no", "yes and no"]:
            if label in vc:
                print(f"{col:<25s} {label:<10s} → {vc[label]}")
        print()
    print("───────────────────────────────────────────────────────")

    # --- Step 1 summary ---
    print("\n────────── Step 1 Summary ──────────")
    print(f"Sequences removed (identity==100 & {global_col}>50, same well, lower read): {total_removed}")
    print(f"Remaining rows: {df.shape[0]}")
    print(f"Processed subset rows: {mask.sum()}")
    print("────────────────────────────────────\n")

    return df