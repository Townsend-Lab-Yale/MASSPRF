# Updated 2025-10-15 by Yide Jin
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
from os.path import isfile, join
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# -------------------------
# Helpers
# -------------------------

RESULTS_ANCHOR = '//Results based on model averaging of gamma using different models'
ABBREV_ANCHOR  = 'Abbreviation: MS=Model Selection; MA=Model Averaging; CI=Confidence Interval'

def find_block(lines):
    """Return start/end indices (exclusive) for the MASS-PRF results table block."""
    start = None
    end = None
    for i, ln in enumerate(lines):
        if RESULTS_ANCHOR in ln:
            start = i + 1  # header line is next
            break
    if start is None:
        raise ValueError("Could not find results header anchor in file.")
    for j in range(start, len(lines)):
        if lines[j].startswith(ABBREV_ANCHOR):
            end = j
            break
    if end is None:
        # Fallback: stop at first completely blank line after header
        for j in range(start, len(lines)):
            if not lines[j].strip():
                end = j
                break
    if end is None:
        end = len(lines)
    return start, end

def parse_table(lines):
    """Parse the MASS-PRF results table into a DataFrame with robust whitespace handling."""
    s, e = find_block(lines)
    block = [ln.rstrip() for ln in lines[s:e] if ln.strip()]  # drop empties
    if not block:
        raise ValueError("Results block is empty.")
    header_line = block[0]
    cols = re.split(r'\s+', header_line.strip())
    rows = []
    for ln in block[1:]:
        toks = re.split(r'\s+', ln.strip())
        # Normalize token count to header length
        if len(toks) < len(cols):
            toks += [''] * (len(cols) - len(toks))
        elif len(toks) > len(cols):
            toks = toks[:len(cols)]
        rows.append(toks)
    df = pd.DataFrame(rows, columns=cols)
    # Coerce numerics where applicable (ignoring non-numeric like R/S/* in the last two cols)
    for c in ['Position', 'DivergentTime', 'Gamma', 'Lower_CI_Gamma', 'Upper_CI_Gamma']:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce')
    return df

def extract_scaling(lines):
    """Infer scaling factor from MASS-PRF text. Defaults to 1 if not found."""
    scaling = None
    for ln in lines:
        if "Scaling by supplied factor of" in ln:
            try:
                scaling = int(ln.split("Scaling by supplied factor of")[-1].strip())
            except Exception:
                pass
        elif "No Scaling requested" in ln:
            scaling = 1
    if scaling is None:
        # default to 1 but don't crash
        scaling = 1
    return scaling

def descale_to_nt(df, scaling):
    """Expand each scaled site row into nucleotide-resolution rows (position 1..N)."""
    # If scaling is 1, we still want a copy to operate uniformly
    out_cols = list(df.columns)
    # Build expanded frame
    expanded = []
    for _, row in df.iterrows():
        for _ in range(int(scaling)):
            expanded.append(row.values.tolist())
    out = pd.DataFrame(expanded, columns=out_cols)
    # Rebuild 1..N position (generic 'Position')
    if 'Position' in out.columns:
        out = out.drop(columns=['Position'])
    out.insert(0, 'Position', list(range(1, len(out) + 1)))
    return out

def write_lists(proc_dir, neg_l, pos_l, str_pos_l, boring_l, failed_l):
    with open(join(proc_dir, "Negative_genes.txt"), "w") as f:
        f.write('\n'.join(neg_l) + '\n')
    with open(join(proc_dir, "Positive_genes.txt"), "w") as f:
        f.write('\n'.join(pos_l) + '\n')
    with open(join(proc_dir, "Strongly_positive_genes.txt"), "w") as f:
        f.write('\n'.join(str_pos_l) + '\n')
    with open(join(proc_dir, "Boring_genes.txt"), "w") as f:
        f.write('\n'.join(boring_l) + '\n')
    with open(join(proc_dir, "Failed_genes.txt"), "w") as f:
        f.write('\n'.join(failed_l) + '\n')

# -------------------------
# Main
# -------------------------

def main():
    wd = os.getcwd()
    file_list = [f for f in os.listdir(wd) if isfile(join(wd, f)) and f.lower().endswith(".txt")]

    proc_dir = './processed_output/'
    plot_dir = './processed_output/plots/'
    os.makedirs(proc_dir, exist_ok=True)
    os.makedirs(plot_dir, exist_ok=True)

    boring_l, neg_l, pos_l, failed_l, str_pos_l = [], [], [], [], []

    for entry in file_list:
        # Expect filenames like NAME_SUFFIX.txt -> use first two tokens for gene name
        name_parts = entry.split('_')
        gene_name = os.path.splitext(os.path.basename(entry))[0]

        try:
            with open(entry, encoding='utf-8', errors='ignore') as fh:
                text = fh.read()
            if 'Mission accomplished' not in text:
                print(f"Failed: {gene_name} (no 'Mission accomplished')")
                failed_l.append(gene_name)
                continue

            lines = text.splitlines()
            # Parse table
            table = parse_table(lines)

            # Scaling: expand to per-nucleotide to unify downstream
            scaling = extract_scaling(lines)
            sel_df = descale_to_nt(table, scaling)

            # Significance filters (same as before)
            pos_df = sel_df[(sel_df['Gamma'] > 4) & (sel_df['Lower_CI_Gamma'] > 0)]
            str_pos_df = sel_df[(sel_df['Gamma'] > 4) & (sel_df['Lower_CI_Gamma'] > 4)]
            neg_df = sel_df[(sel_df['Gamma'] < -1) & (sel_df['Upper_CI_Gamma'] < 0)]

            # Always write the full per-NT CSV
            out_csv = join(proc_dir, f"{gene_name}.csv")
            sel_df.to_csv(out_csv, index=False)

            # Site lists
            if not pos_df.empty:
                pos_l.append(gene_name)
                with open(join(proc_dir, f"{gene_name}_pos_sites.txt"), "w") as f:
                    f.write('\n'.join(map(lambda x: str(int(x)), pos_df['Position'].tolist())) + '\n')
                if not str_pos_df.empty:
                    str_pos_l.append(gene_name)
                    with open(join(proc_dir, f"{gene_name}_str_pos_sites.txt"), "w") as f:
                        f.write('\n'.join(map(lambda x: str(int(x)), str_pos_df['Position'].tolist())) + '\n')
            if not neg_df.empty:
                neg_l.append(gene_name)
                with open(join(proc_dir, f"{gene_name}_neg_sites.txt"), "w") as f:
                    f.write('\n'.join(map(lambda x: str(int(x)), neg_df['Position'].tolist())) + '\n')
            if pos_df.empty and neg_df.empty:
                boring_l.append(gene_name)
                print("Boring: " + gene_name)

            # Plot
            try:
                import numpy as np
                x = sel_df['Position'].to_numpy(dtype=float)
                y = sel_df['Gamma'].to_numpy(dtype=float)
                ylo = sel_df['Lower_CI_Gamma'].to_numpy(dtype=float)
                yhi = sel_df['Upper_CI_Gamma'].to_numpy(dtype=float)
                plt.plot(x, y)
                plt.fill_between(x, ylo, yhi, alpha=0.2)
                plt.axhline(y=0, color='black', linestyle='-')
                plt.xlabel('position')  # generic label for NT or AA
                plt.ylabel('Scaled selection coefficient (gamma)')
                plt.savefig(join(plot_dir, f"{gene_name}.pdf"))
                plt.close()
            except Exception as plot_err:
                print(f"Plot failed for {gene_name}: {plot_err}")

        except Exception as e:
            print(f"Failed: {gene_name} ({e})")
            failed_l.append(gene_name)

    # Final lists
    write_lists(proc_dir, neg_l, pos_l, str_pos_l, boring_l, failed_l)

if __name__ == "__main__":
    main()

