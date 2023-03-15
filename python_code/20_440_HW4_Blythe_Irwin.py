# PSET4
# Blythe Irwin
# Due 3/15/2023

## Import libraries

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from statannot import add_stat_annotation
import seaborn as sns


## Create dataframe for RNA-seq data

df = pd.read_csv('../data/GSE102152_FPKM_normalized_expression.txt.gz', sep= '\t') # Read csv to dataframe
# Rename columns to conditions
df = df.rename(columns={'G1':'Common FW1', 'G2':'Common FW2', 'G3':'Common FW3', 
                   'G4':'SR86 FW1', 'G5':'SR86 FW2', 'G6':'SR86 FW3', 'G7':'SR86 SW1', 
                   'G8':'SR86 SW2', 'G9':'SR86 SW3'})
df['Common FW Avg'] = df[['Common FW1','Common FW2', 'Common FW3']].mean(axis=1) # Add column for common FW mean
df['Common FW Std'] = df[['Common FW1','Common FW2', 'Common FW3']].std(axis=1) # Add column for common FW std
df['SR86 FW Avg'] = df[['SR86 FW1','SR86 FW2', 'SR86 FW3']].mean(axis=1) # Add column for SR86 FW mean
df['SR86 FW Std'] = df[['SR86 FW1','SR86 FW2', 'SR86 FW3']].std(axis=1) # Add column for SR86 FW std
df['SR86 SW Avg'] = df[['SR86 SW1','SR86 SW2', 'SR86 SW3']].mean(axis=1) # Add column for SR86 SW mean
df['SR86 SW Std'] = df[['SR86 SW1','SR86 SW2', 'SR86 SW3']].std(axis=1) # Add column for SR86 SW std


## Define bar plot function

# Defines a function to generate a bar plot for an input gene LOC
# Retrieves mean and std values for each condition at given gene loci from dataframe.
# Generates a bar plot with bonferroni corrected t tests for all three conditions:
# Common FW, SR86 FW, SR86 SW.
# Inputs:
# df - dataframe with RNA-seq data across all conditions
# gene - string, name of gene
# loc - string, gene locus name from dataframe
# Prints:
# Bar plot for RNA-seq data for given gene to png
# Returns:
# gene_avg - list of average FPKM for each condition
# gene_std - list of std of FPKM for each condition
def loc_to_bar(df,gene,loc):
    gene_avg = df.loc[df['gene_id'] == loc, ['Common FW Avg','SR86 FW Avg','SR86 SW Avg']].values.flatten().tolist()
    gene_std = df.loc[df['gene_id'] == loc, ['Common FW Std','SR86 FW Std','SR86 SW Std']].values.flatten().tolist()
    gene_df = pd.DataFrame(columns = ['Condition', 'FPKM'], index = [0, 1, 2, 3, 4, 5, 6, 7, 8])
    gene_df.loc[0] = ['Common FW'] + df.loc[df['gene_id'] == loc, ['Common FW1']].values.flatten().tolist() 
    gene_df.loc[1] = ['Common FW'] + df.loc[df['gene_id'] == loc, ['Common FW2']].values.flatten().tolist()
    gene_df.loc[2] = ['Common FW'] + df.loc[df['gene_id'] == loc, ['Common FW3']].values.flatten().tolist()
    gene_df.loc[3] = ['SR86 FW'] + df.loc[df['gene_id'] == loc, ['SR86 FW1']].values.flatten().tolist() 
    gene_df.loc[4] = ['SR86 FW'] + df.loc[df['gene_id'] == loc, ['SR86 FW2']].values.flatten().tolist()
    gene_df.loc[5] = ['SR86 FW'] + df.loc[df['gene_id'] == loc, ['SR86 FW3']].values.flatten().tolist()
    gene_df.loc[6] = ['SR86 SW'] + df.loc[df['gene_id'] == loc, ['SR86 SW1']].values.flatten().tolist() 
    gene_df.loc[7] = ['SR86 SW'] + df.loc[df['gene_id'] == loc, ['SR86 SW2']].values.flatten().tolist()
    gene_df.loc[8] = ['SR86 SW'] + df.loc[df['gene_id'] == loc, ['SR86 SW3']].values.flatten().tolist()
    gene_df['Condition'] = gene_df['Condition'].astype('category')
    gene_df['FPKM'] = gene_df['FPKM'].astype('float64')
    x = "Condition"
    y = "FPKM"
    order = ['Common FW', 'SR86 FW', 'SR86 SW']
    plt.figure()
    ax = sns.boxplot(data=gene_df, x=x, y=y, order=order)
    test_results = add_stat_annotation(ax, data=gene_df, x=x, y=y, order=order,
                                   box_pairs=[("Common FW", "SR86 FW"), ("SR86 FW", "SR86 SW"), ("Common FW", "SR86 SW")],
                                   test='t-test_ind', text_format='star', loc='inside', verbose=2)
    plt.title(gene + ' Expression')
    plt.savefig('../figures/' + gene + '_Expression_Bar_Plot.png')
    return [gene_avg, gene_std]


## Generate bar plots for genes of interest

Gn1a = loc_to_bar(df,'Gn1a','LOC_Os01g10110') # Generate bar plot for Gn1a gene
MOC1 = loc_to_bar(df,'MOC1','LOC_Os06g40780') # Generate bar plot for MOC1 gene
LOC_Os02g49700 = loc_to_bar(df,'LOC_Os02g49700','LOC_Os02g49700') # Generate bar plot for LOC_Os02g49700 gene
LOC_Os03g28300 = loc_to_bar(df,'LOC_Os03g28300','LOC_Os03g28300') # Generate bar plot for LOC_Os03g28300 gene
OsGAMYB = loc_to_bar(df,'OsGAMYB','LOC_Os01g59660') # Generate bar plot for OsGAMYB gene
OsCTR3 = loc_to_bar(df,'OsCTR3', 'LOC_Os04g52140') # Generate bar plot for OsCTR3 gene





