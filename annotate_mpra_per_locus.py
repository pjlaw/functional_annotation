#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats
from math import log10
import glob
import os.path


def calc_abc_median(annot_scores):
    all_scores=[]
    for i in annot_scores:
        if pd.notna(i):
            allhits=i.split("#")
            for hit in allhits:
                gene, score=hit.split(",")
                all_scores.append(float(score))
    return(np.nanmedian(all_scores))
    #plt.hist(all_scores)
    #plt.show()

def get_abc_scores(this_data, threshold):
    if pd.isna(this_data):
        return(0)
    strong_hit=[]
    allhits=this_data.split("#")
    for hit in allhits:
        gene, score=hit.split(",")
        score=float(score)
        if float(score)>threshold:
            strong_hit.append(True)
        else:
            strong_hit.append(False)

    if any(strong_hit):
        return(2)
    else:
        return(1)

def getscores(this_data, score_threshold):
        if pd.isna(this_data):
            return(0)
        #KLF5,1.20372,+#MAZ,1.20372,+#rst2,1.20372,+#S
        allhits=this_data.split("#")
        strong_hit=[]
        for hit in allhits:
            tfname, score,strand=hit.split(",")
            if float(score)>score_threshold:
                strong_hit.append(True)
            else:
                strong_hit.append(False)

        if any(strong_hit):
            return(2)
        else:
            return(1)

def get_pvals(this_data):
    if pd.isna(this_data):
        return(0)
    #ENSG00000135862_rs3768623#6.01888e-05,ENSG00000224468_rs3768623#1.09509e-05
    allhits=this_data.split(",")
    strong_hit=[]
    weak_hit=[]
    for hit in allhits:
        genersid, pval, genename=hit.split("#")
        pval=float(pval)
        if pval<eqtl_threshold:
            strong_hit.append(True)
        elif pval<0.05:
            weak_hit.append(False)

    if any(strong_hit):
        return(2)
    elif all(weak_hit):
        return(1)
    else:
        return(0)

def get_tf_score_median(tf_scores):
    all_tf_scores=[]
    for i in tf_scores:
        if pd.notna(i):
            allhits=i.split("#")
            for hit in allhits:
                tfname, score,strand=hit.split(",")
                all_tf_scores.append(float(score))

    #plt.hist(all_tf_scores)
    #plt.show()

    #print(np.nanquantile(all_tf_scores, np.linspace(0,1,10)))
    #print(np.nanmedian(all_tf_scores))
    #print(np.median(all_tf_scores))
    score_threshold=np.nanmedian(all_tf_scores)
    return(score_threshold)

def annotate_cellline(annotated_data):
    data_matrix = annotated_data[['rsid','chromHMM']].copy()
    #annotate MPRA
    if 'MPRA_FDR' in annotated_data:
        data_matrix["MPRA"]=0
        data_matrix.loc[annotated_data['MPRA_FDR']<=1e-3, "MPRA"]=2
        data_matrix.loc[annotated_data['MPRA_FDR'].between(1e-3,0.01), "MPRA"]=1

    #annotate ChIP
    data_matrix["ChIP"]=0

    #putative promoter or enhancer
    data_matrix.loc[annotated_data['H3K4me3'].astype("bool") | annotated_data['H3K4me1'].astype("bool"), "ChIP"]=1
    #together with H3K27ac makes it a strong hit
    data_matrix.loc[annotated_data['H3K27ac'].astype("bool") & annotated_data['H3K4me1'].astype("bool"), "ChIP"]=2
    data_matrix.loc[annotated_data['H3K27ac'].astype("bool") & annotated_data['H3K4me3'].astype("bool"), "ChIP"]=2

    #annotate ATAC and CTCF
    data_matrix["ATAC"]=0
    data_matrix.loc[annotated_data['ATAC']==1, "ATAC"]=2

    data_matrix["CTCF"]=0
    data_matrix.loc[annotated_data['CTCF']==1, "CTCF"]=2

    #annotate microC
    #has a contact = hit, contact with gene promoter = strong hit
    if 'microC' in annotated_data:
        data_matrix["microC"]=0
        data_matrix.loc[annotated_data["microC"].notna(), "microC"]=1
        data_matrix.loc[annotated_data["microC_gene"].notna(), "microC"]=2

    #annotate TF binding
    #has a tf = hit, tf with score greater than threshold = strong hit
    #get all the scores

    score_threshold=get_tf_score_median(annotated_data["TF_bound"])
    data_matrix["TF"] = annotated_data.apply(lambda x: getscores(x["TF_bound"], score_threshold), axis=1)

    #annotate ABC
    abc_threshold=calc_abc_median(annotated_data["ABC"])
    data_matrix["ABC"] = annotated_data.apply(lambda x: get_abc_scores(x["ABC"], abc_threshold), axis=1)


    #data_matrix["chromHMM"]=data_matrix["chromHMM"].astype("int")
    return(data_matrix)

#annotation files created by /home/plaw/scripts/mpra/annotate_mpra2.py
#rsid    chr     pos_b37 P_value gwas_snp        LLR     pos_b38 MPRA_P  MPRA_FDR        MPRA_logFC      TF_bound [various_chip_marks ABC microC chromHMM ]

hcec_data = pd.read_table("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/annotated_data/HCEC-1CT_annotated.txt")
hcec_data["MPRA_HCEC-1CT"]=0
hcec_data.loc[hcec_data['MPRA_FDR']<=1e-3, "MPRA_HCEC-1CT"]=2
hcec_data.loc[hcec_data['MPRA_FDR'].between(1e-3,0.01), "MPRA_HCEC-1CT"]=1

score_threshold=get_tf_score_median(hcec_data["TF_bound"])
hcec_data["TF_HCEC-1CT"] = hcec_data.apply(lambda x: getscores(x["TF_bound"], score_threshold), axis=1)

#all_annotated_data = all_annotated_data.drop(columns=["chr","pos_b37","P_value","gwas_snp","LLR","pos_b38" ])
all_annotated_data = hcec_data[["rsid", "MPRA_HCEC-1CT", "TF_HCEC-1CT"]].copy()

cell_lines= [ "C32", "CL11", "CACO2", "HT29", "SW480", "SW403", "SW948" ]
for cell_name in cell_lines:
    print("parsing", cell_name)
    cell_data = pd.read_table(f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/annotated_data/{cell_name}_annotated.txt")
    cell_data = annotate_cellline(cell_data)
    this_columns = cell_data.columns
    new_columns=["rsid"]+[f"{col}_{cell_name}" for col in this_columns[1:]]
    cell_data.columns = new_columns

    all_annotated_data = all_annotated_data.merge(cell_data, on="rsid")

#add in the non-cell line specific data
#file created by /home/plaw/scripts/mpra/annotate_mpra2_noncellspecific.py
#rsid    chr     pos_b37 P_value gwas_snp        LLR     snp     pos_b38 cyto    gwas_with_cyto  gwas_with_cyto_for_sorting      finemap eqtl    SMR     SMR_TCGA      akita
noncellspecific_data = pd.read_table(f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/annotated_data/CRC_MPRA_noncellspecfic_annotated.txt")

print("adding non-cell specific data")
all_annotated_data = all_annotated_data.merge(noncellspecific_data, on="rsid")

#annotate finemapping
#found in finemapping = hit, PIP>0.5 = strong hit
all_annotated_data["finemap_score"]=0
all_annotated_data.loc[all_annotated_data['finemap']>0, "finemap_score"]=1
all_annotated_data.loc[all_annotated_data['finemap']>0.5, "finemap_score"]=2

#annotate eqtl

from math import log10
all_p=[]
for i in all_annotated_data["eqtl"]:
    if pd.notna(i):
        allhits=i.split(",")
        for hit in allhits:
            genersid, pval, genename=hit.split("#")
            pval=float(pval)
            all_p.append(-log10(pval))

#plt.hist(all_p)
#plt.show()

#print(np.nanquantile(all_p, np.linspace(0,1,10)))
#print(np.nanmedian(all_p))
#print(np.median(all_p))

# greater than median
print("calculating scores")
eqtl_threshold=10**(-np.nanmedian(all_p)) #1e-4
all_annotated_data["eqtl_score"] = all_annotated_data.apply(lambda x: get_pvals(x["eqtl"]), axis=1)

annotation_cols =["finemap_score","eqtl_score","SMR","SMR_TCGA","akita"]#, "ABC_corrected"]
multi_annot=["MPRA","ChIP","ATAC","microC","TF", "ABC"]
#get the annotations that come from multiple cell lines
for annot in multi_annot:
    annotation_cols += [col for col in all_annotated_data if col.startswith(annot)]
all_annotated_data['score'] = all_annotated_data[annotation_cols].sum(axis=1)

#stratify
print("statify scores")
all_annotated_data["Tier"]="Tier3"
score_deciles = np.quantile(all_annotated_data['score'], np.linspace(0,1,10))

all_annotated_data.loc[all_annotated_data['score']>=score_deciles[7], "Tier"]="Tier1"
all_annotated_data.loc[all_annotated_data['score'].between(score_deciles[4],score_deciles[7], inclusive="left"), "Tier"]="Tier2"

all_annotated_data = all_annotated_data.drop(columns=["finemap", "eqtl" ])
all_annotated_data.to_csv("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/annotated_data/CRC_all_annotation_scores.txt", sep="\t", index=False)

#all_annotated_data=pd.read_table("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/annotated_data/CRC_all_annotation_scores.txt")


######################################
#plotting
print("plotting")
#unique labels, but keep order
unique_labels=list(dict.fromkeys(all_annotated_data.sort_values("gwas_with_cyto_for_sorting")["gwas_with_cyto"].to_list()))

# boxplot = all_annotated_data.boxplot("score", by="gwas_with_cyto_for_sorting", figsize=(30,10), rot=90, fontsize=8, grid=False)
# boxplot.set_xticks(np.arange(1,len(unique_labels)+1))
# boxplot.set_xticklabels(unique_labels)
# plt.suptitle("")

# plt.hist(all_annotated_data['score'])
# plt.show()

# print(score_deciles[4]) #5th decile, ie median
# print(score_deciles[7]) #8th decile


# mpra_counts = all_annotated_data["MPRA"].value_counts()
# print("Num of MPRA variants by signficance")
# print(mpra_counts)
# mpra_counts.plot.bar(stacked=True)

# variants_per_locus = annotated_data["gwas_snp"].value_counts()
# print("variants per locus")
# variants_per_locus.hist()

# mpra_positive = all_annotated_data.query("MPRA>0")["gwas_snp"].unique()
# total_loci = len(all_annotated_data["gwas_snp"].unique())
# print(len(mpra_positive), "loci with MPRA hit, out of", total_loci)
# tier1_mpra = set(all_annotated_data.query("MPRA==1")["gwas_snp"].to_list())
# tier2_mpra = set(all_annotated_data.query("MPRA==2")["gwas_snp"].to_list())

# print("Tier1 only", len(tier1_mpra-tier2_mpra))
# print("Tier2 only", len(tier2_mpra))

# print("significant hits per locus")
# mpra_loci_counts = all_annotated_data.groupby("gwas_snp")["MPRA"].value_counts()


tier_counts=all_annotated_data.groupby("gwas_with_cyto_for_sorting")["Tier"].value_counts(normalize=True).mul(100)
barplot=tier_counts.unstack().plot.bar(stacked=True, figsize=(30,10), rot=90, fontsize=8,color=["yellowgreen", "bisque", "lightgrey"])

barplot.set_xticks(np.arange(len(unique_labels)))
x=barplot.set_xticklabels(unique_labels)
plt.savefig(f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/annotated_data/figures/CRC_all_tier_counts.svg", facecolor="w", transparent=False, bbox_inches='tight')


gwas_snps = all_annotated_data["gwas_snp"].unique()

plot_annotation_cols =["finemap_score","eqtl_score","SMR","SMR_TCGA","akita"]
multi_annot=["MPRA","ChIP","ATAC","CTCF","microC","TF", "ABC", "chromHMM"]
#get the annotations that come from multiple cell lines
for annot in multi_annot:
    plot_annotation_cols += [col for col in all_annotated_data if col.startswith(annot)]
#plot_annotation_cols =[ "MPRA","ChIP","ATAC","CTCF", "microC_gene","finemap","TF","eqtl","ABC", "chromHMM"] #"ABC_corrected",

# make a color map of fixed colors - extra colours are the chromHMM annotations - will need to rearrange them for each cell line

#cmap = mpl.colors.ListedColormap(['white', 'lightgrey','black'])
#full MSS chromHMM profile
cmap = mpl.colors.ListedColormap([(1,1,1), (0.7,0.7,0.7), (0,0,0),
                             (0.0, 0.4, 0.0),
                             (0.0, 0.8, 0.2),
                             (0.4, 1.0, 0.0),
                             (0.4, 1.0, 0.0),
                             (1.0, 1.0, 0.0),
                             (1.0, 0.2, 0.2),
                             (1.0, 0.0, 0.0),
                             (0.8, 0.8, 0.4),
                             (1.0, 1.0, 1.0),
                             (0.82, 0.82, 0.82),
                             (0.82, 0.82, 0.82),
                             (0.8, 0.6, 1.0),
                             (1.0, 1.0, 0.8),
                             (1.0, 0.2, 0.2),
                             (1.0, 0.2, 0.2)])

#MPRA cell lines
#                                  (0.82, 0.82, 0.82),
#                                    (1.0, 1.0, 1.0),
#                                    (0.8, 0.6, 1.0),
#                                    (1.0, 1.0, 0.8),
#                                    (0.4, 1.0, 0.0),
#                                    (1.0, 1.0, 1.0),
#                                    (0.4, 1.0, 0.0),
#                                    (0.0, 0.4, 0.0),
#                                    (0.0, 0.8, 0.2),
#                                    (0.4, 1.0, 0.0),
#                                    (1.0, 0.2, 0.2),
#                                    (1.0, 0.2, 0.2),
#                                    (1.0, 0.0, 0.0),
#                                    (1.0, 0.2, 0.2),
#                                    (1.0, 1.0, 0.0) ])




bounds=[0, 0.5, 1.5, 4, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

#for positive_snp in mpra_positive: #[:5]:
for gwas_snp in gwas_snps: #[:5]: ["rs6983267,rs7013278,rs4733767"]:#
    print(gwas_snp)
    this_data = all_annotated_data.query(f"gwas_snp == '{gwas_snp}'").sort_values(['pos_b38'])
    #print(this_data)
    snp_names=this_data["rsid"].to_list()
    #print(len(snp_names))
    this_data = this_data[plot_annotation_cols].to_numpy().transpose()

    fig, ax = plt.subplots(figsize=(40,8))
    img = ax.imshow(this_data,interpolation='nearest',
                    cmap = cmap,norm=norm)

    ax.set_xticks(np.arange(len(snp_names)))
    ax.set_xticklabels(snp_names, fontsize="10", rotation=90)
    ax.set_yticks(np.arange(len(plot_annotation_cols)))
    ax.set_yticklabels(plot_annotation_cols, fontsize="12")
    ax.set_title(gwas_snp)

    split_snps = gwas_snp.split(",")
    for snpid in split_snps:
        try:
            i=snp_names.index(snpid)
            plt.gca().get_xticklabels()[i].set_color('red')
        except ValueError:
            continue

    #plt.show()
    fig.savefig(f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/annotated_data/figures/CRC_all_{gwas_snp}_annotation_grid.png", facecolor="w", transparent=False, bbox_inches='tight')
    plt.close()

