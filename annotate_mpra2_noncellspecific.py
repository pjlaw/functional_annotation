'''
Parse the data for the non-cell line specific data sets - eQTL, SMR normal/tumour, finemapping results, chromHMM, (akita?)
'''

import pandas as pd
import glob
import sys
#import sqlite3
import os.path
from math import log10
import numpy as np

#b38 snp positions
results_file=pd.read_table("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/MPRA_LLR_snps.txt")
results_file = results_file.rename(columns={"pos": "pos_b37"})

#adjust for the snps missing the gwas snp (outside of the 500kb window)
results_file.loc[results_file["gwas_snp"].isna(), "gwas_snp"]="rs55810369,rs1445012"
results_file.loc[results_file["rsid"]=="rs62489410", "gwas_snp"]="rs17686932"

#results_file = results_file.drop_duplicates(subset=['rsid', 'gwas_snp'])
results_file = results_file.drop_duplicates(subset=['rsid'])
snp_data = pd.read_table("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/snp_positions_b38.bed", sep="\t", header=None, names=["chr", "start", "end", "snp"])
snp_data_pruned = snp_data[["snp", "start"]]
snp_data_pruned = snp_data_pruned.drop_duplicates()
snp_data_pruned = snp_data_pruned.rename(columns={"start": "pos_b38"})
#add the b38 pos
results_file = pd.merge(results_file, snp_data_pruned, how="outer", left_on="rsid", right_on="snp")

#add the cytoband information
cytofile=pd.read_table("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/reference_data/cytoBand_b38.txt.gz", sep="\t", header=None, names=["chrom", "start", "end", "cyto", "band"])
def find_cytoband(cyto_data, this_chrom, this_pos):
    overlaps = cyto_data.query(f"chrom=='chr{this_chrom}' and start<={this_pos} and end>={this_pos}")
    if overlaps.shape[0]==0:
        return("")
    else:
        #only take the first row
        all_cytobands=overlaps.iloc[:,3].to_list()
        all_cytobands.sort()
        cyto_string=str(this_chrom)+all_cytobands[0]
        return(cyto_string)

results_file["cyto"] = results_file.apply(lambda x: find_cytoband(cytofile, x["chr"], x["pos_b38"]), axis=1)

results_file["gwas_with_cyto"]=results_file["cyto"]+"_"+results_file["gwas_snp"]
#find the cases where the locus falls across multiple cytobands, pick the more common one
snp_counts = results_file.groupby("gwas_snp")["gwas_with_cyto"].nunique()
multi_cyto=snp_counts[snp_counts>1]
for multi in multi_cyto.index:
    max_cyto = results_file.query(f"gwas_snp == '{multi}'").groupby("gwas_with_cyto")["gwas_with_cyto"].count().idxmax()
    results_file.loc[results_file["gwas_snp"]==multi, "gwas_with_cyto"]=max_cyto

results_file["gwas_with_cyto_for_sorting"]=np.where(results_file["chr"].astype(str).str.len()==1, "0"+results_file["gwas_with_cyto"], results_file["gwas_with_cyto"])

#finemapping results
finemap_data=pd.read_table("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/CRC_GWAS/finemapping/baseline_CRC_ASN_regions.txt", sep="\t")
def find_finemap_snps(this_rsid):
    contains=finemap_data["snps"].str.contains(this_rsid)
    #boolean indexing to get the fields that contain the snp of interest
    wanted=finemap_data.iloc[contains[contains].index].to_dict("records")
    output=[]
    for res in wanted:
        #get the PIP for the SNP
        all_snps=res["snps"].split(",")
        all_pips=res["PIPs"].split(",")
        try:
            idx=all_snps.index(this_rsid)
            output.append(all_pips[idx])
        except ValueError:
            #catch for whne the contains finds a partial match to the snp
            continue
    outstring=",".join(output)
    return(outstring)

results_file["finemap"] = results_file.apply(lambda x: find_finemap_snps(x["rsid"]), axis=1)

#eqtl metaed from meta_eqtl.R
eqtl_data=pd.read_table("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/CRC_GWAS/functional_data/complete_meta_intersect_mpra.txt", sep="\t")

def find_eqtl_overlap(this_rsid):
    overlaps = eqtl_data.query("SNP==@this_rsid").to_dict("records")

    #convert the dataframe into a string
    outstring=[]
    for res in overlaps:
        outstring.append(f'{res["SNP"]}_{res["gene"]}#{res["p.value"]}#{res["symbol"]}')
    outstring=",".join(outstring)
    return(outstring)

results_file["eqtl"] = results_file.apply(lambda x: find_eqtl_overlap(x["rsid"]), axis=1)



#SMR
smr_res=pd.read_csv("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/CRC_GWAS/functional_data/all_smr_results_tcga_normal_merged.txt", sep="\t")

#update the gwas snp names - grouped based on gwas snp (smr snp might not be used in the mpra analysis)
gwas_setnames = pd.DataFrame({'gwas_snp':results_file["gwas_snp"].unique()})
def match_gwas_snps(snpname):
    full_gwas_name = gwas_setnames[gwas_setnames["gwas_snp"].str.contains(str(snpname))]
    #print(full_gwas_name)
    if full_gwas_name.shape[0]==0:
        #print(snpname, full_gwas_name)
        return(None)
    else:
        #print(full_gwas_name["gwas_snp"].to_list()[0])
        return(full_gwas_name["gwas_snp"].to_list()[0])
#
smr_res['gwas_snp_new'] = smr_res.apply(lambda x: match_gwas_snps(x["gwas_snp"]), axis=1)
filtered_smr = smr_res.query("`bonf_window_p_SMR.normal`<0.05 & (`p_HEIDI.normal`>0.01 | `p_HEIDI.normal`.isnull())", engine='python')
filtered_smr_tcga = smr_res.query("`bonf_window_p_SMR.tcga`<0.05 & (`p_HEIDI.tcga`>0.01 | `p_HEIDI.tcga`.isnull())", engine='python')

grouped_smr = filtered_smr.groupby("gwas_snp_new")['Gene'].apply(','.join).reset_index()
grouped_smr_tcga = filtered_smr_tcga.groupby("gwas_snp_new")['Gene'].apply(','.join).reset_index()

results_file["SMR"]=0
results_file.loc[results_file['gwas_snp'].isin(grouped_smr["gwas_snp_new"].to_list()), "SMR"]=2
results_file["SMR_TCGA"]=0
results_file.loc[results_file['gwas_snp'].isin(grouped_smr_tcga["gwas_snp_new"].to_list()), "SMR_TCGA"]=2

#from akita analysis - max discruption score>0.5 (median) for any bp within 100bp of SNP
disrupted_snps=['rs2070699', 'rs77148098', 'rs704417', 'rs1554865', 'rs11692435', 'rs34963268', 'rs1773860', 'rs6065668', 'rs10849434', 'rs6584283', 'rs12659017', 'rs6983267', 'rs10978941', 'rs6066825', 'rs1782645', 'rs1800734', 'rs2732875', 'rs80158569', 'rs1791373', 'rs9983528', 'rs7946853', 'rs1800469', 'rs12603526', 'rs983318', 'rs4901473', 'rs4919687', 'rs7810512', 'rs5028523', 'rs56324967', 'rs2001732', 'rs7300312', 'rs2450115', 'rs4444073', 'rs151127921', 'rs17094983', 'rs2527927', 'rs62042090', 'rs2208603', 'rs28840750', 'rs847208', 'rs12427846', 'rs704017', 'rs35564340', 'rs61975764', 'rs1426947', 'rs9271363', 'rs130651', 'rs1446585', 'rs2735940', 'rs11557154', 'rs994308', 'rs6059938', 'rs3087967', 'rs7542665', 'rs17686932', 'rs12979278', 'rs10936599', 'rs55810369', 'rs61776719', 'rs3217810', 'rs3217874', 'rs6012915', 'rs1951864', 'rs4813802', 'rs653178', 'rs10006803', 'rs10817106', 'rs7606562', 'rs113569514', 'rs9614460', 'rs7859362', 'rs1775910', 'rs35204860', 'rs2155065', 'rs7623129', 'rs13831', 'rs7071258', 'rs145997965', 'rs11789898', 'rs6928864', 'rs9924886', 'rs7299936', 'rs116353863', 'rs16892766', 'rs3809570', 'rs16878812', 'rs67550176']
#annotate MPRA
results_file["akita"]=0
results_file.loc[results_file['rsid'].isin(disrupted_snps), "akita"]=2

#output
results_file.to_csv(f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/annotated_data/CRC_MPRA_noncellspecfic_annotated.txt", sep="\t", index=False)


