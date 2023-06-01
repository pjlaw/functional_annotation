#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats
from math import log10
from collections import Counter

# In[6]:


#annotate TF binding

#annotated_data.loc[pd.isna(wanted_mpra_data["TF_bound"]), "TF_bound"all_tf_scoresetscores(this_data):
def get_tf(this_data):
    if pd.isna(this_data):
        return("")
    #KLF5,1.20372,+#MAZ,1.20372,+#rst2,1.20372,+#S
    allhits=this_data.split("#")
    all_tf=set()
    for hit in allhits:
        tfname, score,strand=hit.split(",")
        all_tf.add(tfname)
    all_tf_str="#".join(all_tf)
    return(all_tf_str)

def get_abc_genes(this_data):
    if pd.isna(this_data):
        return("")

    all_genes=set()
    allhits=this_data.split("#")
    for hit in allhits:
        gene, score=hit.split(",")
        all_genes.add(gene)
    all_genes_str="#".join(all_genes)
    return(all_genes_str)

def get_eqtl(this_data):
    # 'rs12146099_ENSG00000135862#8.24204e-09#LAMC1,rs12146099_ENSG00000224468#0.00120194#LAMC1-AS1,rs12146099_ENSG00000116698#0.00343373#SMG7',
    if pd.isna(this_data):
        return("")
    allhits=this_data.split(",")
    all_eqtl=set()
    for hit in allhits:
        eqtl_name, pval, genename=hit.split("#")
        if genename=="nan":
            rsid,ensg=eqtl_name.split("_")
            all_eqtl.add(ensg)
        else:
            all_eqtl.add(genename)

    return(all_eqtl)


from math import isnan

def get_microC_tss(this_data):
    all_genes=set()
    if  pd.isna(this_data):
        return("")
    allhits=this_data.split(",")
    for hit in allhits:
        ensg, genename=hit.split("#")
        if genename == "nan":
            all_genes.add(ensg)
        else:
            all_genes.add(genename)
    all_genes_str="#".join(all_genes)
    return(all_genes_str)

def get_genes(gene_chrom, bin_low, bin_high):
    #get any genes that are near to the other end of the contact
    gene_list = gene_chrom.query(f'start<{bin_high} and end>{bin_low}')
    res=zip(gene_list["ENSG"].to_list(), gene_list["genename"].to_list())
    out=set()
    for geneid, genename in res:
        #print(genename, geneid)
        try:
            x=isnan(genename)
            out.add(geneid)
        except TypeError:
            #if isnan fails, it's because genename is a string, ie has a value
            out.add(genename)
    return(out)

#microC contacts from fithic
def find_microC_overlap(chrom, this_pos, cell_line):
    this_microC_data=pd.read_table(f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/fithic/fithic_juicer_version/merged/{cell_line}_inter_30_chr{chrom}_1000_merged.gz")
    #find the loops where snp of interest is in one of the ends
    #also only cis interactions
    snp_microC = this_microC_data.query(f"chr1==chr2 and ((bin1_low<={this_pos} and bin1_high>={this_pos}) or (bin2_low<={this_pos} and bin2_high>={this_pos}))")
    output_genes=set()
    if snp_microC.shape[0]==0:
        return(output_genes)

    #get the genes on this chromosome
    gene_chrom = gene_data.query(f"chr=='{chrom}'")
    microC_dict=snp_microC.to_dict("records")

    for res in microC_dict:
        #get the genes on the other end of the gwas interaction
        if res["bin1_low"]<=this_pos and res["bin1_high"]>=this_pos:
            this_gene=get_genes(gene_chrom, res["bin2_low"], res["bin2_high"])
        else:
            this_gene=get_genes(gene_chrom, res["bin1_low"], res["bin1_high"])

        output_genes = output_genes | this_gene

    output_genes_str="#".join(output_genes)
    return(output_genes_str)


def cleanup(sets_list):
    sets_list = sets_list.to_list()
    #if all none/empty
    if any(sets_list):
        return( set().union(*sets_list) )
    else:
        return(set())

def merge_all_celline(annot_data):
    annot_df = annot_data.values
    all_genes=[]
    #rows=varaints, cols=cell lines
    for row in annot_df:
        for anno in row:
            #check is not nan
            if anno == anno:
                this_genes = anno.split("#")
                all_genes+=this_genes
    #don't return set for TFs??
    return(set(all_genes))

def get_3d_structure(snp_chr, pos_start, pos_end, tad_data, compartment_data, tad_window=30000):
    this_chr_tads=tad_data.query(f"chrom=='chr{snp_chr}' and is_boundary_{tad_window}")
    #shift the insulation table so it represents the tads instead of the tad boundaries
    tads=pd.DataFrame({"start":[0]+this_chr_tads["end"].to_list(), "end":this_chr_tads["start"].to_list()+[chrom_sizes[f"chr{snp_chr}"]]})
    tads["tad_name"]=f"tad_{snp_chr}_"+tads.index.map(str)

    #snp region is entirely within a tad
    snp_tad=tads.query(f"start<={pos_start} and end>={pos_end}")
    if snp_tad.shape[0]==0:
        #print(this_gwas_snp, "is in a TAD boundary!")
        this_snp_tad="In_boundary"
    else:
        this_snp_tad=snp_tad.iloc[0,2] #get the tad name
        #print(this_snp_tad)

    #Look at E1 - if >0 = active, <0 inactive
    snp_compartment = compartment_data.query(f"chrom=='chr{snp_chr}' and start<={pos_end} and end>={pos_start}")
    if snp_compartment.shape[0]==0:
        #this shouldn't happen
        #print("not in a compartment")
        this_snp_compartment="overlaps"
    else:
        e1_val = snp_compartment["E1"].values
        #dichotonise the values across the bins
        comp_val = [1 if i>=0 else 0 for i in e1_val]
        comp_counts = Counter(comp_val)
        if comp_counts.most_common()[0][0]==1:
            this_snp_compartment="Active"
        else:
            this_snp_compartment="Inactive"

    return((this_snp_tad, this_snp_compartment))

def consensus_annot_calc(all_annots):
    annot_list=all_annots.values
    occurence_count = Counter(annot_list)
    most_common = occurence_count.most_common(1)
    common_n = most_common[0][1]
    common_ele = most_common[0][0]
    #if the most common element is 0, reutn 0
    if common_ele ==0:
        return(0)
    #at least half
    if common_n >= (len(annot_list)/2):
        return(common_ele)
    else:
        return(0)

chrom_sizes={"chr1":248956422,
            "chr2":242193529,
            "chr3":198295559,
            "chr4":190214555,
            "chr5":181538259,
            "chr6":170805979,
            "chr7":159345973,
            "chrX":156040895,
            "chr8":145138636,
            "chr9":138394717,
            "chr11":135086622,
            "chr10":133797422,
            "chr12":133275309,
            "chr13":114364328,
            "chr14":107043718,
            "chr15":101991189,
            "chr16":90338345,
            "chr17":83257441,
            "chr18":80373285,
            "chr20":64444167,
            "chr19":58617616,
            "chrY":57227415,
            "chr22":50818468,
            "chr21":46709983,
            "chrM":16569}

# In[8]:


#file created by /home/plaw/scripts/mpra/annotate_mpra_per_locus_with_consensus.py
#all_annotation_data=pd.read_table("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/annotated_data/CRC_all_annotation_scores_consensus.txt")
#

#
##gene information - start, end, name, tadname for each cell line
#gene_data = pd.read_table("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/annotated_data/gene_TAD_30kb.txt")
gene_data = pd.read_table("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/annotated_data/gene_TAD_30kb_all_MSS.txt")
#
#hcec_data = pd.read_table("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/annotated_data/HCEC-1CT_annotated.txt")
#hcec_data["TF_HCEC-1CT"] = hcec_data.apply(lambda x: get_tf(x["TF_bound"]), axis=1)
#
#all_annotation_data = all_annotation_data.merge(hcec_data[["rsid", "TF_HCEC-1CT"]], on="rsid")
#
cell_lines= [ "C32", "CL11" , "CACO2", "HT29", "SW480", "SW403", "SW948" ]
#for cell_name in cell_lines:
#    print("parsing", cell_name)
#    cell_data = pd.read_table(f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/annotated_data/{cell_name}_annotated.txt")
#    cell_data[f"TF_genes_{cell_name}"] = cell_data.apply(lambda x: get_tf(x["TF_bound"]), axis=1)
#    cell_data[f"ABC_genes_{cell_name}"] = cell_data.apply(lambda x: get_abc_genes(x["ABC"]), axis=1)
#
#    if cell_name != "C32":
#        cell_data[f"microC_TSS_{cell_name}"] = cell_data.apply(lambda x: get_microC_tss(x["microC_gene"]), axis=1)
#        cell_data[f"microC_gene_{cell_name}"] = cell_data.apply(lambda x: find_microC_overlap(x["chr"], x["pos_b38"], cell_name), axis=1)
#
#        insulation_table =  pd.read_csv(f'/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/microC/{cell_name}_10kb_TADs.tsv',sep='\t')
#        compartments = pd.read_csv(f'/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/microC/{cell_name}_50kb_ciseigs.tsv',sep='\t')
#
#        cell_data[f"TAD_{cell_name}"] = cell_data.apply(lambda x: get_3d_structure(x["chr"], x["pos_b38"], insulation_table, compartments), axis=1)
#        cell_data[[f"TAD_{cell_name}",f"compartment_{cell_name}"]] = cell_data[f"TAD_{cell_name}"].str.split("#", expand=True)
#
#        all_annotation_data = all_annotation_data.merge(cell_data[["rsid",f"TF_genes_{cell_name}", f"ABC_genes_{cell_name}",f"microC_TSS_{cell_name}", f"microC_gene_{cell_name}", f"TAD_{cell_name}", f"compartment_{cell_name}"]], on="rsid")
#    else:
#        all_annotation_data = all_annotation_data.merge(cell_data[["rsid",f"TF_genes_{cell_name}", f"ABC_genes_{cell_name}"]], on="rsid")
#
#
#
#all_annotation_data.to_csv("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/annotated_data/CRC_all_gene_annotation_data.tmp", sep="\t", index=False)

all_annotation_data=pd.read_csv("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/annotated_data/CRC_all_gene_annotation_data.tmp", sep="\t", low_memory=False)
print("tier1", all_annotation_data.query("Tier=='Tier1'").shape[0])
print("tier2", all_annotation_data.query("Tier=='Tier2'").shape[0])

#update the gwas snp names

def match_gwas_snps(snpname):
    full_gwas_name = gwas_setnames[gwas_setnames["gwas_snp"].str.contains(str(snpname))]
    #print(full_gwas_name)
    if full_gwas_name.shape[0]==0:
        #print(snpname, full_gwas_name)
        return(None)
    else:
        return(full_gwas_name["gwas_snp"].to_list()[0])

gwas_setnames = pd.DataFrame({'gwas_snp':all_annotation_data["gwas_snp"].unique()})
smr_res=pd.read_csv("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/CRC_GWAS/functional_data/all_smr_results_tcga_normal_merged.txt", sep="\t")
smr_res['gwas_snp_new'] = smr_res.apply(lambda x: match_gwas_snps(x["gwas_snp"]), axis=1)
filtered_smr = smr_res.query("`bonf_window_p_SMR.normal`<0.05 & (`p_HEIDI.normal`>0.05 | `p_HEIDI.normal`.isnull())", engine='python')
filtered_smr_tcga = smr_res.query("`bonf_window_p_SMR.tcga`<0.05 & (`p_HEIDI.tcga`>0.05 | `p_HEIDI.tcga`.isnull())", engine='python')

grouped_smr = filtered_smr.groupby("gwas_snp_new")['Gene'].apply(','.join).reset_index()
grouped_smr_tcga = filtered_smr_tcga.groupby("gwas_snp_new")['Gene'].apply(','.join).reset_index()


# In[9]:

joint_call_tads = pd.read_csv("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/microC/MSS_celllines_10kb_TADs.tsv", sep="\t")
joint_call_compartments = pd.read_csv("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/microC/MSS_celllines_10kb_ciseigs.tsv", sep="\t")


#try guess the causal gene(s)

all_tfs=[]
no_valid_genes=[]
all_valid_genes=[]
snp_in_tad_boundary=0
no_tier1=[]

gene_score_threshold=2 #needs at least 2 sources of annotation

# make a color map of fixed colors
cmap = mpl.colors.ListedColormap(['white','darkgrey', "blue", "red"])
bounds=[0,1,2]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

annot_labels=["microC TSS", "ABC", "SMR normal", "SMR TCGA", "microC","TAD"] # 6=microC_TSS, microC, ABC, eqtl, SMR, SMR_TCGA
n_annot=len(annot_labels)

indiv_snp_regions = pd.read_csv("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/US_UK_ASN_crc_all_gwas_snps_indiv_0.5MB_1e-05.txt", sep="\t")
gwas_snps = indiv_snp_regions["gwas_snp"].unique()
#gwas_snps =  ["rs6983267","rs4733767","rs5028523", "rs61776719"]
#for this_gwas_snp in ["rs6983267,rs7013278,rs4733767"]: #["rs5028523"]:
with open("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/gene_analysis/CRC_gene_predictions.txt", "w") as ofile:
    ofile.write("gwas_snp\tn_T1\t"+"\t".join(annot_labels[:-1])+"\tn_genes\tTF\tTAD\tcompartment\tputative_gene\tweak_gene_prediction\n")
    for this_gwas_snp in gwas_snps:
        print(this_gwas_snp)
        output_line=[this_gwas_snp]
        #wanted_snps = indiv_snp_regions.query("gwas_snp == '{this_gwas_snp}'")["rsid"].unique().values
        #only keep the tier 1 annotations (ie best variants)
        #this_snp_annotation = all_annotation_data.query(f"rsid in {wanted_snps} and Tier=='Tier1'")
        #this_snp_annotation = all_annotation_data.loc[(((all_annotation_data["Tier"]=='Tier1') | (all_annotation_data["Tier"]=='Tier2')) & (all_annotation_data["gwas_snp"].str.contains(this_gwas_snp))),:]
        this_snp_annotation = all_annotation_data.loc[( (all_annotation_data["Tier"]=='Tier1')  & (all_annotation_data["gwas_snp"].str.contains(this_gwas_snp)) ),:]

        num_tier1=this_snp_annotation.shape[0]
        output_line.append(str(num_tier1))

        if num_tier1==0:
            #print("No Tier1 variants for", this_gwas_snp)
            #look at all variants in this region, irrespective of score - is there some other annotation going on?
            this_snp_annotation = all_annotation_data.loc[ (all_annotation_data["gwas_snp"].str.contains(this_gwas_snp)) ,:]
            no_tier1.append(this_gwas_snp)
            no_valid_genes.append(this_gwas_snp)
            #continue

        annots=["microC_TSS", "microC_gene","ABC_genes"]
        all_annots_gene={}
        all_genes= set()
        for annot in annots:
            all_annots_gene[annot] = merge_all_celline(this_snp_annotation.loc[:, this_snp_annotation.columns.str.startswith(annot)])
            all_genes = all_genes | all_annots_gene[annot]
        output_line+=[",".join(all_annots_gene["microC_TSS"]),",".join(all_annots_gene["ABC_genes"]) ]

        #add smr results
        has_smr=grouped_smr["gwas_snp_new"].str.contains(this_gwas_snp)
        if any(has_smr):
            wanted=grouped_smr[has_smr]
            smr_genes=wanted["Gene"].to_list()[0].split(",")
            smr_genes=set(smr_genes)
            all_genes = all_genes | smr_genes
            output_line.append(",".join(smr_genes))
        else:
            smr_genes=set()
            output_line.append("")

        has_smr_tcga=grouped_smr_tcga["gwas_snp_new"].str.contains(this_gwas_snp)
        if any(has_smr_tcga):
            wanted=grouped_smr_tcga[has_smr_tcga]
            smr_tcga_genes=wanted["Gene"].to_list()[0].split(",")
            smr_tcga_genes = set(smr_tcga_genes)
            all_genes = all_genes | smr_tcga_genes
            output_line.append(",".join(smr_tcga_genes))
        else:
            smr_tcga_genes=set()
            output_line.append("")

        output_line.append(",".join(all_annots_gene["microC_gene"]))

        #get the final gene list
        all_genes = list(all_genes)
        all_genes.sort()

        #no genes in region
        n_genes = len(all_genes)
        if n_genes==0:
            no_valid_genes.append(this_gwas_snp)
            #continue

        output_line.append(str(n_genes))

        #get the LD region around the SNP
        region_pos = this_snp_annotation["pos_b38"].to_list()
        #get the tad and compartment it's in
        regions_tad, regions_compartment = get_3d_structure(this_snp_annotation["chr"].to_list()[0], min(region_pos), max(region_pos), joint_call_tads, joint_call_compartments)

        #add the TF data
        all_annots_gene["TF_genes"]=merge_all_celline(this_snp_annotation.loc[:, this_snp_annotation.columns.str.startswith("TF_genes")])
        output_line.append(",".join(all_annots_gene["TF_genes"]))

        tad_region_annotation = set()
        #if any of the cell lines show that gene is in the boundary
        if regions_tad=="In_boundary":
            tad_label="TAD*"
            tad_region_annotation.add("SNP_in_boundary")
        else:
            tad_label="TAD"

        if regions_compartment=="Active":
            tad_label+="+active"
        else:
            tad_label+="+inactive"
        annot_labels[-1]=tad_label

        #compartment consensus
        #compartment_cellline_consensus=this_snp_annotation.loc[:, this_snp_annotation.columns.str.startswith("compartment")].apply(lambda x: consensus_annot_calc(x), axis=1)
        #compartment_consensus=consensus_annot_calc(compartment_cellline_consensus)
        #all_annots_gene["compartment"]=compartment_consensus
        #output_line.append(compartment_consensus)
        all_annots_gene["compartment"]=regions_compartment
        output_line.append(regions_compartment)

        #create annotation matrix
        gene_diff_tad={k:0 for k in all_genes}
        this_output = np.zeros((len(all_genes), n_annot)) # 6=microC_TSS, microC, ABC, eqtl, SMR, SMR_TCGA

        snp_tad_boundary=False
        for i in range(len(all_genes)):
            this_gene=all_genes[i]
            if this_gene in all_annots_gene["microC_TSS"]:
                this_output[i,0]=1
            if this_gene in all_annots_gene["ABC_genes"]:
                this_output[i,1]=1
            if this_gene in smr_genes:
                this_output[i,2]=1
            if this_gene in smr_tcga_genes:
                this_output[i,3]=1
            if this_gene in all_annots_gene["microC_gene"]:
                this_output[i,4]=1

            #check which tad(s) the gene is in
            gene_pos = gene_data.query(f"genename =='{this_gene}'")
            if gene_pos.shape[0]==0:
                #unknown gene name
                #print(this_gene, "not found")
                continue

            this_gene_tad=gene_pos["MSS_celllines_TAD"].values[0]
            if this_gene_tad == regions_tad:
                 #in the same TAD
                 this_output[i,5]=0
            elif this_gene_tad=="In_boundary":
                this_output[i,5]=3
                tad_region_annotation.add("Gene_in_boundary")
            else:
                this_output[i,5]=2
                tad_region_annotation.add("Gene_in_different_TAD")
                #print("annotations differ)

            #same=[]
            #gene_tad_boundary=False

            #for each cell line, check if the gene TAD matches the gwas region TAD
            # for cell_name in cell_lines:
                # if cell_name == "C32":
                    # #no microC data for this cell line
                    # continue
                # #collapse all the TAD annotations for the variants in the region
                # cell_line_region_tad = set(this_snp_annotation[f"TAD_{cell_name}"].values)
                # if "In_boundary" in cell_line_region_tad:
                    # snp_tad_boundary=True

                # #is this gene in the same TAD as the variant
                # this_gene_tad=gene_pos[f"{cell_name}_TAD"].values[0]
                # if this_gene_tad == "In_boundary":
                    # gene_tad_boundary=True

                # #is the gene in the same TAD as the variant
                # if this_gene_tad in cell_line_region_tad:
                    # same.append(True)
                # else:
                    # same.append(False)

            # #collate the results per gene for across the cell lines
            # if sum(same)>=(len(same)/2):
                # #in the same TAD
                # this_output[i,5]=0
            # else:
                # this_output[i,5]=2
                # region_tad.add("SNP_in_different_TAD")
                # #print("annotations differ)

            # #if any of the cell lines show that gene is in the boundary
            # if gene_tad_boundary:
                # this_output[i,5]=3
                # region_tad.add("Gene_in_boundary")


            #is the gwas region in an active compartment
#            if all_annots_gene["compartment"]=="Active":
#                tad_label+="+active"
#            else:
#                tad_label+="+inactive"

#        #if any of the cell lines show that gene is in the boundary
#        if snp_tad_boundary:
#            tad_label="TAD*"
#            region_tad.add("SNP_in_boundary")
#        else:
#            tad_label="TAD"
#        annot_labels[-1]=tad_label

        #add the TAD annotation in second last place
        output_line.insert(-1, ",".join(tad_region_annotation))

        if n_genes > 0:
            #ADJUST SCORES HERE
            #dont include the microc gene or tad column in the scoring
            score_data = this_output[:,:4].copy()
            #up weight the microC score
            #score_data[score_data[:,0]>0,0]=2
            #find the genes that are in a different tad to the variant, lower the score
            #wanted_cols= (this_output[:,5]>0)
            #score_data[wanted_cols] = score_data[wanted_cols]*0.5
            gene_scores=np.sum(score_data, axis=1)
            max_score=np.max(gene_scores)
            out_genes=[]
            weak_genes=[]
            if max_score>=gene_score_threshold:
                max_idx=np.flatnonzero(gene_scores == max_score)
                for max_i in max_idx:
                    all_valid_genes.append(all_genes[max_i])
                    out_genes.append(all_genes[max_i])
            else:
                #no_valid_genes.append(this_gwas_snp)
                max_idx=np.flatnonzero(gene_scores == 1)
                for max_i in max_idx:
                    weak_genes.append(all_genes[max_i])
                out_genes.append("None")

            output_line.append(",".join(out_genes))
            output_line.append(",".join(weak_genes))
            ofile.write("\t".join(output_line)+"\n")

            # tell imshow about color map so that only set colors are used
            fig, ax = plt.subplots(1,1, figsize=(10,30))
            img = ax.imshow(this_output, interpolation='none',
                                cmap = cmap,vmax=3)

            ax.set_xticks(range(n_annot))
            ax.set_xticklabels(annot_labels, rotation=90)
            ax.set_yticks(np.arange(len(all_genes)))
            ax.set_yticklabels(all_genes)
            ax.set_title(this_gwas_snp)
            #plt.show()
            fig.savefig(f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/annotated_data/gene_figures/consensus/CRC_MSS_{this_gwas_snp}_gene_grid.png", facecolor="w", transparent=False, bbox_inches='tight')
            plt.close()
        else:
            #extra tabs for the genes
            ofile.write("\t".join(output_line)+"\t\t\n")

print(len(no_tier1), "with no Tier1 variants")
print(no_tier1)
print(len(no_valid_genes), "GWAS regions with no valid genes")
print(no_valid_genes)

valid_genes = sorted(list(set(all_valid_genes)))
with open("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/gene_analysis/CRC_gene_list.txt", "w") as ofile:
    for vg in valid_genes:
        ofile.write(vg+"\n")

