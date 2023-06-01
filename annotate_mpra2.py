import pandas as pd
import glob
import sys
#import sqlite3
import os.path
from math import log10

cell_line=sys.argv[1]
#cell_line = "HT29"
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


#mpra design file
#mpra_design = pd.read_table("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/CRC_MPRA_design_library_all_variants_v3_newgwas_with_controls", sep="\t")
#mpra_design_pruned = mpra_design.query("test_control == 'test'")
#mpra_design_pruned = mpra_design_pruned[["rsid", "pos", "test_control", "gwas_snp"]]


#mpra data (using mpraR)
#if cell_line=="HT29":
#  results_file=pd.read_table(f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/mpra_r/{cell_line}_mpra_toptable_no3.txt", sep="\t")
#else:
#  results_file=pd.read_table(f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/mpra_r/{cell_line}_mpra_toptable2.txt", sep="\t")
#
##split the rsid_dir into separate columns
#results_file[['snp', 'dir']] = results_file['rsid'].str.split('_', 1, expand=True) #maxsplit=1
#
##add the b38 data
#results_file = pd.merge(results_file, snp_data_pruned, how="left", on="snp")
##add the gwas snp that the mpra probe is associated with, and add b37 data
#results_file = pd.merge(results_file, mpra_design_pruned, how="left", left_on="snp", right_on="rsid")
#
##split off the controls and test data
#control_data = results_file.query('test_control=="control"')
#control_data.to_csv(f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/mpra_r/{cell_line}_mpra_controls.txt", sep="\t", index=False)
#
##only annotate the test data - some of the control data (scramble) don't have positions
#results_file = results_file.query('test_control=="test"')

#mpra data (using MPRAnalyze)
if cell_line in ["SW403", "HT29", "HCEC-1CT"]:
    mpra_file=pd.read_table(f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/mpraflow/{cell_line}/{cell_line}_mpranalyze_res_extraQC_all_scaled.txt", sep="\t")
    mpra_pruned=mpra_file[["rowname", "pval", "fdr", "logFC"]].rename(columns={"pval":"MPRA_P", "fdr":"MPRA_FDR", "logFC":"MPRA_logFC"})

    results_file = pd.merge(results_file, mpra_pruned, how="left", left_on="rsid", right_on="rowname")
    #remove duplicated columns from the joins
    results_file = results_file.drop(columns=["snp", "rowname"])

#chip
chip_path="/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/CRC_GWAS/functional_data/chip_data"
chip_annotations=["ATAC",  "CTCF",  "H3K27ac",  "H3K27me3",  "H3K36me3",  "H3K4me1",  "H3K4me3"]

def find_chip_overlap(chip_data, this_chrom, this_pos):
    overlaps = chip_data.query(f"@chip_data[0]=='chr{this_chrom}' and @chip_data[1]<={this_pos} and @chip_data[2]>={this_pos}")
    return(overlaps.shape[0])

if cell_line != "HCEC-1CT":
    for chip_anno in chip_annotations:
        chip_data = pd.read_table(f"{chip_path}/{cell_line}/{chip_anno}/{chip_anno}.consensus_peaks.bed", sep="\t", header=None)
        results_file[chip_anno] = results_file.apply(lambda x: find_chip_overlap(chip_data, x["chr"], x["pos_b38"]), axis=1)

'''
#extra chip
extra_studies=["GSE133928" , "GSE136888", "GSE156613" ]
for studyname in extra_studies:
    flist=glob.glob(f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/CRC_GWAS/functional_data/extra_chip/{studyname}/intersected/*_merged.bed")
    for fname in flist:
        chip_data = pd.read_table(fname, sep="\t", header=None)
        chip_name = os.path.basename(fname).replace("_merged.bed","")
        results_file[chip_name] = results_file.apply(lambda x: find_chip_overlap(chip_data, x["chr"], x["pos_b37"]), axis=1)

'''
#eqtl metaed from meta_eqtl.R
#eqtl_data=pd.read_table("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/CRC_GWAS/functional_data/complete_meta_intersect_mpra.txt", sep="\t")
#
#def find_eqtl_overlap(this_rsid):
#    overlaps = eqtl_data.query("SNP==@this_rsid").to_dict("records")
#
#    #convert the dataframe into a string
#    outstring=[]
#    for res in overlaps:
#        outstring.append(f'{res["SNP"]}_{res["gene"]}#{res["p.value"]}#{res["symbol"]}')
#    outstring=",".join(outstring)
#    return(outstring)
#
#results_file["eqtl"] = results_file.apply(lambda x: find_eqtl_overlap(x["rsid"]), axis=1)

#finemapping results
#finemap_data=pd.read_table("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/CRC_GWAS/finemapping/baseline_CRC_ASN_regions.txt", sep="\t")
#def find_finemap_snps(this_rsid):
#    contains=finemap_data["snps"].str.contains(this_rsid)
#    #boolean indexing to get the fields that contain the snp of interest
#    wanted=finemap_data.iloc[contains[contains].index].to_dict("records")
#    output=[]
#    for res in wanted:
#        #get the PIP for the SNP
#        all_snps=res["snps"].split(",")
#        all_pips=res["PIPs"].split(",")
#        try:
#            idx=all_snps.index(this_rsid)
#            output.append(all_pips[idx])
#        except ValueError:
#            #catch for whne the contains finds a partial match to the snp
#            continue
#    outstring=",".join(output)
#    return(outstring)
#
#results_file["finemap"] = results_file.apply(lambda x: find_finemap_snps(x["rsid"]), axis=1)


#TOBAIS (TF results)
def find_tf_overlap(tf_dat, this_chrom, this_pos):
    overlaps = tf_dat.query(f"@tf_dat[0]=='chr{this_chrom}' and @tf_dat[1]<={this_pos} and @tf_dat[2]>={this_pos}")
    if overlaps.shape[0]==0:
        return("")
    else:
        tf_overlaps = overlaps[3].str.split('_', n=1, expand=True) #tf_name,jaspar_id
        tf_names=tf_overlaps[0].to_list()
        scores=overlaps[9].to_list()
        strand=overlaps[5].to_list()

        overlap_res=zip(tf_names, scores, strand)
        output=[]
        for tf_name, tf_score, tf_strand in overlap_res:
            #if (tf_strand=="+" and this_dir=="fwd") or (tf_strand=="-" and this_dir=="rev"):
            output.append(f"{tf_name},{tf_score},{tf_strand}")

        tf_str="#".join(output)
        return(tf_str)

#prune results to only motifs found in humans
tf_data = pd.read_table(f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/CRC_GWAS/functional_data/ATAC_TF_footprinting/{cell_line}_all_bound.bed", sep="\t", header=None)
tf_file=open("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/reference_data/JASPAR2022_human_motifs.txt")
human_tf_ids=tf_file.readlines()
human_tf_ids=[s.rstrip() for s in human_tf_ids]
tf_file.close()
human_tf_data=tf_data.query(f"@tf_data[3] in @human_tf_ids")
results_file["TF_bound"] = results_file.apply(lambda x: find_tf_overlap(human_tf_data, x["chr"], x["pos_b38"]), axis=1)



#ABC (enhancer predictions)
def find_abc(abc_data, this_chrom, this_pos):
    if cell_line =="C32":
        score_col="powerlaw.Score"
    else:
        score_col="ABC.Score"
    overlaps = abc_data.query(f"chr=='chr{this_chrom}' and start<={this_pos} and end>={this_pos}")
    if overlaps.shape[0]==0:
        return("")
    else:
        gene_names=overlaps["TargetGene"].to_list()
        scores=overlaps[score_col].to_list()

        overlap_res=zip(gene_names, scores)
        output=[]
        for gname, abc_score in overlap_res:
            output.append(f"{gname},{abc_score}")

        tf_str="#".join(output)
        return(tf_str)


if cell_line == "SW480":
    #no rna seq data for this cell line
    abc_data = pd.read_table(f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/ABC/{cell_line}/Prediction/EnhancerPredictions.txt", sep="\t")
    results_file["ABC"] = results_file.apply(lambda x: find_abc(abc_data, x["chr"], x["pos_b38"]), axis=1)
elif cell_line == "SW403":
    abc_data = pd.read_table(f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw//ABC/ABC_with_rna/{cell_line}_alt/Prediction/EnhancerPredictions.txt", sep="\t")
    results_file["ABC"] = results_file.apply(lambda x: find_abc(abc_data, x["chr"], x["pos_b38"]), axis=1)
    #abc_data = pd.read_table(f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/ABC/{cell_line}_alt/Prediction_inputcorrected/EnhancerPredictions.txt", sep="\t")
    #results_file["ABC_corrected"] = results_file.apply(lambda x: find_abc(abc_data, x["chr"], x["pos_b38"]), axis=1)
elif cell_line != "HCEC-1CT":
    abc_data = pd.read_table(f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/ABC/ABC_with_rna/{cell_line}/Prediction/EnhancerPredictions.txt", sep="\t")
    results_file["ABC"] = results_file.apply(lambda x: find_abc(abc_data, x["chr"], x["pos_b38"]), axis=1)
    #abc_data = pd.read_table(f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/ABC/{cell_line}/Prediction_inputcorrected/EnhancerPredictions.txt", sep="\t")
    #results_file["ABC_corrected"] = results_file.apply(lambda x: find_abc(abc_data, x["chr"], x["pos_b38"]), axis=1)


#tss list
tss_window=500
tss_file=pd.read_table("/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/reference_data/Homo_sapiens.GRCh38.p13_gene_TSS_canonical.bed", sep="\t")
tss_file["TSS_min"]=tss_file["TSS"]-tss_window
tss_file["TSS_max"]=tss_file["TSS"]+tss_window

def get_genes(tss_chrom, bin_low, bin_high):
    #get any genes that are near to the other end of the contact
    gene_list = tss_chrom.query(f'TSS_min<{bin_high} and TSS_max>{bin_low}')
    res=zip(gene_list["geneID"].to_list(), gene_list["genename"].to_list())
    out=set()
    for geneid, genename in res:
        out.add(f"{geneid}#{genename}")
    return(out)

#microC contacts from fithic
def find_microC_overlap(chrom, this_pos):
    #if cell_line == "SW403":
    #    this_microC_data=pd.read_table(f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/fithic/fithic_cooler/merged{cell_line}_alt_inter_30_chr{chrom}_1000_merged.gz")
    #else:
    #    this_microC_data=pd.read_table(f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/fithic/fithic_cooler/merged/{cell_line}_inter_30_chr{chrom}_1000_merged.gz")
    this_microC_data=pd.read_table(f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/fithic/fithic_juicer_version/merged/{cell_line}_inter_30_chr{chrom}_1000_merged.gz")
    #find the loops where snp of interest is in one of the ends
    #also only cis interactions
    snp_microC = this_microC_data.query(f"chr1==chr2 and ((bin1_low<={this_pos} and bin1_high>={this_pos}) or (bin2_low<={this_pos} and bin2_high>={this_pos}))")
    if snp_microC.shape[0]==0:
        return(["", ""])

    #get the genes on this chromosome
    tss_file_chrom = tss_file.query(f"chromosome=='{chrom}'")

    microC_dict=snp_microC.to_dict("records")
    output=set()
    output_genes=set()
    for res in microC_dict:
        #get the coords of the loops
        #output in hic longrange format - chr1    111 222  chr2:333-444,55

        try:
            this_score=int(-log10(res["fdr"]))
        except ValueError:
            this_score="NA"
        #add both directions - needed for washu hic format
        outstring=f'chr{chrom}\t{res["bin1_low"]}\t{res["bin1_high"]}\tchr{chrom}:{res["bin2_low"]}-{res["bin2_high"]},{this_score}'
        output.add(outstring)

        #add the other end - for washu hic format
        outstring=f'chr{chrom}\t{res["bin2_low"]}\t{res["bin2_high"]}\tchr{chrom}:{res["bin1_low"]}-{res["bin1_high"]},{this_score}'
        output.add(outstring)

        #get the genes on the other end of the gwas interaction
        if res["bin1_low"]<=this_pos and res["bin1_high"]>=this_pos:
            this_gene=get_genes(tss_file_chrom, res["bin2_low"], res["bin2_high"])
        else:
            this_gene=get_genes(tss_file_chrom, res["bin1_low"], res["bin1_high"])

        '''
        #arrange it such that left side is always the end containing GWAS snp
        if res["bin1_low"]<=this_pos and res["bin1_high"]>=this_pos:
            outstring=f'chr{chrom}\t{res["bin1_low"]}\t{res["bin1_high"]}\tchr{chrom}:{res["bin2_low"]}-{res["bin2_high"]},{this_score}'
            this_gene=get_genes(tss_file_chrom, res["bin2_low"], res["bin2_high"])
        else:
            outstring=f'chr{chrom}\t{res["bin2_low"]}\t{res["bin2_high"]}\tchr{chrom}:{res["bin1_low"]}-{res["bin1_high"]},{this_score}'
            this_gene=get_genes(tss_file_chrom, res["bin1_low"], res["bin1_high"])

        '''
        output_genes = output_genes | this_gene

    all_contacts="#".join(output)
    gene_contacts = ",".join(output_genes)

    return([all_contacts, gene_contacts])

if cell_line != "HCEC-1CT" and cell_line != "C32":
    results_file["micro_res"] = results_file.apply(lambda x: find_microC_overlap(x["chr"], x["pos_b38"]), axis=1)
    results_file[['microC','microC_gene']] = pd.DataFrame(results_file.micro_res.tolist(), index= results_file.index)
    results_file.drop(columns=['micro_res'], inplace=True)


def get_chromhmm_annot(chromhmm_data, this_chrom, this_pos):
    overlaps = chromhmm_data.query(f"@chromhmm_data[0]=='chr{this_chrom}' and @chromhmm_data[1]<={this_pos} and @chromhmm_data[2]>{this_pos}")
    if overlaps.shape[0]==0:
        return("")
    else:
        chrom_annot = overlaps[3].to_list()
        #adjust the state number (+10) for the scatterplot
        chrom_annot = [i+10 for i in chrom_annot]
        #just return the first annotation - ignore the others
        return(chrom_annot[0])
        #if len(chrom_annot)>1:
        #    annots=",".join(map(str, set(chrom_annot)))
        #else:
        #    annots=str(chrom_annot[0])
        #return(annots)

if cell_line != "HCEC-1CT":
    chromhmm_file = f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/chromHMM/noCTCF/noCTCF_MSS_canonical_chr_15_states/{cell_line}_15_dense.bed"
    #chromhmm_file = f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/chromHMM/noCTCF/noCTCF_MPRA_canonical_chr_15_states/{cell_line}_15_dense.bed"
    chromhmm_data = pd.read_table(chromhmm_file, sep="\t", header=None, skiprows=1)
    results_file["chromHMM"] = results_file.apply(lambda x: get_chromhmm_annot(chromhmm_data, x["chr"], x["pos_b38"]), axis=1)

#output
results_file.to_csv(f"/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/MPRA/annotated_data/{cell_line}_annotated.txt", sep="\t", index=False)


