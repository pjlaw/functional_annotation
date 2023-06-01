#!/bin/bash
#SBATCH -J abc
#SBATCH -t 5-0
#SBATCH -A prepay-houlston
#SBATCH -p compute
#SBATCH -c 8
#SBATCH -o /home/plaw/scripts/tmp/abc_%A_%a.o
# --array=0-7

# -o /home/plaw/scripts/tmp/basenji_%A_%a.o %j.o
#HT29 SW480 SW403_alt
#DLD1 HCA7 HCT116 
#array=( DLD1  HCT116 CL11 CACO2 SW948 HT29 SW480 SW403_alt  )
#idx=${SLURM_ARRAY_TASK_ID}
#cell_line="${array[$idx]}"
cell_line="SW480"
#cell_line="SW403_alt"

#conda activate /data/scratch/DGE/DUDGE/MOPOPGEN/plaw/conda_envs/ABC
module load bedtools
module load SAMtools
module load java

#outpath=/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/ABC/ABC_with_rna/"$cell_line"/
outpath=/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/ABC/"$cell_line"/
mkdir -p $outpath

ABC_DIR="/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/programs/ABC-Enhancer-Gene-Prediction/"
DATA_DIR="/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/CRC_GWAS/functional_data"

#using the "full" chromosome sizes file here for the peaks that may have aligned there, filtered out in later steps
python $ABC_DIR/src/makeCandidateRegions.py \
	--narrowPeak $DATA_DIR/ATAC/"$cell_line".mRp.clN_summits.bed \
	--bam $DATA_DIR/bams/ATAC/"$cell_line".mRp.clN.sorted.bam \
	--outDir $outpath/Peaks/ \
	--chrom_sizes $ABC_DIR/reference/hg38/chr_sizes_b38_full.txt \
	--regions_blocklist /data/scratch/DGE/DUDGE/MOPOPGEN/plaw/basenji/data/hg38-blacklist.v2.bed \
	--regions_includelist $ABC_DIR/reference/hg38/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed \
	--peakExtendFromSummit 250 \
	--nStrongestPeaks 175000  

python $ABC_DIR/src/run.neighborhoods.py \
	--candidate_enhancer_regions $outpath/Peaks/"$cell_line".mRp.clN_summits.bed.candidateRegions.bed \
	--genes $ABC_DIR/reference/hg38/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.bed \
	--H3K27ac $DATA_DIR/bams/"$cell_line"_H3K27ac_R1.mLb.clN.sorted.bam,$DATA_DIR/bams/"$cell_line"_H3K27ac_R2.mLb.clN.sorted.bam \
	--ATAC $DATA_DIR/bams/ATAC/"$cell_line".mRp.clN.sorted.bam \
	--chrom_sizes $ABC_DIR/reference/hg38/chr_sizes_b38_full.txt \
	--ubiquitously_expressed_genes $ABC_DIR/reference/UbiquitouslyExpressedGenesHG19.txt \
	--cellType $cell_line \
	--outdir $outpath/Neighborhoods 
#	--expression_table /data/scratch/DGE/DUDGE/MOPOPGEN/plaw/CRC_GWAS/functional_data/rnaseq_tpm/"$cell_line"_tpm.txt

#input corrected bams
#--outdir $outpath/Neighborhoods_inputcorrected
#--H3K27ac $DATA_DIR/bams/"$cell_line"_H3K27ac_R1_corrected.bigwig,$DATA_DIR/bams/"$cell_line"_H3K27ac_R2_corrected.bigwig \

python $ABC_DIR/src/predict.py \
	--enhancers $outpath/Neighborhoods/EnhancerList.txt \
	--genes $outpath/Neighborhoods/GeneList.txt \
	--outdir $outpath/Prediction/ \
	--HiCdir /data/scratch/DGE/DUDGE/MOPOPGEN/plaw/ABC/hic_data/"$cell_line"/raw/ \
	--hic_type juicebox \
	--hic_resolution 1000 \
	--scale_hic_using_powerlaw \
	--threshold .02 \
	--cellType $cell_line \
	--chrom_sizes $ABC_DIR/reference/hg38/chr_sizes_b38.txt \
	--run_all_genes \
	--make_all_putative

#--enhancers $outpath/Neighborhoods_inputcorrected/EnhancerList.txt \
#--genes $outpath/Neighborhoods_inputcorrected/GeneList.txt \
#--outdir $outpath/Prediction_inputcorrected/ \
#--chromosomes "chr22" #"all"
#--chromosomes "chr22" #"all"



