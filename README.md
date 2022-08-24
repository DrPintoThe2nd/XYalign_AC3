# XYalign_AC3

This repository stores custom WDLs for our project analyzing data using the Terra platform.

Step outlines:

#Extracting reads previously mapped to GRCh38.

(1) Converting CRAM files to BAM files (samtools, CRAM-to-BAM_multithreaded.wdl).

(2) Stripping reads from the BAM files, pairing them and trimming them, with QC (samtools, bbmap, and trim_galore, StripReadsFromBams.wdl).

#Map trimmed reads to new genome assembly (CHM13).

(3) Convert paired fastqs to CRAM files (bwa and samtools, t2t_alignment.wdl).

#Convert CRAM file to gVCF for each individual (and potentially aggregate gVCFs for downstream genotyping.

(4) Use GATK HaplotypeCaller on chr8 and chrX for each individual (gatk, haplotype_calling_chrom_female.wdl).

#Call variants by joint genotyping (3-steps)

(5.1) Generate a sample map for GATK's GenomicsDBImport function (custom Broad script, generate-sample-map.wdl).

(5.2) Generate a genomicsDB using GATK's GenomicsDBImport function (gatk, t2t_genomics_db.wdl).

(5.3) Joint genotyping and filtering using GATK's GenotypeGVCFs function (gatk, joint_genotyping.wdl).

(5.3.1) Calculate final VCF stats (rtg-tools and bcftools, misc/calc_stats.wdl)

(6.0.1) Generate appropriate salmon references using 6.1-salmon/salmon_index.sh (salmon).

(6.1) Run salmon to generate quant files (salmon, salmon.wdl).

(6.1.1) Extract TPM data from salmon output for each individual from both references and subtract the "raw" from "xy" values into a seperate column (process_quants.wdl).

(6.1.2) In  any order, run three collate_XXX WDLs to combine TPM (or TPM differences) data from all individuals for both reference genomes into a single file for each reference (collate_xy_quants.wdl/collate_raw_quants.wdl/collate_quant_diffs.wdl).

(6.0.2) Generate appropriate HiSat2 reference for RNAseq alignment (hisat2).

(6.2) 


