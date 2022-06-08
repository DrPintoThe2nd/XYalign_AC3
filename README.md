# XYalign_AC3

This repository stores custom WDLs for our project analyzing data using the Terra platform.

Step outlines:

#Extracting reads previously mapped to GRCh38.

(1) Converting CRAM files to BAM files (samtools, CRAM-to-BAM_multithreaded.wdl)

(2) Stripping reads from the BAM files, pairing them and trimming them, with QC (samtools, bbmap, and trim_galore, StripReadsFromBams.wdl)

#Map trimmed reads to new genome assembly (CHM13).

(3) Convert paired fastqs to CRAM files (bwa and samtools, t2t_alignment.wdl)

#Convert CRAM file to gVCF for each individual (and potentially aggregate gVCFs for downstream genotyping.

(4) Use GATK HaplotypeCaller (and MergeVcfs?) on chr8 and chrX for each individual (gatk, haplotype_calling_chrom_female.wdl)

#Call variants by joint genotyping (3-steps)

(5.1) Generate a sample map for GATK's GenomicsDBImport function.

(5.2) Generate a genomicsDB using GATK's GenomicsDBImport function.

(5.3) Joint genotyping and filtering using GATK's GenotypeGVCFs function.

#Calculate VCF stats....

(6) Profit?
