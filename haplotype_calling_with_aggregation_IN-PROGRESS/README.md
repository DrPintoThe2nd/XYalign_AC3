Current status:
```
womtool validate haplotype_calling_chrom_female.wdl
```

Failed to process workflow definition 'haplotype_calling_chrom' (reason 1 of 4): Failed to process 'call MergeGVCFs' 

(reason 1 of 4): Failed to supply input input_chrX_vcfs_indexes = hcX.hcVCF_tbi (reason 1 of 1): Cannot coerce expression of type 'File' to 'Array[File]'

Failed to process workflow definition 'haplotype_calling_chrom' (reason 2 of 4): Failed to process 'call MergeGVCFs' 

(reason 2 of 4): Failed to supply input input_chr8_vcfs = hc8.hcVCF (reason 1 of 1): Cannot coerce expression of type 'File' to 'Array[File]'

Failed to process workflow definition 'haplotype_calling_chrom' (reason 3 of 4): Failed to process 'call MergeGVCFs' 

(reason 3 of 4): Failed to supply input input_chrX_vcfs = hcX.hcVCF (reason 1 of 1): Cannot coerce expression of type 'File' to 'Array[File]'

Failed to process workflow definition 'haplotype_calling_chrom' (reason 4 of 4): Failed to process 'call MergeGVCFs' 

(reason 4 of 4): Failed to supply input input_chr8_vcfs_indexes = hc8.hcVCF_tbi (reason 1 of 1): Cannot coerce expression of type 'File' to 'Array[File]'
