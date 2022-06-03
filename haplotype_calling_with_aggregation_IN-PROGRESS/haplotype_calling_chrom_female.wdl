#written by Samantha Zarate
#edited by Brendan Pinto

version 1.0

workflow haplotype_calling_chrom {
    input {
        File refFasta
        File fastaIndex
        File fastaDict
        File cram
        File cramIndex
        String sampleName
        String sex
		String chr8_output_filename
		String chrX_output_filename
    }

#    call haplotypeCaller as hc1 {
#        input:
#            refFasta = refFasta,
#            fastaIndex = fastaIndex,
#            fastaDict = fastaDict,
#            cram = cram,
#            cramIndex = cramIndex,
#            chrom = "chr1",
#            ploidy = 2,
#            sampleName = sampleName
#    }
#
#    call haplotypeCaller as hc2 {
#        input:
#            refFasta = refFasta,
#            fastaIndex = fastaIndex,
#            fastaDict = fastaDict,
#            cram = cram,
#            cramIndex = cramIndex,
#            chrom = "chr2",
#            ploidy = 2,
#            sampleName = sampleName
#    }
#
#    call haplotypeCaller as hc3 {
#        input:
#            refFasta = refFasta,
#            fastaIndex = fastaIndex,
#            fastaDict = fastaDict,
#            cram = cram,
#            cramIndex = cramIndex,
#            chrom = "chr3",
#            ploidy = 2,
#            sampleName = sampleName
#    }
#
#    call haplotypeCaller as hc4 {
#        input:
#            refFasta = refFasta,
#            fastaIndex = fastaIndex,
#            fastaDict = fastaDict,
#            cram = cram,
#            cramIndex = cramIndex,
#            chrom = "chr4",
#            ploidy = 2,
#            sampleName = sampleName
#    }
#
#    call haplotypeCaller as hc5 {
#        input:
#            refFasta = refFasta,
#            fastaIndex = fastaIndex,
#            fastaDict = fastaDict,
#            cram = cram,
#            cramIndex = cramIndex,
#            chrom = "chr5",
#            ploidy = 2,
#            sampleName = sampleName
#    }
#
#    call haplotypeCaller as hc6 {
#        input:
#            refFasta = refFasta,
#            fastaIndex = fastaIndex,
#            fastaDict = fastaDict,
#            cram = cram,
#            cramIndex = cramIndex,
#            chrom = "chr6",
#            ploidy = 2,
#            sampleName = sampleName
#    }
#
#    call haplotypeCaller as hc7 {
#        input:
#            refFasta = refFasta,
#            fastaIndex = fastaIndex,
#            fastaDict = fastaDict,
#            cram = cram,
#            cramIndex = cramIndex,
#            chrom = "chr7",
#            ploidy = 2,
#            sampleName = sampleName
#    }

    call haplotypeCaller as hc8 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr8",
            ploidy = 2,
            sampleName = sampleName
    }

#    call haplotypeCaller as hc9 {
#        input:
#            refFasta = refFasta,
#            fastaIndex = fastaIndex,
#            fastaDict = fastaDict,
#            cram = cram,
#            cramIndex = cramIndex,
#            chrom = "chr9",
#            ploidy = 2,
#            sampleName = sampleName
#    }
#
#    call haplotypeCaller as hc10 {
#        input:
#            refFasta = refFasta,
#            fastaIndex = fastaIndex,
#            fastaDict = fastaDict,
#            cram = cram,
#            cramIndex = cramIndex,
#            chrom = "chr10",
#            ploidy = 2,
#            sampleName = sampleName
#    }
#
#    call haplotypeCaller as hc11 {
#        input:
#            refFasta = refFasta,
#            fastaIndex = fastaIndex,
#            fastaDict = fastaDict,
#            cram = cram,
#            cramIndex = cramIndex,
#            chrom = "chr11",
#            ploidy = 2,
#            sampleName = sampleName
#    }
#
#    call haplotypeCaller as hc12 {
#        input:
#            refFasta = refFasta,
#            fastaIndex = fastaIndex,
#            fastaDict = fastaDict,
#            cram = cram,
#            cramIndex = cramIndex,
#            chrom = "chr12",
#            ploidy = 2,
#            sampleName = sampleName
#    }
#
#    call haplotypeCaller as hc13 {
#        input:
#            refFasta = refFasta,
#            fastaIndex = fastaIndex,
#            fastaDict = fastaDict,
#            cram = cram,
#            cramIndex = cramIndex,
#            chrom = "chr13",
#            ploidy = 2,
#            sampleName = sampleName
#    }
#
#    call haplotypeCaller as hc14 {
#        input:
#            refFasta = refFasta,
#            fastaIndex = fastaIndex,
#            fastaDict = fastaDict,
#            cram = cram,
#            cramIndex = cramIndex,
#            chrom = "chr14",
#            ploidy = 2,
#            sampleName = sampleName
#    }
#
#    call haplotypeCaller as hc15 {
#        input:
#            refFasta = refFasta,
#            fastaIndex = fastaIndex,
#            fastaDict = fastaDict,
#            cram = cram,
#            cramIndex = cramIndex,
#            chrom = "chr15",
#            ploidy = 2,
#            sampleName = sampleName
#    }
#
#    call haplotypeCaller as hc16 {
#        input:
#            refFasta = refFasta,
#            fastaIndex = fastaIndex,
#            fastaDict = fastaDict,
#            cram = cram,
#            cramIndex = cramIndex,
#            chrom = "chr16",
#            ploidy = 2,
#            sampleName = sampleName
#    }
#
#    call haplotypeCaller as hc17 {
#        input:
#            refFasta = refFasta,
#            fastaIndex = fastaIndex,
#            fastaDict = fastaDict,
#            cram = cram,
#            cramIndex = cramIndex,
#            chrom = "chr17",
#            ploidy = 2,
#            sampleName = sampleName
#    }
#
#    call haplotypeCaller as hc18 {
#        input:
#            refFasta = refFasta,
#            fastaIndex = fastaIndex,
#            fastaDict = fastaDict,
#            cram = cram,
#            cramIndex = cramIndex,
#            chrom = "chr18",
#            ploidy = 2,
#            sampleName = sampleName
#    }
#
#    call haplotypeCaller as hc19 {
#        input:
#            refFasta = refFasta,
#            fastaIndex = fastaIndex,
#            fastaDict = fastaDict,
#            cram = cram,
#            cramIndex = cramIndex,
#            chrom = "chr19",
#            ploidy = 2,
#            sampleName = sampleName
#    }
#
#    call haplotypeCaller as hc20 {
#        input:
#            refFasta = refFasta,
#            fastaIndex = fastaIndex,
#            fastaDict = fastaDict,
#            cram = cram,
#            cramIndex = cramIndex,
#            chrom = "chr20",
#            ploidy = 2,
#            sampleName = sampleName
#    }
#
#    call haplotypeCaller as hc21 {
#        input:
#            refFasta = refFasta,
#            fastaIndex = fastaIndex,
#            fastaDict = fastaDict,
#            cram = cram,
#            cramIndex = cramIndex,
#            chrom = "chr21",
#            ploidy = 2,
#            sampleName = sampleName
#    }
#    
#    call haplotypeCaller as hc22 {
#        input:
#            refFasta = refFasta,
#            fastaIndex = fastaIndex,
#            fastaDict = fastaDict,
#            cram = cram,
#            cramIndex = cramIndex,
#            chrom = "chr22",
#            ploidy = 2,
#            sampleName = sampleName
#    }

        call haplotypeCaller as hcX {
            input:
                refFasta = refFasta,
                fastaIndex = fastaIndex,
                fastaDict = fastaDict,
                cram = cram,
                cramIndex = cramIndex,
                chrom = "chrX",
                ploidy = 2,
                sampleName = sampleName
    }

  call MergeGVCFs {
    input:
      input_chr8_vcfs = hc8.hcVCF,
      input_chr8_vcfs_indexes = hc8.hcVCF_tbi,
      input_chrX_vcfs = hcX.hcVCF,
      input_chrX_vcfs_indexes = hcX.hcVCF_tbi,
      chr8_output_filename = chr8_output_filename,
      chrX_output_filename = chrX_output_filename
  }

    output {
#        File chr1_hcVCF = hc1.hcVCF
#        File chr2_hcVCF = hc2.hcVCF
#        File chr3_hcVCF = hc3.hcVCF
#        File chr4_hcVCF = hc4.hcVCF
#        File chr5_hcVCF = hc5.hcVCF
#        File chr6_hcVCF = hc6.hcVCF
#        File chr7_hcVCF = hc7.hcVCF
         File chr8_hcVCF = hc8.hcVCF
         File chr8_hcVCF_tbi = hc8.hcVCF_tbi
#        File chr9_hcVCF = hc9.hcVCF
#        File chr10_hcVCF = hc10.hcVCF
#        File chr11_hcVCF = hc11.hcVCF
#        File chr12_hcVCF = hc12.hcVCF
#        File chr13_hcVCF = hc13.hcVCF
#        File chr14_hcVCF = hc14.hcVCF
#        File chr15_hcVCF = hc15.hcVCF
#        File chr16_hcVCF = hc16.hcVCF
#        File chr17_hcVCF = hc17.hcVCF
#        File chr18_hcVCF = hc18.hcVCF
#        File chr19_hcVCF = hc19.hcVCF
#        File chr20_hcVCF = hc20.hcVCF
#        File chr21_hcVCF = hc21.hcVCF
#        File chr22_hcVCF = hc22.hcVCF
         File chrX_hcVCF = hcX.hcVCF
         File chrX_hcVCF_tbi = hcX.hcVCF_tbi
         File chr8_output_vcf = MergeGVCFs.output_chr8_vcf
         File chr8_output_vcf_index = MergeGVCFs.output_chr8_vcf_index
         File chrX_output_vcf = MergeGVCFs.output_chrX_vcf
         File chrX_output_vcf_index = MergeGVCFs.output_chrX_vcf_index
    }
}

task haplotypeCaller {
    input {
        File refFasta
        File fastaIndex
        File fastaDict
        File cram
        File cramIndex
        String chrom
        String sampleName
        Int ploidy
    }

    String cramName = '~{basename(cram)}'
    String fastaName = '~{basename(refFasta)}'
    String indexName = '~{basename(fastaIndex)}'
    String dictName = '~{basename(fastaDict)}'

    command <<<
        mv "~{cram}" .
        mv "~{cramIndex}" .
        mv "~{refFasta}" .
        mv "~{fastaIndex}" .
        mv "~{fastaDict}" .

        gatk HaplotypeCaller \
            --java-options "-Xmx24G -XX:+UseParallelGC -XX:ParallelGCThreads=$(nproc) -Djava.io.tmpdir=/dev/shm" \
            -R "~{fastaName}" \
            -I ./"~{cramName}" \
            -L "~{chrom}" \
            -pairHMM AVX_LOGLESS_CACHING \
            -O "~{sampleName}.~{chrom}.hc.vcf.gz" \
            -ERC GVCF \
            -ploidy "~{ploidy}" \
            -A Coverage \
            -A DepthPerAlleleBySample \
            -A DepthPerSampleHC \
            -A InbreedingCoeff \
            -A MappingQualityRankSumTest \
            -A MappingQualityZero \
            -A QualByDepth \
            -A ReadPosRankSumTest \
            -A RMSMappingQuality \
            -A StrandBiasBySample
    >>>

    Int diskGb = ceil(2.0 * size(cram, "G"))

    runtime {
        docker: "drpintothe2nd/ac3_xysupp"
        disks : "local-disk ${diskGb} SSD"
        memory: "24 GB"
        cpu : 4
        preemptible: 2
        maxRetries: 2
    }

    output {
        File hcVCF = "~{sampleName}.~{chrom}.hc.vcf.gz"
        File hcVCF_tbi = "~{sampleName}.~{chrom}.hc.vcf.gz.tbi"
    }
}

task MergeGVCFs {
  input {
    # Command parameters
    Array[File] input_chr8_vcfs
    Array[File] input_chr8_vcfs_indexes
    Array[File] input_chrX_vcfs
    Array[File] input_chrX_vcfs_indexes
    String chr8_output_filename
    String chrX_output_filename
  }
  
  command {
  set -e

    gatk --java-options "-Xmx24G"  \
      MergeVcfs \
      --INPUT ~{sep=' --INPUT ' input_chr8_vcfs} \
      --OUTPUT "~{chr8_output_filename}.gz"

    gatk --java-options "-Xmx24G"  \
      MergeVcfs \
      --INPUT ~{sep=' --INPUT ' input_chrX_vcfs} \
      --OUTPUT "~{chrX_output_filename}.gz"

  }
  runtime {
        docker: "drpintothe2nd/ac3_xysupp"
        disks : "local-disk 100 SSD"
        memory: "24 GB"
        cpu : 1
        preemptible: 2
        maxRetries: 2
  }
  output {
    File output_chr8_vcf = "~{chr8_output_filename}.gz"
    File output_chr8_vcf_index = "~{chr8_output_filename}.gz.tbi"
    File output_chrX_vcf = "~{chrX_output_filename}.gz"
    File output_chrX_vcf_index = "~{chrX_output_filename}.gz.tbi"
  }
}
