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
    }

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

    output {
         File chr8_hcVCF = hc8.hcVCF
         File chr8_hcVCF_tbi = hc8.hcVCF_tbi
         File chrX_hcVCF = hcX.hcVCF
         File chrX_hcVCF_tbi = hcX.hcVCF_tbi
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
        disks : "local-disk ${diskGb} HDD"
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
