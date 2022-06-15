#written by Samantha Zarate
#modified by Brendan J. Pinto

version 1.0

workflow joint_calling {
    input {
        File refFasta
        File refIndex
        File refDict
        File genomicsDBtar
        String chromosome
    }

    call genotypeVCF {
        input:
            refFasta = refFasta,
            refIndex = refIndex,
            refDict = refDict,
            genomicsDBtar = genomicsDBtar,
            chromosome = chromosome
    }

    call filterVCF {
        input:
            refFasta = refFasta,
            refIndex = refIndex,
            refDict = refDict,
            genotypedVCF = genotypeVCF.genotypedVCF,
            genotypedVCF_tbi = genotypeVCF.genotypedVCF_tbi,
            chromosome = chromosome
    }

    output {
        File genotyped_VCF = genotypeVCF.genotypedVCF
        File genotyped_VCF_tbi = genotypeVCF.genotypedVCF_tbi
        File filtered_genotyped_VCF = filterVCF.filtered_genotyped_VCF
        File filtered_genotyped_VCF_tbi = filterVCF.filtered_genotyped_VCF_tbi
    }
}

task genotypeVCF {
    input {
        File refFasta
        File refIndex
        File refDict
        File genomicsDBtar
        String chromosome
    }

    String fastaName='~{basename(refFasta)}'
    String vcf_string = '~{basename(genomicsDBtar, "_genomicsDB.tar.gz")}'

    command <<<

        cp "~{refFasta}" .
        cp "~{refIndex}" .
        cp "~{refDict}" .
        tar -zvxf "~{genomicsDBtar}" -C .

        gatk \
            --java-options -Xmx8G \
            GenotypeGVCFs \
            -R "~{fastaName}" \
            -O "~{chromosome}.~{vcf_string}.genotyped.vcf.gz" \
            -L "~{chromosome}" \
            -V "gendb://~{vcf_string}_genomicsDB" \
            --only-output-calls-starting-in-intervals

    >>>

    runtime {
        docker : "drpintothe2nd/ac3_xysupp"
        disks : "local-disk 200 HDD"
        memory: "12G"
        cpu : 1
        preemptible: 1
        maxRetries: 1
    }

    output {
        File genotypedVCF = "~{chromosome}.~{vcf_string}.genotyped.vcf.gz"
        File genotypedVCF_tbi = "~{chromosome}.~{vcf_string}.genotyped.vcf.gz.tbi"
    }
}

task filterVCF {
    input {
        File refFasta
        File refIndex
        File refDict
        String chromosome
        File genotypedVCF
        File genotypedVCF_tbi
    }

    String fastaName='~{basename(refFasta)}'
    String inputVCF='~{basename(genotypedVCF)}'
    String filter_string = '~{basename(genotypedVCF, ".genotyped.vcf.gz")}'

    command <<<

        cp "~{genotypedVCF}" .
        cp "~{genotypedVCF_tbi}" .
        cp "~{refFasta}" .
        cp "~{refIndex}" .
        cp "~{refDict}" .

        gatk \
            --java-options -Xmx8G \
            SelectVariants \
            -R "~{fastaName}" \
            -V "~{inputVCF}" \
            -L "~{chromosome}" \
            -O "~{filter_string}.filtered.genotyped.vcf.gz" --select-type-to-include SNP \
            --restrict-alleles-to BIALLELIC \
            -select "AN >= 4 && MQ > 40.0 && QD > 7.0 && DP >= 10.0 && DP <= 1000.0"

    >>>

    runtime {
        docker : "drpintothe2nd/ac3_xysupp"
        disks : "local-disk 50 HDD"
        memory: "12G"
        cpu : 1
        preemptible: 1
        maxRetries: 1
    }

    output {
        File filtered_genotyped_VCF = "~{filter_string}.filtered.genotyped.vcf.gz"
        File filtered_genotyped_VCF_tbi = "~{filter_string}.filtered.genotyped.vcf.gz.tbi"
    }
}