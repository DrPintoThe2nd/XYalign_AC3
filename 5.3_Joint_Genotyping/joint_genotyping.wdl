#written by Samantha Zarate
#edited by Brendan J. Pinto

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
            chromosome = chromosome
    }

    output {
        File filtered_genotyped_VCF = filterVCF.filtered_genotyped_VCF
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

    command <<<
        tar -zvxf "~{genomicsDBtar}" -C .

        gatk \
            --java-options -Xmx8G \
            GenotypeGVCFs \
            -R "~{refFasta}" \
            -O "~{chromosome}.genotyped.vcf.gz" \
            -L "~{chromosome}" \
            -V "gendb://~{chromosome}" \
            --only-output-calls-starting-in-intervals
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
        File genotypedVCF = "~{chromosome}.genotyped.vcf.gz"
    }
}

task filterVCF {
    input {
        File refFasta
        File refIndex
        File refDict
        String chromosome
        File genotypedVCF
    }

    command <<<

        gatk \
            --java-options -Xmx8G \
            SelectVariants \
            -R "~{refFasta}" \
            -V "~{genotypedVCF}" \
            -L "~{chromosome}" \
            -O "~{chromosome}.filtered.genotyped.vcf.gz" --select-type-to-include SNP \
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
        File filtered_genotyped_VCF = "~{chromosome}.filtered.genotyped.vcf.gz"
    }
}