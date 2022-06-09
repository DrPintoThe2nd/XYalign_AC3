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

    command <<<

        cp "~{refFasta}" .
        cp "~{refIndex}" .
        cp "~{refDict}" .
        tar -zvxf "~{genomicsDBtar}" -C .

        gatk \
            --java-options -Xmx8G \
            GenotypeGVCFs \
            -R "~{fastaName}" \
            -O "~{chromosome}.genotyped.vcf.gz" \
            -L "~{chromosome}" \
            -V "gendb://~{chromosome}_GenomicsDB" \
            --only-output-calls-starting-in-intervals

#            tabix -p vcf "~{chromosome}.genotyped.vcf.gz"
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
        File genotypedVCF = "~{chromosome}.genotyped.vcf.gz"
        File genotypedVCF_tbi = "~{chromosome}.genotyped.vcf.gz.tbi"
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
            -O "~{chromosome}.filtered.genotyped.vcf.gz" --select-type-to-include SNP \
            --restrict-alleles-to BIALLELIC \
            -select "AN >= 4 && MQ > 40.0 && QD > 7.0 && DP >= 10.0 && DP <= 1000.0"

#        tabix -p vcf "~{chromosome}.filtered.genotyped.vcf.gz"
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
        File filtered_genotyped_VCF_tbi = "~{chromosome}.filtered.genotyped.vcf.gz.tbi"
    }
}