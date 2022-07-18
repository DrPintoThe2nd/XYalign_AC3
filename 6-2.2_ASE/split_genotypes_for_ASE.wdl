#written by Brendan J. Pinto

version 1.0

workflow VCF_split_filter {
    input {
        String sampleName
        File filtered_genotyped_VCF
        File filtered_genotyped_VCF_tbi
        String chromosome
        File refFasta
        File refIndex
        File refDict
    }

    call splitVCF {
        input:
            sampleName = sampleName,
            filtered_genotyped_VCF = filtered_genotyped_VCF,
            filtered_genotyped_VCF_tbi = filtered_genotyped_VCF_tbi,
            chromosome = chromosome
    }

    call filterVCF {
        input:
            sampleName = sampleName,
            indiv_VCF = splitVCF.indiv_VCF,
            indiv_VCF_tbi = splitVCF.indiv_VCF_tbi,
            refFasta = refFasta,
            refIndex = refIndex,
            refDict = refDict,
            chromosome = chromosome
    }

    output {
        File ASE_VCF = filterVCF.ASE_VCF
        File ASE_VCF_tbi = filterVCF.ASE_VCF_tbi
    }
}

task splitVCF {
    input {
        String sampleName
        File filtered_genotyped_VCF
        File filtered_genotyped_VCF_tbi
        String chromosome
    }

        String inputVCF='~{basename(filtered_genotyped_VCF)}'

    command <<<

        cp "~{filtered_genotyped_VCF}" .
        cp "~{filtered_genotyped_VCF_tbi}" .

        bcftools view -Oz -s "~{sampleName}" "~{inputVCF}" > "~{sampleName}_~{chromosome}.vcf.gz";
        tabix -p vcf "~{sampleName}_~{chromosome}.vcf.gz";

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
        File indiv_VCF = "~{sampleName}_~{chromosome}.vcf.gz"
        File indiv_VCF_tbi = "~{sampleName}_~{chromosome}.vcf.gz.tbi"
    }
}

#After subsetting for each individual. In some individuals, the genotypes could be homozygous for the reference. This next task is to remove these sites.

task filterVCF {
    input {
        String sampleName
        File refFasta
        File refIndex
        File refDict
        String chromosome
        File indiv_VCF
        File indiv_VCF_tbi
    }

    String fastaName='~{basename(refFasta)}'
    String inputVCF='~{basename(indiv_VCF)}'

    command <<<

        cp "~{indiv_VCF}" .
        cp "~{indiv_VCF_tbi}" .
        cp "~{refFasta}" .
        cp "~{refIndex}" .
        cp "~{refDict}" .

        gatk SelectVariants -R "~{fastaName}" -V "~{inputVCF}" -O "~{sampleName}_~{chromosome}_ASE.vcf.gz" -select "AC == 1";

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
        File ASE_VCF = "~{sampleName}_~{chromosome}_ASE.vcf.gz"
        File ASE_VCF_tbi = "~{sampleName}_~{chromosome}_ASE.vcf.gz.tbi"
    }
}
