#written by Brendan J. Pinto

version 1.0

workflow ASEreadcounter {
    input {
        String sampleName
        File ASE_VCF
        File ASE_VCF_tbi
        File RNA_bam
        File RNA_bai
        String tissue
        String chromosome
        File refFasta
        File refIndex
        File refDict
    }

    call ASE {
        input:
            sampleName = sampleName,
            ASE_VCF = ASE_VCF,
            ASE_VCF_tbi = ASE_VCF_tbi,
            RNA_bam = RNA_bam,
            RNA_bai = RNA_bai,
            refFasta = refFasta,
            refIndex = refIndex,
            refDict = refDict,
            chromosome = chromosome
    }

    output {
        File ASE_out= ASE.ASE_out
    }
}

task ASE {
    input {
        String sampleName
        File refFasta
        File refIndex
        File refDict
        String chromosome
        String tissue
        File ASE_VCF
        File ASE_VCF_tbi
        File RNA_bam
        File RNA_bai
    }

    String fastaName='~{basename(refFasta)}'
    String inputVCF='~{basename(ASE_VCF)}'
    String inputBam='~{basename(RNA_bam)}'

    command <<<

        cp "~{ASE_VCF}" .
        cp "~{ASE_VCF_tbi}" .
        cp "~{RNA_bam}" .
        cp "~{RNA_bai}" .
        cp "~{refFasta}" .
        cp "~{refIndex}" .
        cp "~{refDict}" .
        
        gatk ASEReadCounter -R "~{fastaName}" -O "~{sampleName}_~{chromosome}_~{tissue}_ASE.tsv" -I "~{inputBam}" -V "~{inputVCF}" \
        -min-depth 10 --min-mapping-quality 10 --min-base-quality 10 

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
        File ASE_out = "~{sampleName}_~{chromosome}_~{tissue}_ASE.tsv"
    }
}
