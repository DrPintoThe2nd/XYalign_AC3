#written by Samantha Zarate
#modified by Brendan Pinto

version 1.0

workflow calculate_vcfstats {
    input {
        File inputVCFgz
        File inputVCFgz_tbi
    }

    call calc_stats {
        input:
            inputVCFgz = inputVCFgz,
            inputVCFgz_tbi = inputVCFgz_tbi
    }

    output {
        File bcfstats = calc_stats.bcfstats
        File rtgstats = calc_stats.rtgstats
    }
}

task calc_stats {
    input {
        File inputVCFgz
        File inputVCFgz_tbi
    }

    String vcfName = '~{basename(inputVCFgz,".vcf.gz")}'
    String vcfFile = '~{basename(inputVCFgz)}'

    command <<<
        cp "~{inputVCFgz}" .
        cp "~{inputVCFgz_tbi}" .

        bcftools stats "~{vcfFile}" > "~{vcfName}.bcftools.stats.txt"
        rtg vcfstats "~{vcfFile}" > "~{vcfName}.rtg.stats.txt"

    >>>

    Int diskGb = ceil(2.0 * size(inputVCFgz, "G"))

    runtime {
        docker : "drpintothe2nd/ac3_xysupp:latest"
        disks : "local-disk ${diskGb} SSD"
        memory: "4G"
        cpu : 1
    }

    output {
        File bcfstats = "~{vcfName}.bcftools.stats.txt"
        File rtgstats = "~{vcfName}.rtg.stats.txt"
    }
}
