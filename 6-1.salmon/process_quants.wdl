#written by Brendan Pinto

version 1.0

workflow process_quants {
    input {
        String sampleName
        File quants
        File chr8_transcripts
        File chrX_transcripts
    }

        call split_chrs {
            input:
                sampleName = sampleName,
                quants = quants,
                chr8 = chr8_transcripts,
                chrX = chrX_transcripts
        }

    output {
        File chr8_quants = split_chrs.chr8_out
        File chrX_quants = split_chrs.chrX_out
    }
}

task split_chrs {
    input {
        String sampleName
        File quants
        File chr8
        File chrX
    }

    command <<<
        head -1 "~{quants}" > chr8_quants.tsv
        grep -f "~{chr8}" "~{quants}" | sort >> chr8_quants.tsv
        
        head -1 "~{quants}" > chrX_quants.tsv
        grep -f "~{chrX}" "~{quants}" | sort >> chrX_quants.tsv

    >>>

    runtime {
        docker : "drpintothe2nd/ac3_xysupp"
        disks : "local-disk 10 HDD"
        memory: "4 GB"
        cpu : 1
    }

    output {
        File chr8_out = "chr8_quants.tsv"
        File chrX_out = "chrX_quants.tsv"
    }
}
