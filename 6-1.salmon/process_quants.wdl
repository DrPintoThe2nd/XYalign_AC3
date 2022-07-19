#written by Brendan Pinto

version 1.0

workflow process_quants {
    input {
        String sampleName
        File quants_raw
        File quants_xy
        File chr8_transcripts
        File chrX_transcripts
    }

        call split_chrs {
            input:
                sampleName = sampleName,
                quants_raw = quants_raw,
                quants_xy = quants_xy,
                chr8_transcripts = chr8_transcripts,
                chrX_transcripts = chrX_transcripts
        }

    output {
        File chr8_raw_quants = split_chrs.chr8_raw_out
        File chrX_raw_quants = split_chrs.chrX_raw_out
        File chr8_xy_quants = split_chrs.chr8_xy_out
        File chrX_xy_quants = split_chrs.chrX_xy_out
        File chr8_TPM = split_chrs.chr8_TPM_out
        File chrX_TPM = split_chrs.chrX_TPM_out
    }
}

task split_chrs {
    input {
        String sampleName
        File quants_raw
        File quants_xy
        File chr8_transcripts
        File chrX_transcripts
    }

    command <<<

        #extract chr8 transcripts for each indvidual (RAW)
        head -1 "~{quants_raw}" > "~{sampleName}_chr8_quants_raw.tsv"
        grep -f "~{chr8_transcripts}" "~{quants_raw}" | sort >> "~{sampleName}_chr8_quants_raw.tsv"
        #extract chr8 transcripts for each indvidual (XY)
        head -1 "~{quants_xy}" > "~{sampleName}_chr8_quants_xy.tsv"
        grep -f "~{chr8_transcripts}" "~{quants_xy}" | sort >> "~{sampleName}_chr8_quants_xy.tsv"

        #combine TPM values from raw and XYalign for each indvidual and calculate the difference between the two
        cut -f1,4 "~{sampleName}_chr8_quants_xy.tsv" > tmp.tsv
        cut -f4 "~{sampleName}_chr8_quants_raw.tsv" | paste tmp.tsv - > tmp2.tsv
        awk '{print $2-$3}' tmp2.tsv | paste tmp2.tsv - > "~{sampleName}_chr8_quants_xy_raw.tsv" 
        #fix column headers
        sed -i "s/TPM/~{sampleName}_TPM/g" "~{sampleName}_chr8_quants_xy_raw.tsv" 
        sed -i "0,/"[0]"/s//~{sampleName}/" "~{sampleName}_chr8_quants_xy_raw.tsv" 

        #extract chrX transcripts for each indvidual (RAW)
        head -1 "~{quants_raw}" > "~{sampleName}_chrX_quants_raw.tsv"
        grep -f "~{chrX_transcripts}" "~{quants_raw}" | sort >> "~{sampleName}_chrX_quants_raw.tsv"
        #extract chrX transcripts for each indvidual (XY)
        head -1 "~{quants_xy}" > "~{sampleName}_chrX_quants_xy.tsv"
        grep -f "~{chrX_transcripts}" "~{quants_xy}" | sort >> "~{sampleName}_chrX_quants_xy.tsv"

        #combine TPM values from raw and XYalign for each indvidual and calculate the difference between the two
        cut -f1,4 "~{sampleName}_chrX_quants_xy.tsv" > tmp.tsv
        cut -f4 "~{sampleName}_chrX_quants_raw.tsv" | paste tmp.tsv - > tmp2.tsv
        awk '{print $2-$3}' tmp2.tsv | paste tmp2.tsv - > "~{sampleName}_chrX_quants_xy_raw.tsv" 
        #fix column headers
        sed -i "s/TPM/~{sampleName}_TPM/g" "~{sampleName}_chrX_quants_xy_raw.tsv" 
        sed -i "0,/"[0]"/s//~{sampleName}/" "~{sampleName}_chrX_quants_xy_raw.tsv" 

    >>>

    runtime {
        docker : "drpintothe2nd/ac3_xysupp"
        disks : "local-disk 10 HDD"
        memory: "4 GB"
        cpu : 1
    }

    output {
        File chr8_raw_out = "~{sampleName}_chr8_quants_raw.tsv"
        File chrX_raw_out = "~{sampleName}_chrX_quants_raw.tsv"
        File chr8_xy_out = "~{sampleName}_chr8_quants_xy.tsv"
        File chrX_xy_out = "~{sampleName}_chrX_quants_xy.tsv"
        File chr8_TPM_out = "~{sampleName}_chr8_quants_xy_raw.tsv"
        File chrX_TPM_out = "~{sampleName}_chrX_quants_xy_raw.tsv"
    }
}
