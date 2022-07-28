version 1.0

workflow collate_data {
    input {
        Array[File] inputFiles
        String filename
        }

        call collate {
            input:
                inputFiles = inputFiles,
                filename = filename
        }

    output {
        File quant_diffs = collate.out
        }
    }

task collate {
    input {
        Array[File] inputFiles
        String filename
        }

    command <<<
        paste ~{sep=' ' inputFiles} | \
        cut -f1,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,104,108,112,116,120,124,128,132,136,140,144,148,152,156,160,164,168,172,176,180,184,188,192,196,200 \
        > "~{filename}"

        #check that the output has 51 columns
        awk '{print NF}' "~{filename}" | sort -nu | tail -n 1 > value.txt
    >>>

    runtime {
        docker : "drpintothe2nd/ac3_xysupp"
        disks : "local-disk 10 HDD"
        memory: "4 GB"
        cpu : 1
    }

    output {
        File out = "~{filename}"
		File value = "value.txt"
    }
}
