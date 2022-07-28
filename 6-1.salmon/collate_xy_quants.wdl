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
        File xy_quants = collate.out
        }
    }

task collate {
    input {
        Array[File] inputFiles
        String filename
        }

    command <<<
        paste ~{sep=' ' inputFiles} | \
        cut -f1,2,6,10,14,18,22,26,30,34,38,42,46,50,54,58,62,66,70,74,78,82,86,90,94,98,102,106,110,114,118,122,126,130,134,138,142,146,150,154,158,162,166,170,174,178,182,186,190,194,198 \
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