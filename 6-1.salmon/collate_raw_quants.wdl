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
        File raw_quants = collate.out
        }
    }

task collate {
    input {
        Array[File] inputFiles
        String filename
        }

    command <<<
        paste ~{sep=' ' inputFiles} | \
        cut -f1,3,7,11,15,19,23,27,31,35,39,43,47,51,55,59,63,67,71,75,79,83,87,91,95,99,103,107,111,115,119,123,127,131,135,139,143,147,151,155,159,163,167,171,175,179,183,187,191,195,199 \
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