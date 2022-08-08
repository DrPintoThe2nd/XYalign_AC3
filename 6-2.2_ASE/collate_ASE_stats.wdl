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
        File ASE_stats = collate.out
        }
    }

task collate {
    input {
        Array[File] inputFiles
        String filename
        }

    command <<<

        cat ~{sep=' ' inputFiles} | awk ' NR==1 || $1 != "sample_id" ' > "~{filename}"

    >>>

    runtime {
        docker : "drpintothe2nd/ac3_xysupp"
        disks : "local-disk 10 HDD"
        memory: "4 GB"
        cpu : 1
    }

    output {
        File out = "~{filename}"
    }
}
