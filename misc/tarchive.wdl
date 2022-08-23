#written by Brendan J. Pinto

version 1.0

workflow make_tar {
    input {
        Array[File] inputFiles
        String directory
        }

        call tarchive {
            input:
                inputFiles = inputFiles,
                directory = directory
        }

    output {
        File tar = tarchive.out
        }
    }

task tarchive {
    input {
        Array[File] inputFiles
        String directory
        }

    command <<<

        mkdir "~{directory}"

        cp ~{sep=' ' inputFiles} "~{directory}"/

        tar -czvf "~{directory}.tar.gz" "~{directory}"/

    >>>

    runtime {
        docker : "drpintothe2nd/ac3_xysupp"
        disks : "local-disk 100 HDD"
        memory: "4 GB"
        cpu : 1
    }

    output {
        File out = "~{directory}.tar.gz"
    }
}
