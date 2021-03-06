#written by Brendan Pinto

version 1.0

workflow salmon {
    input {
        String sampleName
        File inputFastq1
        File inputFastq2
        File salmonIndexTar
    }

        call salmonQuant {
            input:
                read1 = inputFastq1,
                read2 = inputFastq2,
                sampleName = sampleName,
                salmonIndexTar = salmonIndexTar
        }

    output {
        File salmon = salmonQuant.quant
        File log = salmonQuant.log
    }
}

task salmonQuant {
    input {
        File read1
        File read2
        String sampleName
        File salmonIndexTar
    }

#WDL adaptation
#salmon_index_{RAW/XY}.tar.gz

    String fastqName1="~{basename(read1)}"
    String fastqName2="~{basename(read2)}"

    String tarName='~{basename(salmonIndexTar)}'
    String tarDir='~{basename(salmonIndexTar, ".tar.gz")}'

    command <<<
        fastqName="~{read1}"

        cp "~{read1}" .
        cp "~{read2}" .
        cp "~{salmonIndexTar}" .

        tar -zxzf "~{tarName}"

        salmon quant -i "./~{tarDir}" -l IU -p "$(nproc)" -1 "./~{fastqName1}" -2 "./~{fastqName2}" --gcBias --validateMappings -o "~{sampleName}"

    >>>

    runtime {
        docker : "drpintothe2nd/ac3_xysupp"
        disks : "local-disk 50 HDD"
        memory: "12 GB"
        cpu : 4
    }

    output {
        File quant = "~{sampleName}/quant.sf"
        File log = "~{sampleName}/logs/salmon_quant.log"
    }
}
