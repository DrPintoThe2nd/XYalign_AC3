#written by Brendan Pinto

version 1.0

workflow t2t_RNAalignment {
    input {
        File targetRef
        String sampleName
        File inputFastq1
        File inputFastq2
        File hisatIndexTar
    }

        call alignRNA {
            input:
                read1 = inputFastq1,
                read2 = inputFastq2,
                targetRef = targetRef,
                sampleName = sampleName,
                hisatIndexTar = hisatIndexTar
        }

    output {
        File bam = alignRNA.bam
        File bai = alignRNA.bai
    }
}

task alignRNA {
    input {
        File read1
        File read2
        File targetRef
        String sampleName
        File hisatIndexTar
    }

    String bamBase='~{sampleName}'

    String fastaName='~{basename(targetRef, ".fasta.gz")}'
    String fastqName1="~{basename(read1)}"
    String fastqName2="~{basename(read2)}"

    String tarName='~{basename(hisatIndexTar)}'
    String tarDir='~{basename(hisatIndexTar, ".tar.gz")}'


    command <<<

        cp "~{read1}" .
        cp "~{read2}" .
        cp "~{hisatIndexTar}" .

        tar -zxzf "~{tarName}"

        hisat2 -p "$(nproc)" --dta --rg-id ~{sampleName} --rg SM:~{sampleName} --rg LB:~{sampleName} \
        --rg PU:1 --rg PL:illumina -x "./~{tarDir}/~{fastaName}" -1 "./~{fastqName1}" -2 "./~{fastqName2}" \
        | samtools fixmate -@ "$(nproc)" - - | samtools sort -@ "$(nproc)" -O bam - -o "~{bamBase}.bam"

        samtools index "~{bamBase}.bam"
        
    >>>

    runtime {
        docker : "drpintothe2nd/ac3_xysupp"
        disks : "local-disk 200 HDD"
        memory: "24 GB"
        cpu : 4
    }

    output {
        File bam = "~{bamBase}.bam"
        File bai = "~{bamBase}.bam.bai"
    }
}
