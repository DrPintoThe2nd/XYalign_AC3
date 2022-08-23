#written by Samantha Zarate
#edited by Brendan Pinto

version 1.0

workflow t2t_alignment {
    input {
        File targetRef
        String sampleName
        File inputFastq1
        File inputFastq2
        File bwaIndexTar
    }

        call alignLane {
            input:
                read1 = inputFastq1,
                read2 = inputFastq2,
                targetRef = targetRef,
                sampleName = sampleName,
                bwaIndexTar = bwaIndexTar
        }

    call bamtoCram {
        input:
            targetRef = targetRef,
            inputBam = alignLane.bam,
            sampleName = sampleName
    }

    call samtoolsIndex as indexCRAM {
        input:
            alignmentFile = bamtoCram.cram,
            isBAM = false
    }

    call samtoolsStats {
        input:
            inputCram = bamtoCram.cram,
            cramIndex =  indexCRAM.alignmentIndex,
            targetRef = targetRef,
            sampleName = sampleName
    }

    output {
        File cram = bamtoCram.cram
        File cramIndex = indexCRAM.alignmentIndex
        File samtools_stats = samtoolsStats.stats
    }
}

task alignLane {
    input {
        File read1
        File read2
        File targetRef
        String sampleName
        File bwaIndexTar
    }

    String bamBase='~{sampleName}'

    String fastaName='~{basename(targetRef)}'
    String fastqName1="~{basename(read1)}"
    String fastqName2="~{basename(read2)}"

    String tarName='~{basename(bwaIndexTar)}'
    String tarDir='~{basename(bwaIndexTar, ".tar.gz")}'


    command <<<
        fastqName="~{read1}"

        cp "~{read1}" .
        cp "~{read2}" .
        cp "~{bwaIndexTar}" .

        tar -zxzf "~{tarName}"

        bwa mem -Y \
            -K 100000000 \
            -t "$(nproc)" \
            -R "@RG\tID:~{sampleName}\tPL:illumina\tPM:Unknown\tLB:~{sampleName}\tDS:CHM13\tSM:~{sampleName}\tCN:AC3_ASU\tPU:1" \
            "./~{tarDir}/~{fastaName}" \
            "./~{fastqName1}" \
            "./~{fastqName2}" | samblaster | samtools fixmate -@ 4 - - | samtools sort -@ 4 -O bam - -o "~{bamBase}.bam"
    >>>

    runtime {
        docker : "drpintothe2nd/ac3_xysupp"
        disks : "local-disk 1000 HDD"
        memory: "64 GB"
        cpu : 24
    }

    output {
        File bam = "~{bamBase}.bam"
    }
}

task samtoolsIndex {
    input {
        File alignmentFile
        Boolean isBAM
    }

    String alignmentName = basename(alignmentFile)

    String outputSuffix = if (isBAM) then "bai" else "crai"

    command <<<
        samtools index "~{alignmentFile}" "~{alignmentName}.~{outputSuffix}"
    >>>

    Int diskGb = ceil(2.0 * size(alignmentFile, "G"))

    runtime {
        docker : "drpintothe2nd/ac3_xysupp"
        disks : "local-disk " + "${diskGb}" + " HDD"
        memory: "8 GB"
        cpu : 1
        preemptible: 1
        maxRetries: 1
    }

    output {
        File alignmentIndex = "~{alignmentName}.~{outputSuffix}"
    }
}

task bamtoCram {
    input {
        File targetRef
        File inputBam
        String sampleName
    }

    String fastaName='~{basename(targetRef)}'

    command <<<
        cp "~{targetRef}" .

        samtools view -C -T ./"~{fastaName}" -o "~{sampleName}.cram" "~{inputBam}"
    
	>>>

    Int diskGb = ceil(2.0 * size(inputBam, "G"))

    runtime {
        docker : "drpintothe2nd/ac3_xysupp"
        disks : "local-disk " + "${diskGb}" + " HDD"
        memory: "8 GB"
        cpu : 1
        preemptible: 1
        maxRetries: 1
    }

    output {
        File cram = "~{sampleName}.cram"
    }
}

task samtoolsStats {
    input {
        File inputCram
        File cramIndex
        File targetRef
        String sampleName
    }

    String fastaName='~{basename(targetRef)}'

    command <<<
        cp "~{targetRef}" .

        samtools stats --reference ./"~{fastaName}" "~{inputCram}" > "~{sampleName}.samtools.stats.txt"
    >>>

    Int diskGb = ceil(2.0 * size(inputCram, "G"))

    runtime {
        docker : "drpintothe2nd/ac3_xysupp"
        disks : "local-disk " + 200 + " HDD"
        memory: "8 GB"
        cpu : 1
        preemptible: 1
        maxRetries: 1
    }

    output {
        File stats = "~{sampleName}.samtools.stats.txt"
    }
}
