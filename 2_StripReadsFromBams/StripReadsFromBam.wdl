#written by Brendan Pinto

version 1.0 

#WORKFLOW DEFINITION
workflow StripReadsFromBam {
    input {
    String bam_to_reads_disk_size
    String bam_to_reads_mem_size
	String SampleName
}

#converts BAM to FASTQ (R1 + R2)
call BamToReads {
	input:
	disk_size = bam_to_reads_disk_size,
	mem_size = bam_to_reads_mem_size,
	Name = SampleName
}

#converts FASTQ to paired FASTQ (R1 + R2)
call RepairReads {
	input:
	Name = SampleName,
	inputReads = BamToReads.outputReads
        }

#quality and adapter trims paired reads
call TrimReads {
	input:
	Name = SampleName,
	inputR1 = RepairReads.outputR1,
	inputR2 = RepairReads.outputR2
        }
  
    output {
		File trimR1 = TrimReads.trimR1
		File trimR2 = TrimReads.trimR2
                }
}

#Task Definitions
task BamToReads {
  input {
  File inputBam
  String Name
  String disk_size
  String mem_size
  }
  
#Calls samtools view to do the conversion
command {
#Set -e and -o says if any command I run fails in this script, make sure to return a failure
set -e
set -o pipefail

samtools fastq -c9 -@4 -n -o ${Name}.fq.gz ${inputBam} 

}

#Run time attributes:
runtime {
    docker: "drpintothe2nd/ac3_xysupp"
    memory: mem_size
    cpu: "4"
    disks: "local-disk " + disk_size + " HDD"
	}
    
output {
	File outputReads = "${Name}.fq.gz"
	}
}

#call Task #2: RepairReads (bbmap)
task RepairReads {
    input {
    File inputReads
    String Name
}

command {
#Set -e and -o says if any command I run fails in this script, make sure to return a failure
set -e
set -o pipefail

repair.sh -Xmx40g in=${inputReads} out1=${Name}_R1.fq.gz out2=${Name}_R2.fq.gz outs=${Name}_singletons.fq.gz ziplevel=9

    }

#Run time attributes:
Int diskGb = ceil(10.0 * size(inputReads, "G"))

runtime {
    docker: "drpintothe2nd/ac3_xysupp"
    memory: "40G"
    cpu: "2"
    disks: "local-disk ${diskGb} HDD"
    }
    
output {
	File outputR1 = "${Name}_R1.fq.gz"
	File outputR2 = "${Name}_R2.fq.gz"
	}
}

#call Task #3: TrimReads (Trim Galore!)
task TrimReads {
  input {
  File inputR1
  File inputR2
  String Name
}

command {
#Set -e and -o says if any command I run fails in this script, make sure to return a failure
set -e
set -o pipefail

trim_galore --paired --cores 4 ${inputR1} ${inputR2}

    }

#Run time attributes:
Int diskGb = ceil(10.0 * size(inputR1, "G"))

runtime {
    docker: "drpintothe2nd/ac3_xysupp"
    memory: "8G"
    cpu: "4"
    disks: "local-disk ${diskGb} HDD"
    }
    
#Outputs a BAM and BAI with the same sample name
output {
	File trimR1 = "${Name}_R1_val_1.fq.gz"
	File trimR2 = "${Name}_R2_val_2.fq.gz"
	}
}
