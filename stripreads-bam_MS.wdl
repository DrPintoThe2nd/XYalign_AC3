#RUN WITH:
#cd /mnt/c/Users/brend/Desktop/miniwdl/stripreads-bam
#miniwdl run stripreads-bam.wdl -i stripreads-bam.json --as-me --debug --no-cache

version 1.0 

#WORKFLOW DEFINITION
workflow StripReadsFromBam {
String bam_to_reads_disk_size
String bam_to_reads_mem_size
String SampleName

#converts BAM to FASTQ
call StripReads {
	input:
	disk_size = bam_to_reads_disk_size,
	mem_size = bam_to_reads_mem_size,
	Name = SampleName
        }

#converts FASTQ to paired FASTQ (R1 + R2)
call RepairReads {
	input:
	disk_size = bam_to_reads_disk_size,
	mem_size = bam_to_reads_mem_size,
	Name = SampleName,
	inputReads = StripReads.outputReads
        }

#quality and adapter trims paired reads
call TrimReads {
	input:
	disk_size = bam_to_reads_disk_size,
	mem_size = bam_to_reads_mem_size,
	Name = SampleName,
	inputR1 = RepairReads.outputR1,
	inputR2 = RepairReads.outputR2
        }
  
meta {
    author: "Brendan J. Pinto"
    email: "bjp004@morningside.edu"
    description: "a workflow that strips reads from a bam file, pairs them, and trims them for downstream application"
    }
}

#call Task #1: StripReads (samtools)
task StripReads {
  File InputBam
  String Name
  String disk_size
  String mem_size

command {
#Set -e and -o says if any command I run fails in this script, make sure to return a failure
set -e
set -o pipefail

samtools fastq -c9 -@6 -n -o ${Name}.fq.gz ${InputBam} 

}

#Run time attributes:
#disk_size should equal input size + output size + buffer

runtime {
    docker: "drpintothe2nd/ac3_xysupp"
    memory: mem_size
    cpu: "4"
    disks: "local-disk " + disk_size + " HDD"
    }
    
#Outputs a BAM and BAI with the same sample name
output {
	File outputReads = "${Name}.fq.gz"
	}
}

#call Task #2: RepairReads (bbmap)
task RepairReads {
  File inputReads
  String Name
  String disk_size
  String mem_size

command {
#Set -e and -o says if any command I run fails in this script, make sure to return a failure
set -e
set -o pipefail

repair.sh -Xmx40g in=${inputReads} out1=${Name}_R1.fq.gz out2=${Name}_R2.fq.gz outs=${Name}_singletons.fq.gz ziplevel=9

    }

#Run time attributes:
#disk_size should equal input size + output size + buffer

runtime {
    docker: "drpintothe2nd/ac3_xysupp"
    memory: mem_size
    cpu: "4"
    disks: "local-disk " + disk_size + " HDD"
    }
    
#Outputs a BAM and BAI with the same sample name
output {
	File outputR1 = "${Name}_R1.fq.gz"
	File outputR2 = "${Name}_R2.fq.gz"
	}
}

#call Task #3: TrimReads (Trim Galore!)
task TrimReads {
  File inputR1
  File inputR2
  String Name
  String disk_size
  String mem_size

command {
#Set -e and -o says if any command I run fails in this script, make sure to return a failure
set -e
set -o pipefail

trim_galore --paired --cores 4 ${inputR1} ${inputR2}

    }

#Run time attributes:
#disk_size should equal input size + output size + buffer

runtime {
    docker: "drpintothe2nd/ac3_xysupp"
    memory: mem_size
    cpu: "4"
    disks: "local-disk " + disk_size + " HDD"
    }
    
#Outputs a BAM and BAI with the same sample name
output {
	File TrimR1 = "${Name}_R1_val_1.fq.gz"
	File TrimR2 = "${Name}_R2_val_2.fq.gz"
	}
}
