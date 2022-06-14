#written by Samantha Zarate
#modified by Brendan Pinto

version 1.0

workflow t2t_genomics_db {
    input {
        File sample_map
        String chromosome
    }

    call generateGenomicsDB {
        input:
            sample_map = sample_map,
            chromosome = chromosome,
    }

    output {
        File genomicsDBtar = generateGenomicsDB.genomicsDB
    }
}

task generateGenomicsDB {
    input {
        File sample_map
        String chromosome
    }

    String genomicsDB = '~{basename(sample_map, "_map.sample_map.tsv")}'

    command <<<

        cp "~{sample_map}" .

        gatk --java-options "-Xmx24g -Xms16g" GenomicsDBImport \
            --sample-name-map "~{genomicsDB}_map.sample_map.tsv" \
            --genomicsdb-workspace-path "~{genomicsDB}_GenomicsDB" \
            --reader-threads $(nproc) \
            -L "~{chromosome}" \
            --batch-size 50
        
        tar -zvcf "~{genomicsDB}_genomicsDB.tar.gz" "~{genomicsDB}_GenomicsDB"
    >>>

    runtime {
        docker : "drpintothe2nd/ac3_xysupp"
        disks : "local-disk 100 SSD"
        memory: "24G"
        cpu : 2
    }

    output {
        File genomicsDB = "~{chromosome}_genomicsDB.tar.gz"
    }
}