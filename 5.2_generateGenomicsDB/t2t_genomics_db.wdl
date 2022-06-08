#written by Samantha Zarate
#edited by Brendan Pinto

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

    command <<<

        cp "~{sample_map}" .

        gatk --java-options "-Xmx24g -Xms16g" GenomicsDBImport \
            --sample-name-map "~{chromosome}_map.sample_map.tsv" \
            --genomicsdb-workspace-path "~{chromosome}_GenomicsDB" \
            --reader-threads $(nproc) \
            -L "~{chromosome}" \
            --batch-size 50
        
        tar -zvcf "~{chromosome}_genomicsDB.tar.gz" "~{chromosome}_GenomicsDB"
    >>>

    runtime {
        docker : "drpintothe2nd/ac3_xysupp"
        disks : "local-disk 100 SSD"
        memory: "24G"
        cpu : 1
    }

    output {
        File genomicsDB = "~{chromosome}_genomicsDB.tar.gz"
    }
}