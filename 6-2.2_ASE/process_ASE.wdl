#written by Brendan J. Pinto

version 1.0

workflow process_ASE {
    input {
        File convert_unphased_data_to_bed
        File remove_duplicates
        File compute_allele_balance
        File calc_variants_in_coding_regions
        File ASE_tsv
        File transcripts
        String sampleName
        String tissue
        String chromosome
    }

    call convert_ASE {
        input:
            convert_unphased_data_to_bed = convert_unphased_data_to_bed,
            remove_duplicates = remove_duplicates,
            compute_allele_balance = compute_allele_balance,
            calc_variants_in_coding_regions = calc_variants_in_coding_regions,
            ASE_tsv = ASE_tsv,
            transcripts = transcripts,
            sampleName = sampleName,
            tissue = tissue,
            chromosome = chromosome
    }

    output {
        File ASE_bed = convert_ASE.ASE_bed
        File ASE_transcripts = convert_ASE.ASE_transcripts
        File ASE_rmdups = convert_ASE.ASE_rmdups
        File ASE_allele_balance = convert_ASE.ASE_allele_balance
        File ASE_stats = convert_ASE.ASE_stats
    }
}

task convert_ASE {
    input {
        File convert_unphased_data_to_bed
        File remove_duplicates
        File compute_allele_balance
        File calc_variants_in_coding_regions
        File transcripts
        File ASE_tsv
        String sampleName
        String tissue
        String chromosome
    }

    command <<<

        cp "~{ASE_tsv}" .

        python "~{convert_unphased_data_to_bed}" --unphased_data "~{sampleName}_~{chromosome}_~{tissue}_ASE.tsv" --outfile "~{sampleName}_~{chromosome}_~{tissue}_ASE.bed"

        bedtools intersect -a "~{transcripts}" -b "~{sampleName}_~{chromosome}_~{tissue}_ASE.bed" -wa -wb > "~{sampleName}_~{chromosome}_~{tissue}_ASE_transcripts.bed"

        python "~{remove_duplicates}" "~{sampleName}_~{chromosome}_~{tissue}_ASE_transcripts.bed" "~{sampleName}_~{chromosome}_~{tissue}_ASE_transcripts_rmdup.bed"

        python "~{compute_allele_balance}" "~{sampleName}_~{chromosome}_~{tissue}_ASE_transcripts_rmdup.bed" "~{sampleName}_~{chromosome}_~{tissue}_ASE_transcripts_rmdup_allelebalance.bed"

        python "~{calc_variants_in_coding_regions}"  --sampleID "~{sampleName}" --tissue "~{tissue}" --chrom "~{chromosome}" \
        --infile "~{sampleName}_~{chromosome}_~{tissue}_ASE_transcripts_rmdup_allelebalance.bed" \
        --outfile "~{sampleName}_~{chromosome}_~{tissue}_ASE_transcripts_stats.tsv"

    >>>

    runtime {
        docker : "drpintothe2nd/ac3_xysupp"
        disks : "local-disk 25 HDD"
        memory: "12G"
        cpu : 1
        preemptible: 1
        maxRetries: 1
    }

    output {
        File ASE_bed = "~{sampleName}_~{chromosome}_~{tissue}_ASE.bed"
        File ASE_transcripts = "~{sampleName}_~{chromosome}_~{tissue}_ASE_transcripts.bed"
        File ASE_rmdups = "~{sampleName}_~{chromosome}_~{tissue}_ASE_transcripts_rmdup.bed"
        File ASE_allele_balance = "~{sampleName}_~{chromosome}_~{tissue}_ASE_transcripts_rmdup_allelebalance.bed"
        File ASE_stats = "~{sampleName}_~{chromosome}_~{tissue}_ASE_transcripts_stats.tsv"
    }
}
