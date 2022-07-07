#!bin/bash

gffread GCF_009914755.1_T2T-CHM13v2.0_genomic.gff -g GCA_009914755.4_CHM13_T2T_v2.0_genomic_ChrNamesAdded.fasta -x GCA_009914755.4_CHM13_T2T_v2.0_genomic_ChrNamesAdded_transcripts_RAW.fasta
gffread GCF_009914755.1_T2T-CHM13v2.0_genomic.gff -g GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_ChrNamesAdded.fasta -x GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_ChrNamesAdded_transcripts_XY.fasta

mkdir decoy_RAW decoy_XY salmon_index_RAW salmon_index_XY

bash generateDecoyTranscriptome.sh -b ~/miniconda3/envs/salmon/bin/bedtools -m ~/miniconda3/envs/salmon/bin/mashmap \
-a GCF_009914755.1_T2T-CHM13v2.0_genomic.gff -g GCA_009914755.4_CHM13_T2T_v2.0_genomic_ChrNamesAdded.fasta \
-t GCA_009914755.4_CHM13_T2T_v2.0_genomic_ChrNamesAdded_transcripts_RAW.fasta -o decoy_RAW/ -j 4;

bash generateDecoyTranscriptome.sh -b ~/miniconda3/envs/salmon/bin/bedtools -m ~/miniconda3/envs/salmon/bin/mashmap \
-a GCF_009914755.1_T2T-CHM13v2.0_genomic.gff -g GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_ChrNamesAdded.fasta \
-t GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_ChrNamesAdded_transcripts_XY.fasta -o decoy_XY/ -j 4;

salmon index -t decoy_RAW/gentrome.fa -d decoy_RAW/decoys.txt -i salmon_index_RAW;

salmon index -t decoy_XY/gentrome.fa -d decoy_XY/decoys.txt -i salmon_index_XY;
