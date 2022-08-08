# In this script, we want to tabulate:
# 1. Mean allele balance
# 2. Median allele balance
# 1. The number of variants that are in coding regions
# 2. The number of variants that are in coding regions where allele balance is greater than 0.8
# 3. The proportion of variants that are in coding regions where allele balance is greater than 0.8
# 4. The number of variants that are in coding regions where allele balance is equal to 1
# 5. The proportion of variants that are in coding regions where allele balance is equal to 1

import pandas as pd
import sys
import os
import argparse

parser = argparse.ArgumentParser(description='Tabulate different statistics in coding region')
parser.add_argument("--sampleID",required=True,help="Input the sample ID.")
parser.add_argument("--tissue",required=True,help="Input the tissue name")
parser.add_argument("--chrom",required=True,help="Input the chromosome number")
parser.add_argument("--infile",required=True,help="Input the file. For example: /scratch/tphung3/Cayo_x_inactivation/06_analyze_skewed_variants/results/PT105290/PT105290_adrenal_filter_for_hardfilter_biallelic_hets_rmsameRefAlt_chrX_coding_regions_rmdups_allelebalance.bed")
parser.add_argument("--outfile",required=True,help="Input the outfile")

args = parser.parse_args()

sample = args.sampleID #sample id
tissue = args.tissue #tissue
chrom = args.chrom


data = pd.read_csv(args.infile, sep='\t', header=None,
                   names=['chr_cds', 'cds_start', 'cds_end', 'transcript', 'gene', 'chr', 'start', 'end', 'ref', 'alt', 'ref_count', 'alt_count', 'total_count', 'allele_balance']) #example: /scratch/tphung3/Cayo_x_inactivation/06_analyze_skewed_variants/results/PT105290/PT105290_adrenal_filter_for_hardfilter_biallelic_hets_rmsameRefAlt_chrX_coding_regions_rmdups_allelebalance.bed

mean_allele_balance = data['allele_balance'].mean()
median_allele_balance = data['allele_balance'].median()

n_variants_in_coding = data.shape[0]
n_variants_in_coding_skewed = data[data['allele_balance'] > 0.8].shape[0]
prop_variants_in_coding_skewed = n_variants_in_coding_skewed/n_variants_in_coding
n_variants_in_coding_completely_skewed = data[data['allele_balance'] == 1].shape[0]
prop_variants_in_coding_completely_skewed = n_variants_in_coding_completely_skewed/n_variants_in_coding

output = [sample, tissue, chrom, str(mean_allele_balance), str(median_allele_balance), str(n_variants_in_coding), str(n_variants_in_coding_skewed), str(prop_variants_in_coding_skewed), str(n_variants_in_coding_completely_skewed), str(prop_variants_in_coding_completely_skewed)]

if not os.path.exists(args.outfile):
    outfile = open(args.outfile, 'w')
    header = ['sample_id', 'tissue', 'chr', 'mean_allele_balance', 'median_allele_balance', 'num_var_coding', 'num_var_coding_skewed', 'prop_var_coding_skewed', 'num_var_coding_complete_skewed', 'prop_var_coding_complete_skewed']
    print ('\t'.join(header), file=outfile)
    print ('\t'.join(output), file=outfile)
else:
    outfile = open(args.outfile, 'a')
    print ('\t'.join(output), file=outfile)