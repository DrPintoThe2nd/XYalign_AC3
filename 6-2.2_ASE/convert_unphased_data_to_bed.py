# In this script, we are converting the results from the unphased data to bed file format
import argparse

parser = argparse.ArgumentParser(description='Convert unphased data to bed file format.')
parser.add_argument('--unphased_data',required=True,help='Input the path to the unphased data.')
parser.add_argument('--outfile',required=True,help='Input the path to the output file for unphased data in bed file format')

args = parser.parse_args()

outfile = open(args.outfile, 'w')
with open(args.unphased_data, 'r') as f:
    for line in f:
        if not line.startswith('contig'):
            items = line.rstrip('\n').split('\t')
            chr = items[0]
            position = int(items[1])
            out = [chr, str(position), str(position+1), items[3], items[4], items[5], items[6], items[7]]
            print ('\t'.join(out), file=outfile)