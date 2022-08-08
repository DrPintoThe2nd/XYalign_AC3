# After using bedtools to intersect the cds fasta bad file and the asereadcounter result, there are some duplicated line
import sys
coordinates_set = set()
out = open(sys.argv[2], 'w')
with open(sys.argv[1], 'r') as f:
    for line in f:
        items = line.rstrip('\n').split('\t')
        coordinate = (items[1], items[2])
        if not coordinate in coordinates_set:
            print (line.rstrip('\n'), file=out)
            coordinates_set.add(coordinate)