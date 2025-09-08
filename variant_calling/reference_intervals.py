import sys

# give the reference index/list of chromosomes with lengths as input
ref_chroms = sys.argv[1]

interval_size = sys.argv[2]

ref_chroms_dictlist = []

with open(ref_chroms, 'r') as f:
    for line in f:
        if line.startswith('#') or line.strip() == '':
            continue
        parts = line.strip().split()
        chrom = parts[0]
        length = int(parts[1])
        ref_chroms_dictlist.append({'chrom': chrom, 'length': length})

with open('reference_intervals.txt', 'w') as out:
    for chrom_info in ref_chroms_dictlist:
        chrom = chrom_info['chrom']
        length = chrom_info['length']
        for start in range(1, length, int(interval_size)):
            end = min(start + int(interval_size) -1, length)
            out.write(f"{chrom}\t{start}\t{end}\n")

