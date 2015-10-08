# input custom txt file and extract AGCT motif and user-defined flanking seq

import sys, os, re, argparse, itertools

def parse_args():
    parser = argparse.ArgumentParser(description='Import txt file and extract AGCT and user-defined flanking seq')
    #parser.add_argument("-f", dest = "file", type = str, required = True,
    #                          help = "txt file" )
    parser.add_argument("-d", dest = "dist", type = int, default = 15,
                              help = "flanking seq size" )
    args = parser.parse_args()
    return args

def sequenced_bases():
    '''
    Return allele sequences and sequenced range
    '''
    file_list = ['CDR_constr_2_reference.fas', 'Gpt_reference.fas', 'VB18_passenger_reference.fas', \
                'VH12_reference.fas', 'betaglobin_reference.fas']
    name_list = ['CDRconstr2', 'gpt', 'VB18', 'VH12', 'betaglobin']
    allele_seq = {}   # key: allele name; value: allele seq
    # Load full sequence
    for i in range(0, 5):
        fn = '/Users/zhoudu/Projects/Joyce/motif/lib/sequences/' + file_list[i]
        tmpseq = ''
        for line in open(fn):
            if line.startswith('>'): continue
            else: tmpseq += line.strip()
        allele_seq[name_list[i]] = tmpseq.upper()
    return allele_seq

def mut_load():
    '''
    Return mutation rate for each base in different alleles
    '''
    name_list = ['betaglobin', 'CDRconstr2', 'gpt', 'VH12', 'VB18']
    file_list = ['betaglobin_3-10_stats', 'CDR_constr_2_3-10_stats', 'gpt_3-10_stats', \
                'VH12_3-10_stats', 'VB18_3-10_stats']
    allele_mut = {}
    for i in range(0, 5):
        fn = '../data/%s.txt' % file_list[i]
        name = name_list[i]
        allele_mut[name] = {}
        for line in open(fn):
            if line.startswith('Pos'): continue
            l = line.strip().split()
            if l[-1] == 'NA': continue
            allele_mut[name][int(l[0])] = float(l[-1])
    return allele_mut

def kmer_search(seq, j = 9):
    kmer_num = {}
    for k in range(2, j):
        # enumerate all kmer and store in dict
        bases=['A','T','G','C']
        for kmer in [''.join(p) for p in itertools.product(bases, repeat=k)]:
            kmer_num[kmer] = 0
        # compute how many kmer in each seq
        matches = re.finditer(r'(?=(\w{%d}))' % k, seq)
        for m in [match.group(1) for match in matches]:
            kmer_num[m] += 1
    return kmer_num

def kmer_total(seq_output):
    kmerAll = {}
    kmerEx = {}
    for seq in seq_output:
        kmer_num = kmer_search(seq)
        for kmer in kmer_num:
            if kmer not in kmerAll:
                kmerAll[kmer] = kmer_num[kmer]
            else:
                kmerAll[kmer] += kmer_num[kmer]
    for kmer in kmerAll:
        if kmerAll[kmer] > 0:
            kmerEx[kmer] = kmerAll[kmer] 
    return kmerEx

def parseFILE(args):
    # sequenced range
    flanknum = int(args.dist)
    allele_range = {'betaglobin': [26, 385], 'CDRconstr2': [141, 500], 'gpt': [26, 385], 'VH12': [59, 428], 'VB18': [141, 500]}
    allele_seq = sequenced_bases()
    allele_mut = mut_load()
    seq_output = {}
    for allele in allele_seq:
        for match in re.finditer('AGCT', allele_seq[allele]):
            AGCTsite = match.start()
            if (AGCTsite + 1) in allele_mut[allele] and allele_range[allele][0] <= AGCTsite + 1 <= allele_range[allele][1]:
                flankseq = allele_seq[allele][(AGCTsite - flanknum): (AGCTsite + flanknum + 5)]
                seq_output[flankseq] = [allele, AGCTsite + 2, allele_mut[allele][AGCTsite + 2]] 
    kmerEx = kmer_total(seq_output)
    # print anlaysis results
    print 'allele\tposition\tMutationRate\tSequence\t%s' % '\t'.join([kmer for kmer in kmerEx])
    for seq in seq_output:
        kmer_num = kmer_search(seq)
        print '%s\t%d\t%f\t%s\t%s' % (seq_output[seq][0], seq_output[seq][1], seq_output[seq][2], seq, \
                            '\t'.join(str(float(kmer_num[kmer])/kmerEx[kmer]) for kmer in kmerEx))

def main():
    args = parse_args()
    seq_output = parseFILE(args)

main()
