# -*- coding: utf-8 -*-
# date : 2022/11/16 
# author : wangh
# file : 2.get_contigN50.py
# project : Unified_Yeast_GEMs_Database_from_13pro
# Get the size of genomes, the number of contigs, N50
from Bio import SeqIO
import os

indir = 'data/genome/1900_assembled_genome/'

outfile1 = 'code_standardize/strains_information_collection/outputs/genomeStats.txt'
outfile2 = 'code_standardize/strains_information_collection/outputs/contigSize.txt'

fhand1 = open(outfile1, 'w')
fhand1.write('#genome,genomeSize,contigNums,N50\n')

fhand2 = open(outfile2, 'w')

def getN50(numlist):
    """
    Abstract: Returns the N50 value of the passed list of numbers.
    Usage:    N50(numlist)
    Based on the Broad Institute definition:
    https://www.broad.harvard.edu/crd/wiki/index.php/N50
    """
    numlist.sort()
    # newlist=numlist
    newlist = []
    for x in numlist :
        newlist += [x]*x
    # take the mean of the two middle elements if there are an even number
    # of elements.  otherwise, take the middle element
    if len(newlist)%2 == 0:
        medianpos = int(len(newlist)/2)
        return float(newlist[medianpos] + newlist[medianpos-1])/2
    else:
        medianpos = int((len(newlist)+1)/2)
        return newlist[medianpos]
    # print(medianpos)

for name in os.listdir(indir):
    if name.startswith('.'): continue
    print(name)
    contig_num = 0
    contig_size = list()
    for rec in SeqIO.parse(os.path.join(indir,name),'fasta'):
        L = len(rec.seq)
        contig_num += 1
        contig_size.append(L)
        fhand2.write('{}\n'.format(L))

    N50 = getN50(numlist=contig_size)

    fhand1.write('{},{},{},{}\n'.format(name,sum(contig_size),contig_num,N50))

fhand1.close()
fhand2.close()