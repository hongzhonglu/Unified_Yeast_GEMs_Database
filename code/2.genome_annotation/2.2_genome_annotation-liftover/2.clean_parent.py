# change the head line of genome fasta file

infile = 'data/R64_cleaned.gff'
outfile = 'data/R64_cleaned_parent.gff'

fhand = open(outfile, 'w')
for line in open(infile):
    if line.startswith('#'):
        fhand.write(line)
        continue
    cont = line.split('\t')
    if cont[2] == 'mRNA':
        parent = ';' + cont[-1].strip().split(';')[-1]
        line = line.replace(parent,'')
    fhand.write(line)
fhand.close()
