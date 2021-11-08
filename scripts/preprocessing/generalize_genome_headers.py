#!/usr/bin/python3
#extract high frequent SNPs from target VCF file.
import os
from argparse import ArgumentParser
from Bio import SeqIO, Seq
import csv
parser = ArgumentParser()
parser.add_argument("-g", "--genome", dest="genome_dir",
                    help="Directory containing the fasta files for each chromosome", metavar="DIR")
parser.add_argument("-i", "--header", dest="headerfile",
                    help="file containing the header to chromosome name conversion", metavar="HEADERFILE")

args = vars(parser.parse_args())

genome_dir=args["genome_dir"]
headerfile=args["headerfile"]


with open(headerfile) as handle:
    handle.readline()
    reader=csv.reader(handle,delimiter='\t')
    chrom_list=[]
    acc_list=[]
    for line in reader:
        chrom_list.append(line[0])
        acc_list.append(line[1])

for filename in os.listdir(genome_dir):
    if filename.endswith('.fasta') or filename.endswith('.fa') or filename.endswith('.fna'):
        with open(genome_dir+filename) as handle, open(genome_dir+filename+'_norm','w') as outputhandle:
            for record in SeqIO.parse(handle,"fasta"):
                index=acc_list.index(record.id)
                chrom='chr'+chrom_list[index]
                record.id=chrom
                record.description=chrom
                SeqIO.write(record,outputhandle,'fasta')
