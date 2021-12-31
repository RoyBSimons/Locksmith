#!/usr/bin/python3
#Script to create bed files with CpG IDs for the reference genome
from argparse import ArgumentParser
import csv
from Bio import SeqIO
import re
parser = ArgumentParser()
parser.add_argument("-g", "--genome_files", nargs = '+', dest = "genome_files",
                    help = "All chromosome fasta files of genome", metavar = "GENOMES")
parser.add_argument("-c", "--chr2acc", dest = "chr2acc",
                    help = "file containing the header to chromosome name conversion", metavar = "HEADERFILE")

args = vars(parser.parse_args())

genome_files=args["genome_files"]
chr2acc=args["chr2acc"]

with open(chr2acc) as handle:
    handle.readline()
    reader=csv.reader(handle,delimiter='\t')
    chrom_list=[]
    acc_list=[]
    for line in reader:
        chrom_list.append(line[0])
        acc_list.append(line[1])

start_ID = 1
nr_list = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']

nr_list_files=[]
for file in genome_files:
    nr = file.split('.')[0].split('_')[-1]
    nr_list_files.append(nr)

for nr in nr_list:
    file_index = nr_list_files.index(nr)
    genome_file=genome_files[file_index]
    bed_file=genome_file.split('.')[0]+'.bed'
    CpG_loc_start_list=[]
    with open(bed_file,'w') as outputfile:
        writer=csv.writer(outputfile,delimiter='\t')
        with open(genome_file) as inputfile:
            for record in SeqIO.parse(inputfile,'fasta'):
                pattern=['CG','cg','Cg','cG']
                regex=re.compile('|'.join(pattern))
                CpG_loc_start_list.extend([m.start()+1 for m in regex.finditer(str(record.seq))])
            chrom_index = acc_list.index(record.id)
            chromosome = 'chr' + chrom_list[chrom_index]
        for index,CpG_loc_start in enumerate(CpG_loc_start_list):
            CpG_ID='CpG'+str(start_ID+index)
            CpG_loc_end=CpG_loc_start+1
            writer.writerow([chromosome,CpG_loc_start,CpG_loc_end,CpG_ID])
        start_ID=start_ID+index+1 #continue numbering for next genome file
