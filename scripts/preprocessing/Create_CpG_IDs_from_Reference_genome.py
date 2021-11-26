#!/usr/bin/python3
#Script to create bed files with CpG IDs for the reference genome
#from argparse import ArgumentParser
import csv
from Bio import SeqIO
import re
genome_files=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
start_ID=1

for file in genome_files:
    genome_file='chr'+file+'.fna'
    bed_file='chr'+file+'.bed'
    chromosome='chr'+file
    CpG_loc_start_list=[]
    with open(bed_file,'w') as outputfile:
        writer=csv.writer(outputfile,delimiter='\t')
        with open(genome_file) as inputfile:
            for record in SeqIO.parse(inputfile,'fasta'):
                print(len(record.seq))
                pattern=['CG','cg','Cg','cG']
                regex=re.compile('|'.join(pattern))
#                CpG_loc_start_list.extend([m.start() for m in re.finditer('CG',str(record.seq))])
                CpG_loc_start_list.extend([m.start()+1 for m in regex.finditer(str(record.seq))])
        for index,CpG_loc_start in enumerate(CpG_loc_start_list):
            CpG_ID='CpG'+str(start_ID+index)
            CpG_loc_end=CpG_loc_start+1
            writer.writerow([chromosome,CpG_loc_start,CpG_loc_end,CpG_ID])
        start_ID=start_ID+index+1 #continue numbering for next genome file

