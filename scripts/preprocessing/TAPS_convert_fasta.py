#!/usr/bin/python3
#extract high frequent SNPs from target VCF file.
import os
from argparse import ArgumentParser
from Bio import SeqIO, Seq

parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="filename",
                    help="open FILE", metavar="FILE")
parser.add_argument("-o", "--output", dest="outputname",
                    help="write report to FILE", metavar="OUTPUTFILE")
parser.add_argument("-c", "--cpgbedfile", dest="cpgbedfile",
                    help="write report to FILE", metavar="OUTPUTFILE")

args = vars(parser.parse_args())

filename=args["filename"]
outputname=args["outputname"]
cpgbedfile=args["cpgbedfile"])
import csv

#Create a list of all CpGs locations in the cpgbedfile
fasta_cpg_loc=[]
fasta_target_cpg_loc=[]
fasta_heads_cpg=[]
with open(cpgbedfile) as inputfile:
    input_reader=csv.reader(inputfile,delimiter="\t")
    for i,row in enumerate(input_reader):
        if i%2 == 1: #only do this for the second part of each row in the bedfile. (corresponds to all CpGs in the n-length target sequence around the target CpG)
            fasta_cpg_loc.append(row[3])
        else: #only for the first part of each row in the bedfile (corresponds to target CpGs)
            fasta_heads_cpg.append(row[0]+":"+row[1]+"-"+row[2])
            fasta_target_cpg_loc.append(row[1])

with open(filename) as handle:
    for record in SeqIO.parse(handle,"fasta"):
        indices=[i for i,cpg in enumerate(fasta_heads_cpg) if cpg == record.id]
        start_nt=fasta_target_cpg_loc[indices[0]]
        nt_cpg_list=fasta_cpg_loc[indices[0]:indices[-1]]
        for i,end_nt in enumerate(nt_cpg_list):
            if i == 0:
                new_seq=record.seq[:int(end_nt)-int(start_nt)-1]
            else:
                new_seq+=Seq.Seq("T")+record.seq[previous_end_nt:int(end_nt)-int(start_nt)-1] #for each CpG the C will turn to a T after conversion
            previous_end_nt=int(end_nt)-int(start_nt)
        new_seq+=Seq.Seq("T")+record.seq[previous_end_nt:] #Add last stretch of record.seq
        Converted_sequence=new_seq
#How to finish this script:
# The converted_sequence can be stored in a list.
# This list can be written in a Fasta file.


#Reason why not to finish
#Conversion of the sequences is actually not needed as no CpGs will be incorporated in both arms of the probe
