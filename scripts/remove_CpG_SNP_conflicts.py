#!/usr/bin/python3
#Compute all possible forward and reverse primers.
import os
from argparse import ArgumentParser
from Bio import SeqIO, Seq
import csv
import operator

parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="filename",
                    help="open FILE", metavar="FILE")
parser.add_argument("-o", "--output", dest="outputname",
                    help="write report to FILE", metavar="OUTPUTFILE")
parser.add_argument("-u", "--input_up", dest="inputname_up",
                    help="write report to FILE", metavar="INPUTFILE_UP")
parser.add_argument("-d", "--intput_down", dest="inputname_down",
                    help="write report to FILE", metavar="INPUTFILE_DOWN")
parser.add_argument("-s", "--input_SNP_up", dest="inputname_SNP_up",
                    help="write report to FILE", metavar="INPUTFILESNP_UP")
parser.add_argument("-t", "--intput_SNP_down", dest="inputname_SNP_down",
                    help="write report to FILE", metavar="INPUTFILESNP_DOWN")
parser.add_argument("-l", "--number_of_iteration", dest="iteration",
                    help="The number of wanted CpGs in arms", metavar="ITERATION")


args = vars(parser.parse_args())

filename=args["filename"]
outputname=args["outputname"]
inputname_up=args["inputname_up"]
inputname_down=args["inputname_down"]
inputname_SNP_up=args["inputname_SNP_up"]
inputname_SNP_down=args["inputname_SNP_down"]
iteration=int(args["iteration"])

cpg_list_up=[]
cpg_list_down=[]
snp_list_up=[]
snp_list_down=[]

#CpG
with open(inputname_up) as handle:
    reader=csv.reader(handle,delimiter="\t")
    for row in reader:
        cpg_list_up.append(int(row[-1]))

with open(inputname_down) as handle:
    reader=csv.reader(handle,delimiter="\t")
    for row in reader:
        cpg_list_down.append(int(row[-1]))

#SNP
with open(inputname_SNP_up) as handle:
    reader=csv.reader(handle,delimiter="\t")
    for row in reader:
        snp_list_up.append(int(row[-1]))

with open(inputname_SNP_down) as handle:
    reader=csv.reader(handle,delimiter="\t")
    for row in reader:
        snp_list_down.append(int(row[-1]))

combined_cpg_snp_list=list(map(operator.add,list(map(operator.add,cpg_list_up,cpg_list_down)),list(map(operator.add,snp_list_up,snp_list_down))))

with open(outputname,"w") as handle:
    outputfile=csv.writer(handle,delimiter="\t")
    with open(filename) as handle:
        header=handle.readline().rstrip('\n').split("\t")
        outputfile.writerow(header)
        reader=csv.reader(handle,delimiter="\t")
        for i,row in enumerate(reader):
            if combined_cpg_snp_list[i] == iteration:
                outputfile.writerow(row)
            else:
                pass
