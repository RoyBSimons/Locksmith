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

args = vars(parser.parse_args())

filename=args["filename"]
outputname=args["outputname"]
inputname_up=args["inputname_up"]
inputname_down=args["inputname_down"]

cpg_list_up=[]
cpg_list_down=[]
cpg_free_arms=[]
cpg_rows=[]
cpgs=[]

with open(inputname_up) as handle:
    reader=csv.reader(handle,delimiter="\t")
    for row in reader:
        cpg_list_up.append(int(row[-1]))

with open(inputname_down) as handle:
    reader=csv.reader(handle,delimiter="\t")
    for row in reader:
        cpg_list_down.append(int(row[-1]))

combined_cpg_list=list(map(operator.add,cpg_list_up,cpg_list_down))

with open(outputname,"w") as handle:
    outputfile=csv.writer(handle,delimiter="\t")
    with open(filename) as handle:
        reader=csv.reader(handle,delimiter="\t")
        for i,row in enumerate(reader):
            cpg_rows.append(row)
            if combined_cpg_list[i] == 0:
                pass
            else:
                 outputfile.writerow(row)
