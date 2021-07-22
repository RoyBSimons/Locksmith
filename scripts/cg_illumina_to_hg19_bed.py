#!/usr/bin/python3
#creates a bed file from cglist
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-c", "--cgid_list", dest="cglistfile",
                    help="open FILE", metavar="FILE1")
parser.add_argument("-i", "--illumina_cgid_list", dest="illumina_cgid_listfile",
                    help="open FILE", metavar="FILE2")
parser.add_argument("-o", "--output", dest="outputname",
                    help="write report to FILE", metavar="OUTPUTFILE")

args = vars(parser.parse_args())

cglistfile=args["cglistfile"]
illumina_cgid_listfile=args["illumina_cgid_listfile"]
outputname=args["outputname"]
import csv
cglist=[]

with open(cglistfile) as inputfile1:
    input_reader1=csv.reader(inputfile1)
    for row in input_reader1:
        cglist.append(row)

chr_list=[]
loc_list=[]
with open(illumina_cgid_listfile) as inputfile:
    input_reader=csv.reader(inputfile,delimiter=",")
    for row in input_reader:
        if [row[1]] in cglist:
            chr_list.append(str(row[3]))
            loc_list.append(str(row[4]))
        else:
            pass

for i,loc in enumerate(loc_list):
    os.system("grep chr"+str(chr_list[i])+"'\t'"+str(loc)+"'\t' ./data/CpG_bed_nomenclature/CpG_hg19_chr"+str(chr_list[i])+".bed >> "+outputname+" -h")
