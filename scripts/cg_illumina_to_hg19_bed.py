#!/usr/bin/python3
#creates a bed file from cglist
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-f", "--file", dest="filename",
                    help="open FILE", metavar="FILE")
parser.add_argument("-o", "--output", dest="outputname",
                    help="write report to FILE", metavar="OUTPUTFILE")

args = vars(parser.parse_args())

filename=args["filename"]
outputname=args["outputname"]
import csv
chr_list=[]
loc_list=[]
with open(filename) as inputfile:
    input_reader=csv.reader(inputfile)
    for row in input_reader:
        loc_list.append(int(row[0][2:]))
print(loc_list)
for i,loc in enumerate(loc_list):
    os.system("grep  CpG"+str(loc)+"$  ~/genomes/hg19/CpG_bed/*.bed >> "+outputname+" -h")
