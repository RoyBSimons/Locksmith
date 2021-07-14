#!/usr/bin/python3
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-f", "--file", dest="filename",
                    help="open FILE", metavar="FILE")
parser.add_argument("-o", "--output", dest="outputname",
                    help="write report to FILE", metavar="OUTPUTFILE")
parser.add_argument("-r", "--range", dest="target_region",
                    help="The length of target region", metavar="range")
args = vars(parser.parse_args())

filename=args["filename"]
outputname=args["outputname"]
target_region=int(args["target_region"])
import csv
outputlist=[]
with open(filename) as inputfile:
    input_reader=csv.reader(inputfile)
    for row in input_reader:
        splitrow=row[0].split("\t")
        chrom=splitrow[0]
        region_start=int(splitrow[1])-target_region
        region_stop=int(splitrow[2])+target_region
        cpg_id=splitrow[3]
        outputlist.append([chrom,region_start,region_stop,cpg_id])

with open(outputname, 'w') as handle:
    outputfile=csv.writer(handle,delimiter="\t")
    for row in outputlist:
        outputfile.writerow(row)
