#!/usr/bin/python3
#extract high frequent SNPs from target VCF file.
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="filename",
                    help="open FILE", metavar="FILE")
parser.add_argument("-o", "--output", dest="outputname",
                    help="write report to FILE", metavar="OUTPUTFILE")
parser.add_argument("-f", "--frequency_threshold", dest="freq_threshold",
                    help="write report to FILE", metavar="OUTPUTFILE")

args = vars(parser.parse_args())

filename=args["filename"]
outputname=args["outputname"]
freq_threshold=1.0-float(args["freq_threshold"])
import csv
chr_list=[]
loc_list=[]
with open(filename) as inputfile:
    input_reader=csv.reader(inputfile)
    for row in input_reader:
        if "FREQ" in row[0]:
            freq=float(row[0].split(";")[-1].split(":")[-1])
            if freq < freq_threshold:
                os.system("grep " + row[0].split("\t")[2] + " " + filename  + " >>" +outputname)
