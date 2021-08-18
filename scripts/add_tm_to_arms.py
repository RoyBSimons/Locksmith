#!/usr/bin/python3
#Compute all possible forward and reverse primers.
import os
from argparse import ArgumentParser
import csv

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

tm_up=[]
tm_down=[]

with open(inputname_up) as handle:
    reader=csv.reader(handle)
    for row in reader:
        tm_up.append(row[0])

with open(inputname_down) as handle:
    reader=csv.reader(handle)
    for row in reader:
        tm_down.append(row[0])

with open(outputname, "w") as handle_out:
    outputfile=csv.writer(handle_out,delimiter="\t")
    with open(filename) as handle:
        header=handle.readline().rstrip('\n').split('\t')
        header.append('tm_down')
        header.append('tm_up')
        outputfile.writerow(header)
        reader=csv.reader(handle,delimiter="\t")
        for i,row in enumerate(reader):
            row.append(tm_down[i])
            row.append(tm_up[i])
            outputfile.writerow(row)
