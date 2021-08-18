#!/usr/bin/python3
#Compute all possible forward and reverse primers.
import os
from argparse import ArgumentParser
import json
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="filename",
                    help="open FILE", metavar="FILE")
parser.add_argument("-o", "--output", dest="outputname",
                    help="write report to FILE", metavar="OUTPUTFILE")
parser.add_argument("-c", "--backbone_sequence", dest="Backbone_sequence",
                    help="A string containing the backbone sequence", metavar="BACKBONE_SEQUENCE")

args = vars(parser.parse_args())

filename=args["filename"]
outputname=args["outputname"]
Backbone_sequence=args["Backbone_sequence"]

probe_list=[]
with open(filename) as handle:
    handle.readline().rstrip('\n').split('\t')
    reader=csv.reader(handle,delimiter='\t')
    for row in reader:
        upstream_arm=row[1]
        downstream_arm=row[0]
        id_name=row[4]
        probe_sequence=downstream_arm+Backbone_sequence+upstream_arm
        record=SeqRecord(Seq(probe_sequence),id=id_name)
        probe_list.append(record)
SeqIO.write(probe_list, outputname,'fasta')

