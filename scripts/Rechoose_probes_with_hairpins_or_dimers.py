#!/usr/bin/python3
#Form a list of targets to for which a new probe must be designed due to hairpin or dimer formation.
from argparse import ArgumentParser
import json
import csv
parser = ArgumentParser()
parser.add_argument("-c", "--conflicting_cpg_file", dest="conflicting_cpg_file",
                    help="This is a text file containing the cg_id+arms_id of a probe that forms a hairpin or dimer", metavar="CONFLICTING_CPG_FILE")
parser.add_argument("-o", "--output", dest="outputname",
                    help="write report to FILE", metavar="OUTPUTFILE")
parser.add_argument("-F", "--files", nargs="+", type=str, required=True,
                    help="whitespace delimited list of files...")

args = vars(parser.parse_args())
conflicting_cpg_file=args["conflicting_cpg_file"]
#Read the file with all of the CpG_id's that form hairpins or dimers
conflicting_cpg_list=[]
with open(conflicting_cpg_file) as handle:
    reader=csv.reader(handle)
    for row in reader:
        conflicting_cpg_list.append(row)
#This is raw code for opening multiple files with the -F flag 
anc=[]
for path in args.files:
    with open(path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if "i>" in line:
                anc.append(line)


