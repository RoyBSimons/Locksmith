#!/usr/bin/python3
import os
import csv
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-f", "--file", dest="filename",
                    help="open FILE", metavar="FILE")
parser.add_argument("-o", "--output", dest="outputname",
                    help="write report to FILE", metavar="OUTPUTFILE")
parser.add_argument("-r", "--range", dest="target_region",
                    help="The length of target region", metavar="range")
parser.add_argument("-p", "--existing_panel_path", dest="existing_panel_csv_path",
                    help="Path to existing panel", metavar="panel")
args = vars(parser.parse_args())

filename=args["filename"]
outputname=args["outputname"]
target_region=int(args["target_region"])
existing_panel_csv_path = args["existing_panel_csv_path"]
existing_panel_cpg_id_list = []
if existing_panel_csv_path == 'False':
    pass
else:
    with open(existing_panel_csv_path, 'r') as handle:
        reader = csv.reader(handle, delimiter = ',')
        next(reader, None)
        for row in reader:
            if len(row) > 2: #if there is a probe in the panel for this locus
                existing_panel_cpg_id_list.append(row[1])
            else:
                pass

outputlist=[]
with open(filename) as inputfile:
    input_reader=csv.reader(inputfile)
    for row in input_reader:
        splitrow=row[0].split("\t")
        chrom=splitrow[0]
        region_start=int(splitrow[1])-target_region
        region_stop=int(splitrow[2])+target_region
        cpg_id=splitrow[3]
        if cpg_id in existing_panel_cpg_id_list:
            pass    # If there already exists a probe for this target do not design a new one.
        else:
            outputlist.append([chrom,region_start,region_stop,cpg_id])

with open(outputname, 'w') as handle:
    outputfile=csv.writer(handle,delimiter="\t")
    for row in outputlist:
        outputfile.writerow(row)
