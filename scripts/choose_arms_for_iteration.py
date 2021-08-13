#!/usr/bin/python3
#
from argparse import ArgumentParser
from operator import itemgetter
import csv
from itertools import compress
parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="filename",
                    help="open FILE", metavar="FILE")
parser.add_argument("-o", "--output", dest="outputname",
                    help="write report to FILE", metavar="OUTPUTFILE")
parser.add_argument("-n", "--input_not_selected", dest="inputname_not_selected",
                    help="tsv file containing the CpGs for which no arms were designed in the previous iteration", metavar="INPUTFILE_NOT_SELECTED")

args = vars(parser.parse_args())

filename=args["filename"]
outputname=args["outputname"]
inputname_not_selected=args["inputname_not_selected"]

if inputname_not_selected == 'output/Arm_selection_initialization_file':
    with open(outputname, 'w') as handle:
        writer=csv.writer(handle)
        with open(filename) as handle2:
            reader=csv.reader(handle2)
            for row in reader:
                writer.writerow(row)
else:
    cpg_list_non_selected=[]
    with open(inputname_not_selected) as handle:
        reader=csv.reader(handle)
        for row in reader:
            cpg_list_non_selected.append(row[0])
        
    cpg_id_list=[]
    arm_list=[]
    with open(filename) as handle:
        reader=csv.reader(handle,delimiter='\t')
        for row in reader:
            CpG_id=row[-2].split('_')[0]
            cpg_id_list.append(CpG_id)
            arm_list.append(row)
    
    Boolean_list=[cpg in cpg_list_non_selected for cpg in cpg_id_list]
    selected_arms=list(compress(arm_list,Boolean_list))
    
    
    with open(outputname, "w") as handle:
        writer=csv.writer(handle, delimiter="\t")
        for row in selected_arms:
            writer.writerow(row)

