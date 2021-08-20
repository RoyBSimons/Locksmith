#!/usr/bin/python3
#Form a list of targets to for which a new probe must be designed due to hairpin or dimer formation.
from argparse import ArgumentParser
import json
import csv
parser = ArgumentParser()
parser.add_argument("-p", "--hairpin_json_file", dest="hairpin_file",
                    help="open FILE", metavar="HAIRPIN_JSON_FILE")
parser.add_argument("-d", "--dimer_json_file", dest="dimer_file",
                    help="open FILE", metavar="HAIRPIN_JSON_FILE")
parser.add_argument("-o", "--output", dest="outputname",
                    help="write report to FILE", metavar="OUTPUTFILE")

args = vars(parser.parse_args())

hairpin_file=args["hairpin_file"]
dimer_file=args["dimer_file"]
outputname=args["outputname"]
with open(hairpin_file) as jsonFile:
    hairpinlist = json.load(jsonFile)
    jsonFile.close()
with open(dimer_file) as jsonFile:
    dimerlist = json.load(jsonFile)
    jsonFile.close()
conflicting_cpg_list1=[]
conflicting_cpg_list2=[]
if hairpinlist is None:
    pass
else:
    for hairpin in hairpinlist:
        cpg_id=hairpin['Seq']['ID']
        conflicting_cpg_list1.append(cpg_id)
if dimerlist is None:
    pass
else:
    for dimer in dimerlist:  
        cpg_id=dimer['S1']['ID']
        cpg_id2=dimer['S2']['ID']
        conflicting_cpg_list1.append(cpg_id)
        conflicting_cpg_list2.append(cpg_id2)
conflicting_cpg_list1=list(dict.fromkeys(conflicting_cpg_list1))
#conflicting_cpg_list2=list(dict.fromkeys(conflicting_cpg_list2))

#Return a file with all of the CpG_id's that form hairpins or dimers
#only return the first mentioned probe in a dimer
with open(outputname, 'w') as handle:
    writer=csv.writer(handle)
    for cpg_id in conflicting_cpg_list1:
        writer.writerow([cpg_id])

