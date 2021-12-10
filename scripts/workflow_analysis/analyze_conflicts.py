#!/usr/bin/python3
#analyzing probe dimer conflicts
import csv
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="input_name",
                    help="input json file (mfeprimer dimer output)", metavar="INPUTNAME")
parser.add_argument("-o", "--output", dest="output_name",
                    help="Fasta file which will contain the chosen probes as output", metavar="OUTPUTNAME")
args=vars(parser.parse_args())
input_name=args['input_name']
output_name=args['output_name']
iterations=[]
with open(input_name,'r') as handle:
    reader=csv.reader(handle,delimiter='\t')
    for row in reader:
        iterations.append(row)
last_iteration=iterations[-1]
cpg_list=[]
conflicts_per_cpg_list=[]
cpg_list_highest=[]

for conflict in last_iteration:
    if conflict.split(':')[0] in cpg_list:
        index=cpg_list.index(conflict.split(':')[0])
        if conflict in conflicts_per_cpg_list[index]:
            pass
        else:
            conflicts_per_cpg_list[index].append(conflict)
            highest=max(int(conflict.split(':')[1]),int(cpg_list_highest[index]))
            cpg_list_highest[index]=highest
    else:
        cpg_list.append(conflict.split(':')[0])
        conflicts_per_cpg_list.append([conflict])
        highest=int(conflict.split(':')[1])
        cpg_list_highest.append(highest)

for i,cpg in enumerate(cpg_list):
    if len(conflicts_per_cpg_list[i])==cpg_list_highest[i]-1:
        print(str(i)+ ' all cpgs are conflicted')
    else:
        print('Not here: '+str(len(conflicts_per_cpg_list[i]))+' of '+str(cpg_list_highest[i]-1)+' give conflicts')
