#!/usr/bin/python3
#analyzing probe dimer conflicts
import csv
from argparse import ArgumentParser
import json
import matplotlib
import matplotlib.pyplot as plt
parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="input_name",
                    help="input json file (mfeprimer dimer output)", metavar="INPUTNAME")
parser.add_argument("-o", "--output", dest="output_name",
                    help="Fasta file which will contain the chosen probes as output", metavar="OUTPUTNAME")
probe_list=[]
probe_list_count=[]
Tm_list=[]
args=vars(parser.parse_args())
input_name=args['input_name']
output_name=args['output_name']
with open(input_name) as handle:
    data=json.load(handle)
    for row in data:
        Tm=row['Tm']
        probe=row['S1']['ID']
        Tm_list.append(Tm)
        probe2=row['S2']['ID']
        if probe in probe_list:
            probe_list_count[probe_list.index(probe)]+=1
        else:
            probe_list.append(probe)
            probe_list_count.append(1)
        if probe2 in probe_list:
            probe_list_count[probe_list.index(probe2)]+=1
        else:
            probe_list.append(probe2)
            probe_list_count.append(1)
matplotlib.use('Agg')
plt.hist(Tm_list)
plt.xlabel('Tm (C)')
plt.ylabel('Occurence')
plt.savefig(output_name+'.png')
nr_of_times=max(probe_list_count)
print('Total amount of dimer forming probes is ',len(probe_list))
sorted_nr_list=[]
sorted_probe_list=[]
while nr_of_times >0:
    for i,nr in enumerate(probe_list_count):
        if nr == nr_of_times:
            sorted_probe_list.append(probe_list[i])
            sorted_nr_list.append(nr)
    nr_of_times=nr_of_times-1




with open(output_name,'w') as outputfile:
    outputfile.write('Probe_name\tnr\n')
    writer=csv.writer(outputfile,delimiter='\t')
    for i,probe in enumerate(sorted_probe_list):
        writer.writerow([probe,sorted_nr_list[i]])
