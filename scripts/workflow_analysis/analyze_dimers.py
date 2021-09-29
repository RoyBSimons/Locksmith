#!/usr/bin/python3
#analyzing probe dimer conflicts
import csv
probe_list=[]
probe_list_count=[]
with open('output/selection_0/conflicting_probes_dimers.tsv') as handle:
    handle.readline()
    reader=csv.reader(handle,delimiter='\t')
    for row in reader:
        probe=row[0]
        if probe in probe_list:
            probe_list_count[probe_list.index(probe)]+=1
        else:
            probe_list.append(probe)
            probe_list_count.append(1)
nr_of_times=max(probe_list_count)
sorted_nr_list=[]
sorted_probe_list=[]
while nr_of_times >0:
    for i,nr in enumerate(probe_list_count):
        if nr == nr_of_times:
            sorted_probe_list.append(probe_list[i])
            sorted_nr_list.append(nr)
    nr_of_times=nr_of_times-1




with open('Analysis_probe_conflicts.tsv','w') as outputfile:
    outputfile.write('Probe_name\tnr\n')
    writer=csv.writer(outputfile,delimiter='\t')
    for i,probe in enumerate(sorted_probe_list):
        writer.writerow([probe,sorted_nr_list[i]])
