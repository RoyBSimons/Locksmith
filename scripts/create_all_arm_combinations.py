#!/usr/bin/python3
#Compute all possible forward and reverse primers.
import os
from argparse import ArgumentParser
import json
from Bio import SeqIO, Seq
from Bio.SeqUtils import GC
import subprocess

parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="filename",
                    help="open FILE", metavar="FILE")
parser.add_argument("-o", "--output", dest="outputname",
                    help="write report to FILE", metavar="OUTPUTFILE")
parser.add_argument("-c", "--config_file", dest="config_file",
                    help="Config file with probe specifics", metavar="JSONfile")
parser.add_argument("-u", "--output_up", dest="outputname_up",
                    help="write report to FILE", metavar="OUTPUTFILE_UP")
parser.add_argument("-d", "--output_down", dest="outputname_down",
                    help="write report to FILE", metavar="OUTPUTFILE_DOWN")
parser.add_argument("-b", "--bed", dest="bedfile",
                    help="Bed file containing the targets", metavar="BED")

args = vars(parser.parse_args())

filename=args["filename"]
bedfile=args["bedfile"]
outputname=args["outputname"]
outputname_up=args["outputname_up"]
outputname_down=args["outputname_down"]
config_file=args["config_file"]
with open(config_file) as jsonFile:
    configObject = json.load(jsonFile)
    jsonFile.close()


import csv
probe_specifics=configObject['probe_specifics'][0]
min_arm_length=probe_specifics['min_arm_length']
max_arm_length=probe_specifics['max_arm_length']
min_target_length=probe_specifics['min_target_length']
max_target_length=probe_specifics['max_target_length']
target_range=configObject['target_range']
cpg_flanks=probe_specifics['cpg_flanks']
max_CG_percentage=float(probe_specifics['max_CG_percentage'])
min_CG_percentage=float(probe_specifics['min_CG_percentage'])

mid_loc=int(target_range)+1 #The middle nucleotide is the target range +1 because the total target is created by adding the target range on both sides.



downstream_arms=[]
upstream_arms=[]
target_length_list=[]
cg_list=[]
downstream_tm=[]
upstream_tm=[]
downstream_ids=[]
upstream_ids=[]
record_list=[]
arm_nr_list=[]

with open(outputname, 'w') as handle:
    outputfile=csv.writer(handle, delimiter="\t")
    with open(filename) as handle:
        for record in SeqIO.parse(handle,"fasta"):
            #reverse Convert the record
            i=0
            target=record.reverse_complement()
            for target_length in range(min_target_length,max_target_length+1): #loop over range of target lengths, including the maximum
                for arm_length_downstream in range(min_arm_length,max_arm_length+1):
                    start_loc_downstream=mid_loc-(target_length-cpg_flanks)
                    end_loc_downstream=mid_loc-cpg_flanks
                    for start_loc in range(start_loc_downstream,end_loc_downstream):
                        start_loc_upstream=start_loc+target_length
                        new_d_arm=record.seq[start_loc-arm_length_downstream:start_loc]
                        if GC(new_d_arm)> max_CG_percentage or GC(new_d_arm)<min_CG_percentage:
                            break
                        else:
                            pass
                        downstream_id=record.id.split(':')[0]+":"+str(int(record.id.split(':')[1].split("-")[0])+start_loc-arm_length_downstream)+"-"+str(int(record.id.split(':')[1].split("-")[0])+start_loc-1)
                        for arm_length_upstream in range(min_arm_length,max_arm_length+1):
                            new_u_arm=record.seq[start_loc_upstream:start_loc_upstream+arm_length_upstream]
                            if GC(new_u_arm)> max_CG_percentage or GC(new_d_arm)<min_CG_percentage:
                                break
                            else:
                                pass
                            upstream_id=record.id.split(':')[0]+":"+str(int(record.id.split(':')[1].split("-")[0])+start_loc_upstream)+"-"+str(int(record.id.split(':')[1].split("-")[0])+start_loc_upstream+arm_length_upstream-1)
                            downstream_arms.append(new_d_arm)
                            downstream_ids.append(downstream_id)
                            upstream_arms.append(new_u_arm)
                            upstream_ids.append(upstream_id)
                            target_length_list.append(target_length)
                            record_list.append(record.id)
                            arm_nr_list.append(i)
                            i+=1
cg_id_list_in=[]
record_id_list_in=[]
with open(bedfile) as handle:
    bedreader=csv.reader(handle,delimiter="\t")
    for row in bedreader:
        cg_id_list_in.append(row[3])
        record_id_list_in.append(row[0]+":"+row[1]+"-"+row[2])
cg_id_list=[]
for record_id in record_list:
    cg_id=cg_id_list_in[record_id_list_in.index(record_id)]
    cg_id_list.append(cg_id)

with open(outputname, 'w') as handle:
    outputfile=csv.writer(handle, delimiter="\t")
    for i,upstream_arm in enumerate(upstream_arms):
        outputfile.writerow([upstream_arm,downstream_arms[i],upstream_ids[i],downstream_ids[i],cg_id_list[i]+"_arm_"+str(arm_nr_list[i]),target_length_list[i]])

with open(outputname_up,"w") as handle:
    outputfile=csv.writer(handle,delimiter="\t")
    for i,row in enumerate(upstream_ids):
        chrom=row.split(":")[0]
        start_id=row.split(":")[1].split("-")[0]
        end_id=row.split(":")[1].split("-")[1]
        outputfile.writerow([chrom,start_id,end_id,cg_id_list[i]+"_arm_"+str(arm_nr_list[i])])
with open(outputname_down,"w") as handle:
    outputfile=csv.writer(handle,delimiter="\t")
    for i,row in enumerate(downstream_ids):
        chrom=row.split(":")[0]
        start_id=row.split(":")[1].split("-")[0]
        end_id=row.split(":")[1].split("-")[1]
        outputfile.writerow([chrom,start_id,end_id,cg_id_list[i]+"_arm_"+str(arm_nr_list[i])])
