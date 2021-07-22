#!/usr/bin/python3
#Compute all possible forward and reverse primers.
import os
from argparse import ArgumentParser
import json
from Bio import SeqIO, Seq
import subprocess

parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="filename",
                    help="open FILE", metavar="FILE")
parser.add_argument("-o", "--output", dest="outputname",
                    help="write report to FILE", metavar="OUTPUTFILE")
parser.add_argument("-c", "--config_file", dest="config_file",
                    help="Config file with probe specifics", metavar="JSONfile")

args = vars(parser.parse_args())

filename=args["filename"]
outputname=args["outputname"]
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


mid_loc=int(target_range)+1 #The middle nucleotide is the target range +1 because the total target is created by adding the target range on both sides.



downstream_arms=[]
upstream_arms=[]
target_length_list=[]
cg_list=[]
downstream_tm=[]
upstream_tm=[]
downstream_ids=[]
upstream_ids=[]

with open(outputname, 'w') as handle:
    outputfile=csv.writer(handle, delimiter="\t")
    with open(filename) as handle:
        for record in SeqIO.parse(handle,"fasta"):
            #reverse Convert the record
            target=record.reverse_complement()
            for target_length in range(min_target_length,max_target_length+1): #loop over range of target lengths, including the maximum
                for arm_length_downstream in range(min_arm_length,max_arm_length+1):
                    start_loc_downstream=mid_loc-(target_length-2*cpg_flanks)
                    end_loc_downstream=mid_loc-cpg_flanks-arm_length_downstream
                    for start_loc in range(start_loc_downstream,end_loc_downstream):
                        start_loc_upstream=start_loc+target_length+arm_length_downstream
                        for arm_length_upstream in range(min_arm_length,max_arm_length+1):
                            new_d_arm=record.seq[start_loc:start_loc+arm_length_downstream]
                            downstream_arms.append(new_d_arm)
                            downstream_id=record.id.split(':')[0]+":"+str(int(record.id.split(':')[1].split("-")[0])+start_loc)+"-"+str(int(record.id.split(':')[1].split("-")[0])+start_loc+arm_length_downstream-1)
                            downstream_ids.append(downstream_id)
                            #proc=subprocess.Popen("~/opt/primer3/src/oligotm "+str(new_d_arm),shell=True, stdout=subprocess.PIPE)
                            #oligo_tm_d=proc.communicate()[0]
                            #downstream_tm.append(oligo_tm_d)
                            new_u_arm=record.seq[start_loc_upstream:start_loc_upstream+arm_length_upstream]
                            upstream_arms.append(new_u_arm)
                            upstream_id=record.id.split(':')[0]+":"+str(int(record.id.split(':')[1].split("-")[0])+start_loc_upstream)+"-"+str(int(record.id.split(':')[1].split("-")[0])+start_loc_upstream+arm_length_upstream-1)
                            upstream_ids.append(upstream_id)
                            target_length_list.append(target_length)
                            cg_list.append(record.id)
                            #proc=subprocess.Popen("~/opt/primer3/src/oligotm "+str(new_u_arm),shell=True, stdout=subprocess.PIPE)
                            #oligo_tm_u=proc.communicate()[0]
                            #upstream_tm.append(oligo_tm_u)
                            outputfile.writerow([new_u_arm,new_d_arm,upstream_id,downstream_id,target_length])
                            #outputfile.writerow([new_u_arm,new_d_arm,oligo_tm_u,oligo_tm_d,record.id,target_length]) 
#with open(outputname, 'w') as handle:
#    outputfile=csv.writer(handle, delimiter="\t")
#    for i,row in enumerate(upstream_arms):
#        outputfile.writerow([row,downstream_arms[i],upstream_tm[i],downstream_tm[i],cg_list[i],target_length_list[i]])

with open(outputname[:-4]+"_upstream_arms.bed","w") as handle:
    outputfile=csv.writer(handle,delimiter="\t")
    for i,row in enumerate(upstream_ids):
        chrom=row.split(":")[0]
        start_id=row.split(":")[1].split("-")[0]
        end_id=row.split(":")[1].split("-")[1]
        outputfile.writerow([chrom,start_id,end_id])
with open(outputname[:-4]+"_downstream_arms.bed","w") as handle:
    outputfile=csv.writer(handle,delimiter="\t")
    for i,row in enumerate(downstream_ids):
        chrom=row.split(":")[0]
        start_id=row.split(":")[1].split("-")[0]
        end_id=row.split(":")[1].split("-")[1]
        outputfile.writerow([chrom,start_id,end_id])

#return fasta 


