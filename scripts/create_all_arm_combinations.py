#!/usr/bin/python3
#Compute all possible forward and reverse primers.
import os
from argparse import ArgumentParser
import json
from Bio import SeqIO, Seq, SeqRecord
from Bio.SeqUtils import GC
import subprocess
import multiprocessing as mp

parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="filename",
                    help="open FILE", metavar="FILE")
parser.add_argument("-o", "--output", dest="outputname",
                    help="write report to FILE", metavar="OUTPUTFILE")
parser.add_argument("-c", "--config_file", dest="config_file",
                    help="Config file with probe specifics", metavar="JSONfile")
parser.add_argument("-b", "--bed", dest="bedfile",
                    help="Bed file containing the targets", metavar="BED")
parser.add_argument("-t", "--cores", dest="cores",
                    help="Passed amount of cores to script", metavar="Cores")
args = vars(parser.parse_args())

filename=args["filename"]
bedfile=args["bedfile"]
outputname=args["outputname"]
config_file=args["config_file"]
nr_of_cores=int(args["cores"])
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
backbone_sequence=configObject["Backbone_sequence"]
mid_loc=int(target_range)+1 #The middle nucleotide is the target range +1 because the total target is created by adding the target range on both sides.

def create_all_possible_arms(record,min_target_length,max_target_length,min_arm_length,max_arm_length,min_CG_percentage,max_CG_percentage,cpg_flanks,mid_loc):
    probe_list=[]
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
                    probe=[new_u_arm,new_d_arm,upstream_id,downstream_id,i,target_length]
                    i+=1
                    probe_list.append(probe)
    return probe_list

def create_all_possible_arms_both_strands(record,min_target_length,max_target_length,min_arm_length,max_arm_length,min_CG_percentage,max_CG_percentage,cpg_flanks,mid_loc):
    probe_list=[]
    i=0
    rev_record=record.reverse_complement()
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
                    probe=[new_u_arm,new_d_arm,upstream_id,downstream_id,i,target_length,'+']
                    i+=1
                    probe_list.append(probe)
    for target_length in range(min_target_length,max_target_length+1): #loop over range of target lengths, including the maximum
        for arm_length_downstream in range(min_arm_length,max_arm_length+1):
            start_loc_downstream=mid_loc-(target_length-cpg_flanks)
            end_loc_downstream=mid_loc-cpg_flanks
            for start_loc in range(start_loc_downstream,end_loc_downstream):
                start_loc_upstream=start_loc+target_length
                new_d_arm_rev=rev_record.seq[start_loc-arm_length_downstream:start_loc]
                if GC(new_d_arm_rev)> max_CG_percentage or GC(new_d_arm_rev)<min_CG_percentage:
                    break
                else:
                    pass
                upstream_id_rev=record.id.split(':')[0]+":"+str(int(record.id.split(':')[1].split("-")[0])+start_loc-arm_length_downstream)+"-"+str(int(record.id.split(':')[1].split("-")[0])+start_loc-1)
                for arm_length_upstream in range(min_arm_length,max_arm_length+1):
                    new_u_arm_rev=rev_record.seq[start_loc_upstream:start_loc_upstream+arm_length_upstream]
                    if GC(new_u_arm_rev)> max_CG_percentage or GC(new_u_arm_rev)<min_CG_percentage:
                        break
                    else:
                        pass
                    downstream_id_rev=record.id.split(':')[0]+":"+str(int(record.id.split(':')[1].split("-")[0])+start_loc_upstream)+"-"+str(int(record.id.split(':')[1].split("-")[0])+start_loc_upstream+arm_length_upstream-1)
                    probe=[new_u_arm_rev,new_d_arm_rev,upstream_id_rev,downstream_id_rev,i,target_length,'-']
                    i+=1
                    probe_list.append(probe)
    return probe_list
with open(filename) as handle:
    for record in SeqIO.parse(handle,"fasta"):
        pass
#
def get_delta_Tm(probe_arms): #Wallace rule from (Abdul-Latiff et al., 2017)
#    print(probe)
#    print('above this line is the probe')
    upstream_arm=probe_arms[0]
    downstream_arm=probe_arms[1]
#    print(upstream_arm)
    A=upstream_arm.count('A')
    T=upstream_arm.count('T')
    G=upstream_arm.count('G')
    C=upstream_arm.count('C')
    Tm_up=64.9+41*(G+C-16.4)/(A+T+C+G)
    A=downstream_arm.count('A')
    T=downstream_arm.count('T')
    G=downstream_arm.count('G')
    C=downstream_arm.count('C')
    Tm_down=64.9+41*(G+C-16.4)/(A+T+C+G)
    delta_Tm=abs(Tm_up-Tm_down)
    return delta_Tm

def report_CpGs_in_arms(probe_arms):
    upstream_arm=probe_arms[0]
    downstream_arm=probe_arms[1]
    counts=upstream_arm.count('CG')+downstream_arm.count('CG')
    return counts

def add_backbone(probe_arms,backbone_sequence):
    upstream_arm=probe_arms[0]
    downstream_arm=probe_arms[1]
    probe=upstream_arm+backbone_sequence+downstream_arm
    return str(probe)

def check_probe_for_hairpin_score(probe_list,fasta_name,outputname_json):
    seq_list=[]
    with open(fasta_name,'w') as output_handle:
        for i,probe_arms in enumerate(probe_list):
            for j,probe in enumerate(probe_arms):
                seq_list.append(SeqRecord.SeqRecord(Seq.Seq(probe),id=str(i)+'-'+str(j)))
        SeqIO.write(seq_list,output_handle,'fasta')
    os.system('mfeprimer hairpin --in '+fasta_name+' -j -o '+outputname_json) #Remove/move the output files
    with open(outputname_json+'.json') as jsonFile:
        hairpinlist = json.load(jsonFile)
        jsonFile.close()
    os.system('rm '+outputname_json+'.json')
    os.system('rm '+outputname_json)
    os.system('rm '+fasta_name)
    hairpin_scores=[[0 for probe in probe_arms] for probe_arms in probe_list]
    if hairpinlist is None:
        pass
    else:
        for hairpin in hairpinlist:
            cpg_id=hairpin['Seq']['ID']
            score=hairpin['Score']
            i=int(cpg_id.split('-')[0])
            j=int(cpg_id.split('-')[1])
            hairpin_scores[i][j]=score
    return hairpin_scores

#----------------------------------------------------------------------------------------------------------

with open(filename) as handle:
    pool=mp.Pool(nr_of_cores)#mp.cpu_count())
    possible_arm_combinations_all_targets=[pool.apply(create_all_possible_arms_both_strands,args=(record,min_target_length,max_target_length,min_arm_length,max_arm_length,min_target_length,max_CG_percentage,cpg_flanks,mid_loc)) for record in SeqIO.parse(handle,"fasta")]
    pool.close()

pool=mp.Pool(nr_of_cores)
tms=[[pool.apply(get_delta_Tm, args=([probe_arms])) for probe_arms in possible_arm_combinations] for possible_arm_combinations in possible_arm_combinations_all_targets]
pool.close()

pool=mp.Pool(nr_of_cores)
CpG_conflicts=[[pool.apply(report_CpGs_in_arms, args=([probe_arms])) for probe_arms in possible_arm_combinations] for possible_arm_combinations in possible_arm_combinations_all_targets]
pool.close()

pool=mp.Pool(nr_of_cores)
possible_probes_all_targets=[[pool.apply(add_backbone, args=(probe_arms,backbone_sequence)) for probe_arms in possible_arm_combinations] for possible_arm_combinations in possible_arm_combinations_all_targets]
pool.close()

fasta_name='temp_test_fasta.fa'
outputname_json='temp_test_json'
hairpin_scores=check_probe_for_hairpin_score(possible_probes_all_targets,fasta_name,outputname_json)
#to implement
    #SNPs in arms
    #hairpin-score

#test print first probe. First arm combinations from first target site.
print(tms[0][0])
print(possible_arm_combinations_all_targets[0][0])
print(CpG_conflicts[0][0])
print(possible_probes_all_targets[0][0])
print(hairpin_scores)
