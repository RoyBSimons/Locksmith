#!/usr/bin/python3
#Create panels and choose the best
import os
from argparse import ArgumentParser
import json
from Bio import SeqIO, Seq, SeqRecord
import multiprocessing as mp
from numpy.random import choice
import csv
import numpy as np
import pickle
parser = ArgumentParser()
parser.add_argument("-o", "--output", dest="output_name",
                    help="Fasta file which will contain the chosen probes as output", metavar="OUTPUTNAME")
parser.add_argument("-c", "--config_file", dest="config_file",
                    help="Config file with probe specifics", metavar="JSONfile")
parser.add_argument("-m", "--tms", dest="tms",
                    help="csv file containing probe binding temperatures", metavar="TMS")
parser.add_argument("-g", "--cpg", dest="cpg",
                    help="csv file containing amount of cpgs in probe arms", metavar="CPG")
parser.add_argument("-p", "--probes", dest="probes",
                    help="csv file containing the full probes", metavar="PROBES")
parser.add_argument("-a", "--hairpins", dest="hairpins",
                    help="csv file containing hairpins of probes", metavar="HAIRPINS")
parser.add_argument("-n", "--snp", dest="snps",
                    help="csv file containing amount of snps in probe arms", metavar="SNPS")
parser.add_argument("-r", "--targets", dest="targets",
                    help="Bedfile containing the targets to capture", metavar="BED")
parser.add_argument("-b", "--arms", dest="probe_arms_file",
                    help="csv file containing probe arms", metavar="ARMS")
parser.add_argument("-t", "--cores", dest="cores",
                    help="Passed amount of cores to script", metavar="Cores")
args = vars(parser.parse_args())

output_name=args["output_name"]
config_file=args["config_file"]
tms_file=args['tms']
cpgs=args['cpg']
probes=args['probes']
hairpins=args['hairpins']
snps=args['snps']
targets=args['targets']
probe_arms_file=args['probe_arms_file']
nr_of_cores=int(args["cores"])
with open(config_file) as jsonFile:
    configObject = json.load(jsonFile)
    jsonFile.close()
permutations=int(configObject['Permutations'])
backbone_length=len(configObject["Backbone_sequence"][0]["Reverse_complement_universal_forward_primer"])+len(configObject["Backbone_sequence"][0]["Common_region"])+len(configObject["Backbone_sequence"][0]["Universal_reverse_primer"])

with open(targets,'r') as handle:
    reader=csv.reader(handle,delimiter='\t')
    probe_cpg_id_list=[]
    for row in reader:
        probe_cpg_id_list.append(row[3])

with open(probes, 'r') as handle:
    reader=csv.reader(handle)
    possible_probes_all_targets=[]
    for i,row in enumerate(reader):
        possible_probes_all_targets.append([])
        for probe in row:
            possible_probes_all_targets[i].append(probe)

with open(tms_file, 'rb') as handle:
    tms=pickle.load(handle)
    tms=np.array(tms,dtype=object)
with open(cpgs, 'r') as handle:
    reader=csv.reader(handle)
    CpG_conflicts=[]
    for i,row in enumerate(reader):
        cpg_conflict_list=[]
        for probe in row:
            cpg_conflict_list.append(int(probe))
        CpG_conflicts.append(np.array(cpg_conflict_list,dtype=object))
    CpG_conflicts=np.array(CpG_conflicts,dtype=object)

with open(snps, 'r') as handle:
    reader=csv.reader(handle)
    SNP_conflicts=[]
    for i,row in enumerate(reader):
        snp_conflict_list=[]
        for probe in row:
            snp_conflict_list.append(int(probe))
        SNP_conflicts.append(np.array(snp_conflict_list,dtype=object))
    SNP_conflicts=np.array(SNP_conflicts,dtype=object)

with open(hairpins, 'r') as handle:
    reader=csv.reader(handle)
    hairpin_scores=[]
    for i,row in enumerate(reader):
        hairpin_score_list=[]
        for probe in row:
            hairpin_score_list.append(int(probe))
        hairpin_scores.append(np.array(hairpin_score_list,dtype=object))
    hairpin_scores=np.array(hairpin_scores,dtype=object)
with open(probe_arms_file,'r') as handle:
    reader=csv.reader(handle,delimiter='\t')
    probe_arm_list=[]
    arm_upstream_loc_list=[]
    arm_downstream_loc_list=[]
    probe_id_list=[]
    arm_upstream_list=[]
    arm_downstream_list=[]
    for i,row in enumerate(reader):
        probe_arm_list.append([])
        arm_upstream_loc_list.append([])
        arm_downstream_loc_list.append([])
        probe_id_list.append([])
        arm_upstream_list.append([])
        arm_downstream_list.append([])
        for probe in row:
            probe=probe.strip('][').split(',')
            probe_arm_list[i].append(probe)
            arm_upstream_loc_list[i].append(probe[2])
            arm_downstream_loc_list[i].append(probe[3])
            probe_id_list[i].append(probe[7][2:-1]+':'+probe[4][1:])
            arm_upstream_list[i].append(probe[0])
            arm_downstream_list[i].append(probe[1])

def get_probe_scores_array(tms,CpG_conflicts,SNP_conflicts,hairpin_array):
    probe_score=np.add(np.add(np.add(np.multiply(hairpin_array,int(10^10)),CpG_conflicts),SNP_conflicts),tms)
    return probe_score

probe_scores=get_probe_scores_array(tms,CpG_conflicts,SNP_conflicts,hairpin_scores)
print('\tProbe scores computed')
def choose_probes_from_scores(probe_scores,possible_arm_combinations,n,counter,probe_id_list,seed):
    np.random.seed(seed)
    if probe_scores.size==0:
        return
    else:
        sum_score=sum(probe_scores)
        if sum_score==0: #If all scores are 0, there should still be a probability to choose one of the probes
            probability_distribution=[1/len(probe_scores) for i in probe_scores]
        else:
            probability_distribution=[i/sum_score for i in probe_scores]
        probe= choice(possible_arm_combinations,n,p=probability_distribution)
        index=possible_arm_combinations.index(probe[counter])
    return [probe[counter],probe_id_list[index]]

def find_dimer_forming_probes(chosen_probes):
    #Create a fasta file with the chosen probes
    passed_list=[]
    with open('tmp_fasta_chosen_probes.fasta','w') as handle:
        for i,probe in enumerate(chosen_probes):
            if probe ==None:
                passed_list.append(i)
            else:
                handle.write('>'+str(i)+'\n')
                handle.write(probe+'\n')
    #test the chosen probes set for dimers
    os.system('mfeprimer dimer -i tmp_fasta_chosen_probes.fasta -j -o tmp_dimers_chosen_probes')
    os.system('rm tmp_fasta_chosen_probes.fasta tmp_dimers_chosen_probes')
    dimers=[False for i in range(len(chosen_probes))]
    with open('tmp_dimers_chosen_probes.json') as handle:
        reader = json.load(handle)
        jsonFile.close()
        if reader==None:
            pass
        else:
            for row in reader:
                index1=int(row['S1']['ID'])
                dimers[index1]=True
                index2=int(row['S2']['ID'])
                dimers[index2]=True
    os.system('rm tmp_dimers_chosen_probes.json')
    #obtain which probes form dimers --> No dimer=False, dimer=True
    return dimers
###################################
###################################
###################################
def find_dimer_forming_probes_iterative(chosen_probes):
    #Create a fasta file with the chosen probes
    passed_list=[]
    with open('tmp_fasta_chosen_probes.fasta','w') as handle:
        for i,probe in enumerate(chosen_probes):
            if probe ==None:
                passed_list.append(i)
            else:
                handle.write('>'+probe[1]+'\n')
                handle.write(probe[0]+'\n')
    #test the chosen probes set for dimers
    os.system('cp tmp_fasta_chosen_probes.fasta ~/opt/test.fasta')
    os.system('mfeprimer dimer -i tmp_fasta_chosen_probes.fasta -j -o tmp_dimers_chosen_probes')
    os.system('rm tmp_fasta_chosen_probes.fasta tmp_dimers_chosen_probes')
    dimer_file='tmp_dimers_chosen_probes.json'
    with open(dimer_file) as jsonFile:
        dimerlist = json.load(jsonFile)
        jsonFile.close()
    conflicting_cpg_list_dimer=[]
    conflicting_cpg_list_dimer2=[]
    conflict_range_dimer=[]
    conflict_range_dimer2=[]
    if dimerlist is None:
        pass
    else:
        for dimer in dimerlist:
                cpg_id=dimer['S1']['ID']
                cpg_id2=dimer['S2']['ID']
                dimer_structure=dimer['Aseq']
                S1_length=len(dimer['S1']['Seq'])
                S2_length=len(dimer['S2']['Seq'])
                #Obtain the location of the dimer on the two sequences
                if dimer['S1Dangling']== 0:
                    dimer_start_num=dimer['AseqDangling']
                    dimer_end_num=dimer_start_num+len(dimer_structure)
                    dimer_start_num2=dimer['AseqDangling']-dimer['S2Dangling']
                    dimer_end_num2=dimer_start_num2+len(dimer_structure)
                else:
                    dimer_start_num=dimer['AseqDangling']-dimer['S1Dangling']
                    dimer_end_num=dimer_start_num+len(dimer_structure)
                    dimer_start_num2=dimer['S1Dangling']+dimer['AseqDangling']-dimer['S1Dangling']
                    dimer_end_num2=dimer_start_num2+len(dimer_structure)
                conflicting_cpg_list_dimer.append(cpg_id)
                conflicting_cpg_list_dimer2.append(cpg_id2)
                probe_index=probe_cpg_id_list.index(cpg_id.split(':')[0])
                nested_index=probe_id_list[probe_index].index(cpg_id)
                probe_index2=probe_cpg_id_list.index(cpg_id2.split(':')[0])
                nested_index2=probe_id_list[probe_index2].index(cpg_id2)
                S1_upstream_arm_length=int(arm_upstream_loc_list[probe_index][nested_index].split(':')[1].split('-')[1][:-1])-int(arm_upstream_loc_list[probe_index][nested_index].split(':')[1].split('-')[0])+1
                S1_downstream_arm_length=int(arm_downstream_loc_list[probe_index][nested_index].split(':')[1].split('-')[1][:-1])-int(arm_downstream_loc_list[probe_index][nested_index].split(':')[1].split('-')[0])+1
                S2_upstream_arm_length=int(arm_upstream_loc_list[probe_index2][nested_index2].split(':')[1].split('-')[1][:-1])-int(arm_upstream_loc_list[probe_index2][nested_index2].split(':')[1].split('-')[0])+1
                S2_downstream_arm_length=int(arm_downstream_loc_list[probe_index2][nested_index2].split(':')[1].split('-')[1][:-1])-int(arm_downstream_loc_list[probe_index2][nested_index2].split(':')[1].split('-')[0])+1
                #check whether the S1 and S2 sequences are 5' to 3' or the other way around.
        #        print('\n')
        #        print(dimer['S1']['Seq'][:S1_upstream_arm_length])
        #        print(arm_upstream_loc_list[probe_id_list.index(cpg_id)])
                if arm_upstream_list[probe_index][nested_index] in dimer['S1']['Seq'][:S1_upstream_arm_length]:
        #            print('yes')
                    dimer_start_num_corrected=dimer_start_num
                    dimer_end_num_corrected=dimer_end_num

                else: #In this situation the S1 sequence is 3'-5', therefore we need to correct the indices.
                    dimer_start_num_corrected=S1_downstream_arm_length-(S1_length-dimer_start_num)
                    dimer_end_num_corrected=S1_downstream_arm_length-(S1_length-dimer_end_num)
                if arm_upstream_list[probe_index2][nested_index2] in dimer['S2']['Seq'][:S2_upstream_arm_length]:
        #            print('yes2')
                    dimer_start_num2_corrected=dimer_start_num2
                    dimer_end_num2_corrected=dimer_end_num2
                else: #In this situation the S2 sequence is 3'-5', therefore we need to correct the indices.
                    dimer_start_num2_corrected=S2_length-dimer_end_num2
                    dimer_end_num2_corrected=S2_length-dimer_start_num2

        #The next two if statements only include the situation if dimer['S1Dangling']== 0: (or the other way around?)
        #By this I mean that dimer_range should always give the arm at the start of the designed probe and dimer_range2 should always give the arm at the end of the probe. Although this is only the case if dimer['S1Dangling']== 0 AND both the S1 and S2 sequences are 5'-3'. Look into this.
                #obtain the genomic locations of the dimer forming regions for S1
                if arm_upstream_list[probe_index][nested_index] in dimer['S1']['Seq'][:S1_upstream_arm_length]: #Check whether the dimer was from 3'-5'
        #            print('yes3')
                    if dimer_start_num_corrected < backbone_length: #dimer is in upstream arm
                        chrom=arm_upstream_loc_list[probe_index][nested_index].split(':')[0]
                        start=int(arm_upstream_loc_list[probe_index][nested_index].split(':')[1].split('-')[0].strip("'"))
                        end=int(arm_upstream_loc_list[probe_index][nested_index].split(':')[1].split('-')[1].strip("'"))
                        dimer_start_range=max(start,start+dimer_start_num_corrected)
                        dimer_end_range=min(end,start+dimer_end_num_corrected)
                        dimer_range=str(chrom)+':'+str(dimer_start_range)+'-'+str(dimer_end_range)
                    else:   #dimer is in downstream arm
                        chrom=arm_downstream_loc_list[probe_index][nested_index].split(':')[0]
                        start=int(arm_downstream_loc_list[probe_index][nested_index].split(':')[1].split('-')[0]).strip("'")
                        end=int(arm_downstream_loc_list[probe_index][nested_index].split(':')[1].split('-')[1].strip("'"))
                        dimer_start_range=max(start,start+dimer_start_num_corrected-(S1_length-(end-start+1)))
                        dimer_end_range=min(end,start+dimer_end_num_corrected-(S1_length-(end-start+1)))
                        dimer_range=str(chrom)+':'+str(dimer_start_range)+'-'+str(dimer_end_range)
                else:
                    if dimer_start_num_corrected < backbone_length: #dimer is in upstream arm
                        chrom=arm_downstream_loc_list[probe_index][nested_index].split(':')[0]
                        start=int(arm_downstream_loc_list[probe_index][nested_index].split(':')[1].split('-')[0].strip("'"))
                        end=int(arm_downstream_loc_list[probe_index][nested_index].split(':')[1].split('-')[1].strip("'"))
                        dimer_start_range=max(start,start+dimer_start_num_corrected)
                        dimer_end_range=min(end,start+dimer_end_num_corrected)
                        dimer_range=str(chrom)+':'+str(dimer_start_range)+'-'+str(dimer_end_range)
                    else: #dimer is in downstream arm
                        chrom=arm_downstream_loc_list[probe_index][nested_index].split(':')[0]
                        start=int(arm_downstream_loc_list[probe_index][nested_index].split(':')[1].split('-')[0].strip("'"))
                        end=int(arm_downstream_loc_list[probe_index][nested_index].split(':')[1].split('-')[1].strip("'"))
                        dimer_start_range=max(start,start+dimer_start_num_corrected-(S1_length-(end-start+1)))
                        dimer_end_range=min(end,start+dimer_end_num_corrected-(S1_length-(end-start+1)))
                        dimer_range=str(chrom)+':'+str(dimer_start_range)+'-'+str(dimer_end_range)
            #obtain the genomic locations of the dimer forming regions for S2
                if arm_upstream_list[probe_index2][nested_index2] in dimer['S2']['Seq'][:S2_upstream_arm_length]: #Check whether the dimer was from 3'-5'
        #            print('yes4')
                    if dimer_start_num2_corrected < backbone_length: #dimer is in upstream arm
                        chrom=arm_upstream_loc_list[probe_index2][nested_index2].split(':')[0]
                        start=int(arm_upstream_loc_list[probe_index2][nested_index2].split(':')[1].split('-')[0].strip("'"))
                        end=int(arm_upstream_loc_list[probe_index2][nested_index2].split(':')[1].split('-')[1].strip("'"))
                        dimer_start_range=max(start,start+dimer_start_num2_corrected)
                        dimer_end_range=min(end,start+dimer_end_num2_corrected)
                        dimer_range2=str(chrom)+':'+str(dimer_start_range)+'-'+str(dimer_end_range)
                    else: #dimer is in downstream arm
                        chrom=arm_downstream_loc_list[probe_index2][nested_index2].split(':')[0]
                        start=int(arm_downstream_loc_list[probe_index2][nested_index2].split(':')[1].split('-')[0].strip("'"))
                        end=int(arm_downstream_loc_list[probe_index2][nested_index2].split(':')[1].split('-')[1].strip("'"))
                        dimer_start_range=max(start,start+dimer_start_num2_corrected-(S2_length-(end-start+1)))
                        dimer_end_range=min(end,start+dimer_end_num2_corrected-(S2_length-(end-start+1)))
                        dimer_range2=str(chrom)+':'+str(dimer_start_range)+'-'+str(dimer_end_range)
                else: #If the sequence was obtained in 3'-5' from the dimer file, use downstream_loc_list
                    if dimer_start_num2_corrected < backbone_length: #dimer is in upstream arm
                        chrom=arm_upstream_loc_list[probe_index2][nested_index2].split(':')[0]
                        start=int(arm_upstream_loc_list[probe_index2][nested_index2].split(':')[1].split('-')[0].strip("'"))
                        end=int(arm_upstream_loc_list[probe_index2][nested_index2].split(':')[1].split('-')[1].strip("'"))
                        dimer_start_range=max(start,start+dimer_start_num2_corrected)
                        dimer_end_range=min(end,start+dimer_end_num2_corrected)
                        dimer_range2=str(chrom)+':'+str(dimer_start_range)+'-'+str(dimer_end_range)
                    else: #dimer is in downstream arm
                        chrom=arm_downstream_loc_list[probe_index2][nested_index2].split(':')[0]
                        start=int(arm_downstream_loc_list[probe_index2][nested_index2].split(':')[1].split('-')[0].strip("'"))
                        end=int(arm_downstream_loc_list[probe_index2][nested_index2].split(':')[1].split('-')[1].strip("'"))
                        dimer_start_range=max(start,start+dimer_start_num2_corrected-(S2_length-(end-start+1)))
                        dimer_end_range=min(end,start+dimer_end_num2_corrected-(S2_length-(end-start+1)))
                        dimer_range2=str(chrom)+':'+str(dimer_start_range)+'-'+str(dimer_end_range)
                conflict_range_dimer.append(dimer_range)
                conflict_range_dimer2.append(dimer_range2)

    #Only report the probe that formed the most amount of dimers
    probe_list=[]
    probe_list_count=[]
    for probe in conflicting_cpg_list_dimer: #only look into the first mentioned probes when forming a dimer
        if probe in probe_list:
            probe_list_count[probe_list.index(probe)]+=1
        else:
            probe_list.append(probe)
            probe_list_count.append(1)
    if probe_list_count==[]:
        probe_list_count=[0]
    nr_of_times=max(probe_list_count)
    sorted_nr_list=[]
    sorted_probe_list=[]
    while nr_of_times >0:
        for i,nr in enumerate(probe_list_count):
            if nr == nr_of_times:
                sorted_probe_list.append(probe_list[i])
                sorted_nr_list.append(nr)
        nr_of_times=nr_of_times-1
    most_conflicting_probe_list=sorted_probe_list[0:int(len(sorted_probe_list)/2)] #exclude the top half most dimer-forming proibes
    most_conflicting_dimer_range_list=[]
    for probe in most_conflicting_probe_list: #report the conflicting locus at first occurrence of dimer
        index=conflicting_cpg_list_dimer.index(probe)
        most_conflicting_dimer_range_list.append(conflict_range_dimer[index])

    #find all probes that have common regions with the most_conflicting_dimer_range and report those as conflicting probes
    all_conflicting_probe_list=[]
    conflicting_targets=set([])
    for dimer_range in most_conflicting_dimer_range_list:
        chrom=dimer_range.split(':')[0]
        for i,row in enumerate(arm_downstream_loc_list):
            for j,loc_d in enumerate(row):
                chrom_loc_d=loc_d.split(':')[0]
                if chrom_loc_d==chrom:
                    start=int(loc_d.split(':')[1].split('-')[0].strip("'"))
                    end=int(loc_d.split(':')[1].split('-')[1].strip("'"))
                    dimer_start=int(dimer_range.split(':')[1].split('-')[0].strip("'"))
                    dimer_end=int(dimer_range.split(':')[1].split('-')[1].strip("'"))
                    #when there is more than 6 overlap report the probe as conflicting
                    if start>=dimer_start and start<=dimer_end:
                        overlap=dimer_end-start+1
                    elif end>=dimer_start and end<=dimer_end:
                        overlap=end-dimer_start+1
                    else: #there is no overlap between the two
                        overlap=0
                    if overlap>6: #get this from config file
                        all_conflicting_probe_list.append(probe_id_list[i][j])
                        conflicting_targets.add(probe_id_list[i][0].split(':')[0])
        for i,row in enumerate(arm_upstream_loc_list):
            for j,loc_u in enumerate(row):
                chrom_loc_u=loc_u.split(':')[0]
                if chrom_loc_u==chrom:
                    start=int(loc_u.split(':')[1].split('-')[0].strip("'"))
                    end=int(loc_u.split(':')[1].split('-')[1].strip("'"))
                    dimer_start=int(dimer_range.split(':')[1].split('-')[0].strip("'"))
                    dimer_end=int(dimer_range.split(':')[1].split('-')[1].strip("'"))
                    #when there is more than 6 overlap report the probe as conflicting
                    if start>=dimer_start and start<=dimer_end:
                        overlap=dimer_end-start+1
                    elif end>=dimer_start and end<=dimer_end:
                        overlap=end-dimer_start+1
                    else: #there is no overlap between the two
                        overlap=0
                    if overlap>6: #get this from config file
                        all_conflicting_probe_list.append(probe_id_list[i][j])
                        conflicting_targets.add(probe_id_list[i][0].split(':')[0])
    return all_conflicting_probe_list,conflicting_targets


##############################
##############################
##############################

#-----------------------------------------------------------------------
#Iteratively exclude the most dimer forming probes.
chosen_probes_lists=[]
probes_with_dimers_lists=[]
counter=0
seed=11 #move to config file
min_dimers=len(possible_probes_all_targets)
while counter<permutations and min_dimers>0:
    #Choose random probe set
    # put the below lines in a function that can be called ~100 times to be able to find the best subset of probes.
    pool=mp.Pool(nr_of_cores)
    chosen_probes=[]
    chosen_probes=[pool.apply(choose_probes_from_scores,args=(probe_scores[i],possible_arm_combinations,permutations,counter,probe_id_list[i],seed)) for i,possible_arm_combinations in enumerate(possible_probes_all_targets)] #function that choses probes by the probability which is based on the score
    pool.close()
    print('\tRound '+str(counter)+' :Probes chosen')
    
    chosen_probes_lists.append(chosen_probes)
    probes_with_dimers,conflicting_targets=find_dimer_forming_probes_iterative(chosen_probes)
    nr_of_dimer_probes=len(conflicting_targets)
    probes_with_dimers_lists.append(nr_of_dimer_probes)
    print('\t'+str(nr_of_dimer_probes)+' out of '+str(len(chosen_probes))+' probes form a dimer')
    min_dimers=min(nr_of_dimer_probes,min_dimers)
    
    #Find hairpins
    

    #exclude probes that form hairpins from being chosen
    for probe in probes_with_dimers:
        probe_index=probe_cpg_id_list.index(probe.split(':')[0])
        nested_index=probe_id_list[probe_index].index(probe)
        probe_scores[probe_index][nested_index]+=10^10
    #exclude probes that form dimers from being chosen (try around with e.g. top 25% or top 1)
    

    counter+=1
chosen_set_index=probes_with_dimers_lists.index(min(probes_with_dimers_lists)) #Also check the sum of all probe scores?
chosen_set=chosen_probes_lists[chosen_set_index]
seq_list=[]

#------------------------------------------------------------------------
#write output

print('minimal number of dimers is '+str(min_dimers))
with open(targets,'r') as handle:
    reader=csv.reader(handle,delimiter='\t')
    description_list=[]
    for row in reader:
        description_list.append(str(row[0]+'\t'+row[-1]))

with open(output_name,'w') as handle:
    for i,probe in enumerate(chosen_set):
        if probe == None:
            pass
        else:
            probename=str(i)
            probe_description=description_list[i]
            seq_list.append(SeqRecord.SeqRecord(Seq.Seq(probe[0]),id=probename,description=probe[1]))
    SeqIO.write(seq_list,handle,'fasta')




