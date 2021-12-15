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
parser.add_argument("-u", "--output_dir", dest="output_dir",
                    help="Directory for output", metavar="OUTPUTDIR")
parser.add_argument("-f", "--conflicting_output_file", dest="conflicting_file",
                    help="Output filename for conflicting probes per iteration", metavar="CONFLICTS")
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
outputdir=args["output_dir"]+'/'
conflicting_file=args["conflicting_file"]
with open(config_file) as jsonFile:
    configObject = json.load(jsonFile)
    jsonFile.close()
permutations=int(configObject['Permutations'])
backbone_length=len(configObject["Backbone_sequence"][0]["Reverse_complement_universal_forward_primer"])+len(configObject["Backbone_sequence"][0]["Common_region"])+len(configObject["Backbone_sequence"][0]["Universal_reverse_primer"])
exclusion_factor=float(configObject["Dimer_Exclusion_Factor"])

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
with open(probe_arms_file,'rb') as handle:
    possible_arm_combinations_all_targets=pickle.load(handle)
    probe_arm_list=[]
    arm_upstream_loc_list=[]
    arm_downstream_loc_list=[]
    probe_id_list=[]
    arm_upstream_list=[]
    arm_downstream_list=[]
    for i, row in enumerate(possible_arm_combinations_all_targets):
        probe_arm_list.append([])
        arm_upstream_loc_list.append([])
        arm_downstream_loc_list.append([])
        probe_id_list.append([])
        arm_upstream_list.append([])
        arm_downstream_list.append([])
        for probe in row:
            probe_arm_list[i].append(probe)
            arm_upstream_loc_list[i].append(probe[2])
            arm_downstream_loc_list[i].append(probe[3])
            probe_id_list[i].append(str(probe[7])+':'+str(probe[4]))
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

def get_conflict_ranges_S1(conflicting_cpg_list_dimer,conflict_range_dimer):
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
    amount_to_exlude=int(int(len(sorted_probe_list)*exclusion_factor))
    if amount_to_exlude==0:
        most_conflicting_probe_list=sorted_probe_list
    else:
        most_conflicting_probe_list=sorted_probe_list[0:int(len(sorted_probe_list)*exclusion_factor)] #exclude the top portion (top 10% if exclusion factor== 0.1
    most_conflicting_dimer_range_list=[]
    for probe in most_conflicting_probe_list: #report the conflicting locus at first occurrence of dimer
        index=conflicting_cpg_list_dimer.index(probe)
        most_conflicting_dimer_range_list.append(conflict_range_dimer[index])
    return most_conflicting_dimer_range_list

def get_conflict_ranges_S1_and_S2(conflicting_cpg_list_dimer,conflict_range_dimer,conflicting_cpg_list_dimer2,conflict_range_dimer2):
    most_conflicting_dimer_range_list=get_conflict_ranges_S1(conflicting_cpg_list_dimer,conflict_range_dimer)
    most_conflicting_dimer_range_list2=get_conflict_ranges_S1(conflicting_cpg_list_dimer2,conflict_range_dimer2)
    combined_most_conflicting_dimer_range_list=most_conflicting_dimer_range_list+most_conflicting_dimer_range_list2
    return combined_most_conflicting_dimer_range_list

def get_dimer_range_upstream(arm_loc_list,probe_index,nested_index,dimer_start_num,dimer_end_num,dimer):
    chrom=arm_upstream_loc_list[probe_index][nested_index].split(':')[0]
    start=int(arm_upstream_loc_list[probe_index][nested_index].split(':')[1].split('-')[0].strip("'"))
    end=int(arm_upstream_loc_list[probe_index][nested_index].split(':')[1].split('-')[1].strip("'"))
    template_strand=possible_arm_combinations_all_targets[probe_index][nested_index][-2]
    dimer_start_range=max(start,start+dimer_start_num)
    if template_strand=='+':
        dimer_end_range=min(end,start+dimer_end_num)
    else:
        dimer_end_range=max(end,start+dimer_end_num)
    dimer_range=str(chrom)+':'+str(dimer_start_range)+'-'+str(dimer_end_range)
    if dimer_start_range > dimer_end_range:
        print('error '+str([dimer_start_range,dimer_end_range,dimer]))
    return [dimer_range,template_strand]

def get_dimer_range_downstream(arm_loc_list,probe_index,nested_index,dimer_start_num,dimer_end_num,S_length,dimer):
    chrom=arm_downstream_loc_list[probe_index][nested_index].split(':')[0]
    start=int(arm_downstream_loc_list[probe_index][nested_index].split(':')[1].split('-')[0].strip("'"))
    end=int(arm_downstream_loc_list[probe_index][nested_index].split(':')[1].split('-')[1].strip("'"))
    template_strand=possible_arm_combinations_all_targets[probe_index][nested_index][-2]
    dimer_start_range=max(start,start+dimer_start_num-(S_length-(end-start+1)))
    if template_strand=='+':
        dimer_end_range=min(end,start+dimer_end_num-(S_length-(end-start+1)))
    else:
        dimer_end_range=max(end,start+dimer_end_num-(S_length-(end-start+1)))
    dimer_range=str(chrom)+':'+str(dimer_start_range)+'-'+str(dimer_end_range)
    if dimer_start_range > dimer_end_range:
        print('error '+str([dimer_start_range,dimer_end_range,dimer]))
    return [dimer_range,template_strand]

def write_nested_loci_to_bedfile(output_name,loci_list):
    to_write_list=[[locus.split(':')[0],locus.split(':')[1].split('-')[0],locus.split('-')[1],'','',probe_arm_list[i][j][-2],str(i)+':'+str(j)] for i,row in enumerate(loci_list) for j,locus in enumerate(row)] #[chrom,start,end,i:j] e.g. [1,105386,105396,0:304] for chr1:105386-105396 in the locus at loci_list[0][304]
    #to_write_list=[[locus.split(':')[0][4:],locus.split(':')[1].split('-')[0],locus.split('-')[1]] for i,row in enumerate(loci_list) for j,locus in enumerate(row)]
    with open(output_name,'w') as handle:
        writer=csv.writer(handle,delimiter='\t')
        writer.writerows(to_write_list)
    return output_name

def write_loci_to_bedfile(output_name,loci_list):
    to_write_list=[[locus[0].split(':')[0],locus[0].split(':')[1].split('-')[0],locus[0].split('-')[1],'','',locus[1],str(i)] for i,locus in enumerate(loci_list)] #[chrom,start,end,i] e.g. [1,105386,105396,0] for chr1:105386-105396 in the locus at loci_list[0]
    #to_write_list=[[locus.split(':')[0][4:],locus.split(':')[1].split('-')[0],locus.split('-')[1]] for i,locus in enumerate(loci_list)]
    with open(output_name,'w') as handle:
        writer=csv.writer(handle,delimiter='\t')
        writer.writerows(to_write_list)
    return output_name

def create_conflicting_indices_list_bedtools(loci_list_up,loci_list_down,loci_up_bedfile_name,loci_down_bedfile_name,dimer_range_list,dimer_bedfile_name,bedfile_intersect_name):
    write_nested_loci_to_bedfile(loci_up_bedfile_name,loci_list_up)
    write_nested_loci_to_bedfile(loci_down_bedfile_name,loci_list_down)
    write_loci_to_bedfile(dimer_bedfile_name,dimer_range_list)
    #combine loci_down_bedfile_name and loci_up_bedfile_name
    os.system('cat '+loci_up_bedfile_name+' '+loci_down_bedfile_name+' > '+outputdir+'combined.bed')
    #Do bedtools intersect here on combined_loci_bedfile_name and dimer_bedfile_name
    os.system('bedtools intersect -wa -s -a '+outputdir+'combined.bed -b '+dimer_bedfile_name+' -f 7E-9 > '+bedfile_intersect_name)
    
    #obtain identifier from bedfile_intersect_name
    new_conflicting_indices_list=[]
    with open(bedfile_intersect_name,'r') as handle:
        reader=csv.reader(handle,delimiter='\t')
        for row in reader:
            identifier=row[-1]
            i=identifier.split(':')[0]
            j=identifier.split(':')[1]
            new_conflicting_indices_list.append([i,j])
    return new_conflicting_indices_list


def rescore_dimer_forming_probes_iterative(chosen_probes):#Create a fasta file with the chosen probes
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
    os.system('mfeprimer dimer -i tmp_fasta_chosen_probes.fasta -j -o tmp_dimers_chosen_probes -s 7 -t 10') #score cut-off is 7, temperature minimum is 10 degrees
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
        for dindex,dimer in enumerate(dimerlist):
                cpg_id=dimer['S1']['ID']
                cpg_id2=dimer['S2']['ID']
                dimer_structure=dimer['Aseq']
                S1_length=len(dimer['S1']['Seq'])
                S2_length=len(dimer['S2']['Seq'])
                #Obtain the location of the dimer on the two sequences
                S1Dangling=dimer['S1Dangling']
                AseqDangling=dimer['AseqDangling']
                S2Dangling=dimer['S2Dangling']
                dimer_start_num=AseqDangling-S1Dangling
                dimer_end_num=dimer_start_num+len(dimer_structure)
                dimer_start_num2=AseqDangling-S2Dangling
                dimer_end_num2=dimer_start_num2+len(dimer_structure)
                conflicting_cpg_list_dimer.append(cpg_id)
                conflicting_cpg_list_dimer2.append(cpg_id2)
                probe_index=probe_cpg_id_list.index(cpg_id.split(':')[0])
                nested_index=probe_id_list[probe_index].index(cpg_id)
                probe_index2=probe_cpg_id_list.index(cpg_id2.split(':')[0])
                nested_index2=probe_id_list[probe_index2].index(cpg_id2)
                #obtain the genomic locations of the dimer forming regions for S1
                #check whether the S1 and S2 sequences are 5' to 3' or the other way around.)
                if arm_upstream_list[probe_index][nested_index].strip("'") in dimer['S1']['Seq'][:len(arm_upstream_list[probe_index][nested_index].strip("'"))]:
                    if dimer_start_num < len(arm_upstream_list[probe_index][nested_index].strip("'")): #dimer is in upstream arm
                        dimer_range=get_dimer_range_upstream(arm_upstream_loc_list,probe_index,nested_index,dimer_start_num,dimer_end_num,dimer)
                    elif dimer_end_num >= S1_length-len(arm_downstream_list[probe_index][nested_index].strip("'")):   #dimer is in downstream arm. If dimer end has no overlap with the downsteam or upstream arm is should be in the backbone.
                        dimer_range=get_dimer_range_downstream(arm_downstream_loc_list,probe_index,nested_index,dimer_start_num,dimer_end_num,S1_length,dimer)
                    else: #dimer forming region is in backbone, only report the other probe that forms the dimer together with this probe.
                        pass
                else: #In this situation the S1 sequence is 3'-5', therefore we need to correct the indices.
                    dimer_start_num,dimer_end_num = S1_length-dimer_end_num , S1_length-dimer_start_num
                    if dimer_start_num < len(arm_upstream_list[probe_index][nested_index].strip("'")): #dimer is in upstream arm
                        dimer_range=get_dimer_range_upstream(arm_upstream_loc_list,probe_index,nested_index,dimer_start_num,dimer_end_num,dimer)
                    elif dimer_end_num >= S1_length-len(arm_downstream_list[probe_index][nested_index].strip("'")): #dimer is in downstream arm If dimer end has no overlap with the downsteam or upstream arm is should be in the backbone.
                        dimer_range=get_dimer_range_downstream(arm_downstream_loc_list,probe_index,nested_index,dimer_start_num,dimer_end_num,S1_length,dimer)
                    else:#dimer forming region is in backbone, only report the other probe that forms the dimer together with this probe.
                        pass
                if arm_upstream_list[probe_index2][nested_index2].strip("'") in dimer['S2']['Seq'][:len(arm_upstream_list[probe_index2][nested_index2].strip("'"))]:
                    if dimer_start_num2 < len(arm_upstream_list[probe_index2][nested_index2].strip("'")): #dimer is in upstream arm
                        dimer_range2=get_dimer_range_upstream(arm_upstream_loc_list,probe_index2,nested_index2,dimer_start_num2,dimer_end_num2,dimer)
                    elif dimer_end_num2 >= S2_length-len(arm_downstream_list[probe_index2][nested_index2].strip("'")):   #dimer is in downstream arm. If dimer end has no overlap with the downsteam or upstream arm is should be in the backbone.
                        dimer_range2=get_dimer_range_downstream(arm_downstream_loc_list,probe_index2,nested_index2,dimer_start_num2,dimer_end_num2,S2_length,dimer)
                    else:#dimer forming region is in backbone, only report the other probe that forms the dimer together with this probe.
                        pass
                else: #In this situation the S2 sequence is 3'-5', therefore we need to correct the indices.
                    dimer_start_num2,dimer_end_num2 = S2_length-dimer_end_num2 , S2_length-dimer_start_num2
                    if dimer_start_num2 < len(arm_upstream_list[probe_index2][nested_index2].strip("'")): #dimer is in upstream arm
                        dimer_range2=get_dimer_range_upstream(arm_upstream_loc_list,probe_index2,nested_index2,dimer_start_num2,dimer_end_num2,dimer)
                    elif dimer_end_num2 >= S2_length-len(arm_downstream_list[probe_index2][nested_index2].strip("'")):   #dimer is in downstream arm. If dimer end has no overlap with the downsteam or upstream arm is should be in the backbone.
                        dimer_range2=get_dimer_range_downstream(arm_downstream_loc_list,probe_index2,nested_index2,dimer_start_num2,dimer_end_num2,S2_length,dimer)
                    else:#dimer forming region is in backbone, only report the other probe that forms the dimer together with this probe.
                        pass
                conflict_range_dimer.append(dimer_range)
                conflict_range_dimer2.append(dimer_range2)

    #Only report the probe that formed the most amount of dimers
    most_conflicting_dimer_range_list=get_conflict_ranges_S1_and_S2(conflicting_cpg_list_dimer,conflict_range_dimer,conflicting_cpg_list_dimer2,conflict_range_dimer2) #obtain conflict ranges from dimer forming probes
    
    #find all probes that have common regions with the most_conflicting_dimer_range and report those as conflicting probes
    all_conflicting_probe_list=[]
    conflicting_targets=set([])
    print('Dimer ranges obtained')
    loci_up_bedfile_name=outputdir+'loci_up.bed'
    loci_down_bedfile_name=outputdir+'loci_down.bed'
    dimer_bedfile_name=outputdir+'conflicting_dimer.bed'
    bedfile_intersect_name=outputdir+'conflicting_loci.bed'
    new_conflicting_indices_list=create_conflicting_indices_list_bedtools(arm_upstream_loc_list,arm_downstream_loc_list,loci_up_bedfile_name,loci_down_bedfile_name,most_conflicting_dimer_range_list,dimer_bedfile_name,bedfile_intersect_name)
    print('Conflicts found')
    all_conflicting_indices=[indices for indices_list in new_conflicting_indices_list for indices in indices_list]
    for indexes in new_conflicting_indices_list:
        i=int(indexes[0])
        j=int(indexes[1])
        all_conflicting_probe_list.append(probe_id_list[i][j])
        conflicting_targets.add(probe_id_list[i][0].split(':')[0])
        probe_scores[i][j]+=10^10 #Adjust score of conflicting probe
    return all_conflicting_probe_list,conflicting_targets


#-----------------------------------------------------------------------
#Iteratively exclude the most dimer forming probes.
chosen_probes_lists=[]
probes_with_dimers_lists=[]
conflicting_probe_list=[]
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
    probes_with_dimers,conflicting_targets=rescore_dimer_forming_probes_iterative(chosen_probes)
    nr_of_dimer_probes=len(conflicting_targets)
    probes_with_dimers_lists.append(nr_of_dimer_probes)
    conflicting_probe_list.append(probes_with_dimers)
    print('\t'+str(nr_of_dimer_probes)+' out of '+str(len(chosen_probes))+' probes form a dimer')
    min_dimers=min(nr_of_dimer_probes,min_dimers)
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
with open(conflicting_file,'w') as handle:
    writer=csv.writer(handle,delimiter='\t')
    for row in conflicting_probe_list:
        writer.writerow(row)
