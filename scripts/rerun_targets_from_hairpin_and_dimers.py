#!/usr/bin/python3
#Form a list of targets to for which a new probe must be designed due to hairpin or dimer formation.
from argparse import ArgumentParser
import json
import csv
import glob
import re
parser = ArgumentParser()
parser.add_argument("-p", "--hairpin_json_file", dest="hairpin_file",
                    help="open FILE", metavar="HAIRPIN_JSON_FILE")
parser.add_argument("-d", "--dimer_json_file", dest="dimer_file",
                    help="open FILE", metavar="HAIRPIN_JSON_FILE")
parser.add_argument("-o", "--output_hairpin", dest="outputname_hairpin",
                    help="write report to FILE", metavar="OUTPUTFILE")
parser.add_argument("-c", "--output_dimer", dest="outputname_dimer",
                    help="write report to FILE", metavar="OUTPUTFILE")
parser.add_argument("-i", "--config_file", dest="config_file",
                    help="Configuration file in JSON format", metavar="CONFIGFILE")
args = vars(parser.parse_args())

hairpin_file=args["hairpin_file"]
dimer_file=args["dimer_file"]
outputname_hairpin=args["outputname_hairpin"]
outputname_dimer=args["outputname_dimer"]
config_file=args["config_file"]
with open(config_file) as jsonFile:
    configObject = json.load(jsonFile)
    jsonFile.close()
Backbone_sequence=configObject['Backbone_sequence']
backbone_length=len(Backbone_sequence)
probe_id_list=[]
arm_upstream_loc_list=[]
arm_upstream_list=[]
arm_downstream_loc_list=[]
arm_downstream_list=[]
output_dir=outputname_dimer.split('/')[0]
path=output_dir+'/iteration_*/cpg_snp_*_arms.tsv'#Open all output/iteration_*/cpg_snp_*_arms.tsv
for filename in glob.glob(path):
    with open(filename) as handle:
        handle.readline() #read the header
        reader=csv.reader(handle,delimiter='\t')
        for row in reader:
            probe_id_list.append(row[4])            #Create a list of all CpGId_Arm_ids
            arm_upstream_loc_list.append(row[2])    #Create a list of all downstream arm ranges
            arm_downstream_loc_list.append(row[3])  #Create a list of all upstream arm ranges
            arm_upstream_list.append(row[0])
            arm_downstream_list.append(row[1])
with open(hairpin_file) as jsonFile:
    hairpinlist = json.load(jsonFile)
    jsonFile.close()
with open(dimer_file) as jsonFile:
    dimerlist = json.load(jsonFile)
    jsonFile.close()
conflicting_cpg_hairpin_list=[]
conflicting_cpg_hairpin_start_range=[]
conflicting_cpg_hairpin_end_range=[]
if hairpinlist is None:
    pass
else:
    for hairpin in hairpinlist:
        cpg_id=hairpin['Seq']['ID']
        hairpin_seq=hairpin['Seq']['Seq']
        hairpin_A_seq=hairpin['Aseq']
        chrom=arm_upstream_loc_list[probe_id_list.index(cpg_id)].split(':')[0]
        #find the genomic position of the sequence that forms the start of the hairpin sequence
        hairpin_conflict_start=[m.start() for m in re.finditer('/',hairpin_A_seq)][0]
        hairpin_conflict_end=[m.start() for m in re.finditer('/',hairpin_A_seq)][-1]
        arm_loc=arm_upstream_loc_list[probe_id_list.index(cpg_id)].split(':')[1]
        arm_range=[arm_loc.split('-')[0],arm_loc.split('-')[1]]
        hairpin_start_range=int(arm_range[0])+hairpin_conflict_start        #Convert conflict start to range
        hairpin_end_range=int(arm_range[0])+hairpin_conflict_end        #convert conflict end to range
        if hairpin_start_range < int(arm_range[0]): #if the hairpin extends into the backbone, return the end location of the arm, as another arm will be chosen.
            hairpin_start_range= int(arm_range[0])
        else:
            pass

        #perform the same for the rigth side of the hairpin
        total_length_probe=len(hairpin_seq)
        arm_loc2=arm_downstream_loc_list[probe_id_list.index(cpg_id)].split(':')[1]
        arm_range2=[arm_loc2.split('-')[0],arm_loc2.split('-')[1]]
        arm_length2=int(arm_range2[1])-int(arm_range2[0])+1
        arm_start2=total_length_probe-arm_length2
        hairpin_conflict_start2=[m.start() for m in re.finditer('\\\\',hairpin_A_seq)][0]-arm_start2
        hairpin_conflict_end2=[m.start() for m in re.finditer('\\\\',hairpin_A_seq)][-1]-arm_start2
        hairpin_start_range2=int(arm_range2[0])+hairpin_conflict_start2        #Convert conflict start to range
        hairpin_end_range2=int(arm_range2[0])+hairpin_conflict_end2        #convert conflict end to range

        if hairpin_start_range2 < int(arm_range2[0]): #if the hairpin start from the backbone, return the start location of the arm, as another arm will be chosen.
            hairpin_start_range2= int(arm_range2[0])
        else:
            pass
        conflicting_cpg_hairpin_list.append(cpg_id)
        if hairpin_start_range > int(arm_range[1]): #The start of the hairpin is not in the upstream arm.
            conflicting_cpg_hairpin_start_range.append('Outside arm range')
        else:
            conflicting_cpg_hairpin_start_range.append(chrom+':'+str(hairpin_start_range)+'-'+str(hairpin_end_range))
        if hairpin_end_range2 < int(arm_range2[0]): #The end of the hairpin is not in the downstream arm.
            conflicting_cpg_hairpin_end_range.append('Outside arm range')
        else:
            conflicting_cpg_hairpin_end_range.append(chrom+':'+str(hairpin_start_range2)+'-'+str(hairpin_end_range2))
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
        S1_upstream_arm_length=int(arm_upstream_loc_list[probe_id_list.index(cpg_id)].split(':')[1].split('-')[1])-int(arm_upstream_loc_list[probe_id_list.index(cpg_id)].split(':')[1].split('-')[0])+1
        S1_downstream_arm_length=int(arm_downstream_loc_list[probe_id_list.index(cpg_id)].split(':')[1].split('-')[1])-int(arm_downstream_loc_list[probe_id_list.index(cpg_id)].split(':')[1].split('-')[0])+1
        S2_upstream_arm_length=int(arm_upstream_loc_list[probe_id_list.index(cpg_id2)].split(':')[1].split('-')[1])-int(arm_upstream_loc_list[probe_id_list.index(cpg_id2)].split(':')[1].split('-')[0])+1
        S2_downstream_arm_length=int(arm_downstream_loc_list[probe_id_list.index(cpg_id2)].split(':')[1].split('-')[1])-int(arm_downstream_loc_list[probe_id_list.index(cpg_id2)].split(':')[1].split('-')[0])+1
        #check whether the S1 and S2 sequences are 5' to 3' or the other way around.
#        print('\n')
#        print(dimer['S1']['Seq'][:S1_upstream_arm_length])
#        print(arm_upstream_loc_list[probe_id_list.index(cpg_id)])
        if arm_upstream_list[probe_id_list.index(cpg_id)] in dimer['S1']['Seq'][:S1_upstream_arm_length]:
#            print('yes')
            dimer_start_num_corrected=dimer_start_num
            dimer_end_num_corrected=dimer_end_num

        else: #In this situation the S1 sequence is 3'-5', therefore we need to correct the indices.
            dimer_start_num_corrected=S1_downstream_arm_length-(S1_length-dimer_start_num)
            dimer_end_num_corrected=S1_downstream_arm_length-(S1_length-dimer_end_num)
        if arm_upstream_list[probe_id_list.index(cpg_id2)] in dimer['S2']['Seq'][:S2_upstream_arm_length]:
#            print('yes2')
            dimer_start_num2_corrected=dimer_start_num2
            dimer_end_num2_corrected=dimer_end_num2
        else: #In this situation the S2 sequence is 3'-5', therefore we need to correct the indices.
            dimer_start_num2_corrected=S2_length-dimer_end_num2
            dimer_end_num2_corrected=S2_length-dimer_start_num2
        
#The next two if statements only include the situation if dimer['S1Dangling']== 0: (or the other way around?)
#By this I mean that dimer_range should always give the arm at the start of the designed probe and dimer_range2 should always give the arm at the end of the probe. Although this is only the case if dimer['S1Dangling']== 0 AND both the S1 and S2 sequences are 5'-3'. Look into this.
        #obtain the genomic locations of the dimer forming regions for S1
        if arm_upstream_list[probe_id_list.index(cpg_id)] in dimer['S1']['Seq'][:S1_upstream_arm_length]: #Check whether the dimer was from 3'-5'
#            print('yes3')
            if dimer_start_num_corrected < backbone_length: #dimer is in upstream arm
                chrom=arm_upstream_loc_list[probe_id_list.index(cpg_id)].split(':')[0]
                start=int(arm_upstream_loc_list[probe_id_list.index(cpg_id)].split(':')[1].split('-')[0])
                end=int(arm_upstream_loc_list[probe_id_list.index(cpg_id)].split(':')[1].split('-')[1])
                dimer_start_range=max(start,start+dimer_start_num_corrected)
                dimer_end_range=min(end,start+dimer_end_num_corrected)
                dimer_range=str(chrom)+':'+str(dimer_start_range)+'-'+str(dimer_end_range)
            else:   #dimer is in downstream arm
                chrom=arm_downstream_loc_list[probe_id_list.index(cpg_id)].split(':')[0]
                start=int(arm_downstream_loc_list[probe_id_list.index(cpg_id)].split(':')[1].split('-')[0])
                end=int(arm_downstream_loc_list[probe_id_list.index(cpg_id)].split(':')[1].split('-')[1])
                dimer_start_range=max(start,start+dimer_start_num_corrected-(S1_length-(end-start+1)))
                dimer_end_range=min(end,start+dimer_end_num_corrected-(S1_length-(end-start+1)))
                dimer_range=str(chrom)+':'+str(dimer_start_range)+'-'+str(dimer_end_range)
        else:
            if dimer_start_num_corrected < backbone_length: #dimer is in upstream arm 
                chrom=arm_downstream_loc_list[probe_id_list.index(cpg_id)].split(':')[0]
                start=int(arm_downstream_loc_list[probe_id_list.index(cpg_id)].split(':')[1].split('-')[0])
                end=int(arm_downstream_loc_list[probe_id_list.index(cpg_id)].split(':')[1].split('-')[1])
                dimer_start_range=max(start,start+dimer_start_num_corrected)
                dimer_end_range=min(end,start+dimer_end_num_corrected)
                dimer_range=str(chrom)+':'+str(dimer_start_range)+'-'+str(dimer_end_range)
            else: #dimer is in downstream arm
                chrom=arm_downstream_loc_list[probe_id_list.index(cpg_id)].split(':')[0]
                start=int(arm_downstream_loc_list[probe_id_list.index(cpg_id)].split(':')[1].split('-')[0])
                end=int(arm_downstream_loc_list[probe_id_list.index(cpg_id)].split(':')[1].split('-')[1])
                dimer_start_range=max(start,start+dimer_start_num_corrected-(S1_length-(end-start+1)))
                dimer_end_range=min(end,start+dimer_end_num_corrected-(S1_length-(end-start+1)))
                dimer_range=str(chrom)+':'+str(dimer_start_range)+'-'+str(dimer_end_range)
    #obtain the genomic locations of the dimer forming regions for S2
        if arm_upstream_list[probe_id_list.index(cpg_id2)] in dimer['S2']['Seq'][:S2_upstream_arm_length]: #Check whether the dimer was from 3'-5'
#            print('yes4')
            if dimer_start_num2_corrected < backbone_length: #dimer is in upstream arm
                chrom=arm_upstream_loc_list[probe_id_list.index(cpg_id2)].split(':')[0]
                start=int(arm_upstream_loc_list[probe_id_list.index(cpg_id2)].split(':')[1].split('-')[0])
                end=int(arm_upstream_loc_list[probe_id_list.index(cpg_id2)].split(':')[1].split('-')[1])
                dimer_start_range=max(start,start+dimer_start_num2_corrected)
                dimer_end_range=min(end,start+dimer_end_num2_corrected)
                dimer_range2=str(chrom)+':'+str(dimer_start_range)+'-'+str(dimer_end_range)
            else: #dimer is in downstream arm
                chrom=arm_downstream_loc_list[probe_id_list.index(cpg_id2)].split(':')[0]
                start=int(arm_downstream_loc_list[probe_id_list.index(cpg_id2)].split(':')[1].split('-')[0])
                end=int(arm_downstream_loc_list[probe_id_list.index(cpg_id2)].split(':')[1].split('-')[1])
                dimer_start_range=max(start,start+dimer_start_num2_corrected-(S2_length-(end-start+1)))
                dimer_end_range=min(end,start+dimer_end_num2_corrected-(S2_length-(end-start+1)))
                dimer_range2=str(chrom)+':'+str(dimer_start_range)+'-'+str(dimer_end_range)
        else: #If the sequence was obtained in 3'-5' from the dimer file, use downstream_loc_list
            if dimer_start_num2_corrected < backbone_length: #dimer is in upstream arm
                chrom=arm_upstream_loc_list[probe_id_list.index(cpg_id2)].split(':')[0]
                start=int(arm_upstream_loc_list[probe_id_list.index(cpg_id2)].split(':')[1].split('-')[0])
                end=int(arm_upstream_loc_list[probe_id_list.index(cpg_id2)].split(':')[1].split('-')[1])
                dimer_start_range=max(start,start+dimer_start_num2_corrected)
                dimer_end_range=min(end,start+dimer_end_num2_corrected)
                dimer_range2=str(chrom)+':'+str(dimer_start_range)+'-'+str(dimer_end_range)
            else: #dimer is in downstream arm
                chrom=arm_downstream_loc_list[probe_id_list.index(cpg_id2)].split(':')[0]
                start=int(arm_downstream_loc_list[probe_id_list.index(cpg_id2)].split(':')[1].split('-')[0])
                end=int(arm_downstream_loc_list[probe_id_list.index(cpg_id2)].split(':')[1].split('-')[1])
                dimer_start_range=max(start,start+dimer_start_num2_corrected-(S2_length-(end-start+1)))
                dimer_end_range=min(end,start+dimer_end_num2_corrected-(S2_length-(end-start+1)))
                dimer_range2=str(chrom)+':'+str(dimer_start_range)+'-'+str(dimer_end_range)
        conflict_range_dimer.append(dimer_range)
        conflict_range_dimer2.append(dimer_range2)
#        if cpg_id == 'CpG24941838_arm_8496' or cpg_id2=='CpG24941838_arm_8496':
#            print(dimer)
#            print(cpg_id)
#            print('dimer start num is '+str(dimer_start_num_corrected))
#            print('dimer end num is '+str(dimer_end_num_corrected))
#            print(dimer_range)
#            print(cpg_id2)
#            print('dimer start num is '+str(dimer_start_num2_corrected))
#            print('dimer end num is '+str(dimer_end_num2_corrected))
#            print(dimer_range2)
#        print('\n')
#Return a file with all of the CpG_id's that form hairpins or dimers
#only return the first mentioned probe in a dimer
with open(outputname_hairpin, 'w') as handle:
    writer=csv.writer(handle,delimiter='\t')
    header=['CpG_ID','hairpin_range1','hairpin_range2']
    writer.writerow(header)
    for i,cpg_id in enumerate(conflicting_cpg_hairpin_list):
        writer.writerow([cpg_id,conflicting_cpg_hairpin_start_range[i],conflicting_cpg_hairpin_end_range[i]])

with open(outputname_dimer, 'w') as handle:
    writer=csv.writer(handle,delimiter='\t')
    header=['CpG_ID1','dimer_range1','CpG_ID2','dimer_range2']
    writer.writerow(header)
    for i,cpg_id in enumerate(conflicting_cpg_list_dimer):
        writer.writerow([cpg_id,conflict_range_dimer[i],conflicting_cpg_list_dimer2[i],conflict_range_dimer2[i]])
