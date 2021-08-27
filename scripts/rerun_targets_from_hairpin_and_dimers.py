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

args = vars(parser.parse_args())

hairpin_file=args["hairpin_file"]
dimer_file=args["dimer_file"]
outputname_hairpin=args["outputname_hairpin"]
outputname_dimer=args["outputname_dimer"]

probe_id_list=[]
arm_upstream_loc_list=[]
arm_downstream_loc_list=[]
path='output/iteration_*/cpg_snp_*_arms.tsv'#Open all output/iteration_*/cpg_snp_*_arms.tsv
for filename in glob.glob(path):
    with open(filename) as handle:
        handle.readline() #read the header
        reader=csv.reader(handle,delimiter='\t')
        for row in reader:
            probe_id_list.append(row[4])            #Create a list of all CpGId_Arm_ids
            arm_upstream_loc_list.append(row[2])    #Create a list of all downstream arm ranges
            arm_downstream_loc_list.append(row[3])  #Create a list of all upstream arm ranges
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

if dimerlist is None:
    pass
else:
    for dimer in dimerlist:  
        cpg_id=dimer['S1']['ID']
        cpg_id2=dimer['S2']['ID']
        conflicting_cpg_list_dimer.append(cpg_id)
        conflicting_cpg_list_dimer2.append(cpg_id2)
conflicting_cpg_list_dimer=list(dict.fromkeys(conflicting_cpg_list_dimer))
#conflicting_cpg_list2=list(dict.fromkeys(conflicting_cpg_list2))

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
    for cpg_id in conflicting_cpg_list_dimer:
        writer.writerow([cpg_id])

