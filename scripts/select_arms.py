#!/usr/bin/python3
#Compute all possible forward and reverse primers.
import os
from argparse import ArgumentParser
from operator import itemgetter
import csv
import numpy as np
parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="filename",
                    help="open FILE", metavar="FILE")
parser.add_argument("-o", "--output", dest="outputname",
                    help="write report to FILE", metavar="OUTPUTFILE")
parser.add_argument("-t", "--target_file", dest="target_file",
                    help="target_file", metavar="TARGET_FILE")
parser.add_argument("-n", "--output_not_selected", dest="outputname_not_selected",
                    help="write report to FILE", metavar="OUTPUTFILE_NOT_SELECTED")
parser.add_argument("-y", "--iteration", dest="iteration",
                    help="input value of iteration", metavar="ITERATION")
parser.add_argument("-c", "--conflicts_hairpin", dest="conflicts_hairpin",
                    help="file containing cpgid+arm_id of probes forming hairpins or dimers", metavar="CONFLICTS")
parser.add_argument("-d", "--conflicts_dimer", dest="conflicts_dimer",
                    help="file containing cpgid+arm_id of probes forming hairpins or dimers", metavar="CONFLICTS")
parser.add_argument("-s", "--selection_round", dest="selection_round",
                    help="input value of selection_round", metavar="SELECTION_ROUND")
parser.add_argument("-a", "--input_up", dest="input_up",
                    help="input_up", metavar="input_up")
parser.add_argument("-b", "--input_down", dest="input_down",
                    help="input_down", metavar="input_down")
args = vars(parser.parse_args())
filename=args["filename"]
outputname=args["outputname"]
target_file=args["target_file"]
outputname_not_selected=args["outputname_not_selected"]
iteration=int(args["iteration"])
selection_round=int(args["selection_round"])
conflict_file_hairpin=args["conflicts_hairpin"]
conflict_file_dimer=args["conflicts_dimer"]
inputname_up=args["input_up"]
inputname_down=args["input_down"]
tm_up=[]
tm_down=[]
with open(inputname_up) as handle:
    reader=csv.reader(handle)
    for row in reader:
        tm_up.append(float(row[0]))
with open(inputname_down) as handle:
    reader=csv.reader(handle)
    for row in reader:
        tm_down.append(float(row[0]))
preferred_tm=np.mean([np.max(tm_down),np.max(tm_up),np.min(tm_down),np.min(tm_up)])

conflicting_cpg_arm_id_list=[]
conflicting_hairpins_start_range=[]
conflicting_hairpins_end_range=[]
conflicting_cpg_id_list=[]
#Construct list of conflicting arms due to hairpins and dimers
with open(conflict_file_hairpin) as handle:
    reader=csv.reader(handle,delimiter='\t')
    handle.readline()
    for row in reader:
        conflicting_cpg_arm_id_list.append(row[0])
        conflicting_cpg_id_list.append(row[0].split('_')[0])
        conflicting_hairpins_start_range.append(row[1])
        conflicting_hairpins_end_range.append(row[2])
with open(conflict_file_dimer) as handle:
    reader=csv.reader(handle,delimiter='\t')
    handle.readline()
    for row in reader:
        conflicting_cpg_arm_id_list.append(row[0])
        conflicting_cpg_id_list.append(row[0].split('_')[0])
        conflicting_hairpins_start_range.append(row[1])
        conflicting_hairpins_end_range.append(row[1])

#        conflicting_cpg_arm_id_list.append(row[2])
#        conflicting_cpg_id_list.append(row[2].split('_')[0])
#        conflicting_hairpins_start_range.append(row[3])
#        conflicting_hairpins_end_range.append('Outside arm range')
arm_list=[]
tm_list=[]
cpg_list=[]
last_cpg_id=0
selected_arms=[]
with open(filename) as handle:
    header=handle.readline().rstrip('\n').split('\t')
    reader=csv.reader(handle,delimiter="\t")
    for i,row in enumerate(reader):
        passed_boolean=False
        cpg_arm_id=row[4]
        cpg_id=cpg_arm_id.split('_')[0]
        if cpg_arm_id in conflicting_cpg_arm_id_list:
            passed_boolean=True
        elif cpg_id in conflicting_cpg_id_list:
            part_of_conflicting_cpg_arm_list=[x for index, x in enumerate(conflicting_cpg_arm_id_list) if x.split('_')[0]==cpg_id]
            part_of_conflicting_cpg_arm_list_i=[index for index, x in enumerate(conflicting_cpg_arm_id_list) if x.split('_')[0]==cpg_id]
            for index, conflicting_cpg_arm_id in enumerate(part_of_conflicting_cpg_arm_list):#loop over all possible arms that are conflicting (CpG_id_arms)
                #Add restrictions on which ranges from both arms cannot be combined for hairpins 
                list_index=part_of_conflicting_cpg_arm_list_i[index]
                if conflicting_hairpins_start_range[list_index] == 'Outside arm range':
                    hairpin_region_downstream=0
                    conflict_start_range='Outside arm range'
                    arm_downstream_loc=row[2]
                    downstream_range=range(int(arm_downstream_loc.split(':')[1].split('-')[0]),int(arm_downstream_loc.split(':')[1].split('-')[1]))
                else:
                    arm_downstream_loc=row[2]
                    downstream_range=range(int(arm_downstream_loc.split(':')[1].split('-')[0]),int(arm_downstream_loc.split(':')[1].split('-')[1]))
                    conflict_start_range1=conflicting_hairpins_start_range[list_index]
                    conflict_start_range=range(int(conflict_start_range1.split(':')[1].split('-')[0]),int(conflict_start_range1.split(':')[1].split('-')[1]))
                    hairpin_region_downstream=len(set(downstream_range).intersection(conflict_start_range))
                if hairpin_region_downstream>6: #if hairpin region in downstream arm is longer than 1 nucleotide: pass
                    passed_boolean=True
                    break
                else:
                    if conflicting_hairpins_end_range[list_index] == 'Outside arm range':
                        hairpin_region_upstream=0
                        conflict_end_range='Outside arm range'
                        arm_upstream_loc=row[3]
                        upstream_range=range(int(arm_upstream_loc.split(':')[1].split('-')[0]),int(arm_upstream_loc.split(':')[1].split('-')[1]))
                    else:
                        arm_upstream_loc=row[3]
                        upstream_range=range(int(arm_upstream_loc.split(':')[1].split('-')[0]),int(arm_upstream_loc.split(':')[1].split('-')[1]))
                        conflict_end_range1=conflicting_hairpins_end_range[list_index]
                        conflict_end_range=range(int(conflict_end_range1.split(':')[1].split('-')[0]),int(conflict_end_range1.split(':')[1].split('-')[1]))
                        hairpin_region_upstream=len(set(upstream_range).intersection(conflict_end_range))
                    if hairpin_region_upstream>6: ##if hairpin region in upstream arm is longer than 1 nucleotide: pass
                        passed_boolean=True
                        break
                    else: #In this situation a new probe will be chosen which has no significant hairpins.
                        pass
        else:
            pass
        if passed_boolean==False:
            tm_down=float(row[6])
            tm_up=float(row[7])
            delta_tm_score=abs(preferred_tm-tm_up)**2+abs(preferred_tm-tm_down)**2+abs(tm_up-tm_down)**2
            if cpg_id==last_cpg_id or i==0 or len(tm_list)==0:
                arm_list.append(row)
                tm_list.append(delta_tm_score)
                cpg_list.append(cpg_id)
            else:
                #find best arm from previous loop
                indexx=min(enumerate(tm_list),key=itemgetter(1))[0]
                selected_arms.append(arm_list[indexx]) #select the arms with the lowest Tm
                #empty all lists for next loop
                arm_list=[row]
                tm_list=[delta_tm_score]
                cpg_list.append(cpg_id)
        else:
            pass
        last_cpg_id=cpg_id
 
#Also select the arms for the last CpG
if len(cpg_list) == 0:
    pass
elif passed_boolean==False: 
    index=min(enumerate(tm_list),key=itemgetter(1))[0]
    selected_arms.append(arm_list[index]) #select the arms with the lowest Tm
else:
    pass

#Write a file returning the selected arms
with open(outputname, "w") as handle:
    writer=csv.writer(handle, delimiter="\t")
    header.append('iteration')
    writer.writerow(header)
    for row in selected_arms:
        row.append(str(iteration))
        writer.writerow(row)


for i in range(0,iteration):
    previous_selected_cg_ids=[]
    previous_selected_filename="output/iteration_"+str(i)+"/selected_arms_"+str(i)+"_"+str(selection_round)+".tsv"
    with open(previous_selected_filename) as handle:
        handle.readline()
        reader=csv.reader(handle,delimiter='\t')
        for row in reader:
            cpg_id=row[4].split('_')[0]
            cpg_list.append(cpg_id)

cpg_target_list=[]
#Return a file with all of the CpG_id's for which no arms were selected.
with open(target_file) as handle:
    reader=csv.reader(handle,delimiter='\t')
    for row in reader:
        cpg_id=row[-1]
        cpg_target_list.append(cpg_id)
with open(outputname_not_selected, 'w') as handle:
    writer=csv.writer(handle)
    not_selected=set(cpg_target_list).symmetric_difference(cpg_list)
    for row in not_selected:
        writer.writerow([row])
