#!/usr/bin/python3
#Compute all possible forward and reverse primers.
import os
from argparse import ArgumentParser
import json
from Bio import SeqIO, Seq, SeqRecord
from Bio.SeqUtils import GC
import multiprocessing as mp
import csv

#parse all files
parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="filename",
                    help="open FILE", metavar="FILE")
parser.add_argument("-o", "--output", dest="output_dir",
                    help="Directory for output", metavar="OUTPUTDIR")
parser.add_argument("-c", "--config_file", dest="config_file",
                    help="Config file with probe specifics", metavar="JSONfile")
parser.add_argument("-b", "--bed", dest="bedfile",
                    help="Bed file containing the targets", metavar="BED")
parser.add_argument("-t", "--cores", dest="cores",
                    help="Passed amount of cores to script", metavar="Cores")
args = vars(parser.parse_args())

filename=args["filename"]
bedfile=args["bedfile"]
outputdir=args["output_dir"]+'/'
config_file=args["config_file"]
nr_of_cores=int(args["cores"])

#import config file
with open(config_file) as jsonFile:
    configObject = json.load(jsonFile)
    jsonFile.close()
ftp_path_snp_database=configObject['ftp_path_snp_database']
probe_specifics=configObject['probe_specifics'][0]
min_arm_length=probe_specifics['min_arm_length']
max_arm_length=probe_specifics['max_arm_length']
min_target_length=probe_specifics['min_target_length']
max_target_length=probe_specifics['max_target_length']
target_range=configObject['target_range']
cpg_flanks=probe_specifics['cpg_flanks']
max_CG_percentage=float(probe_specifics['max_CG_percentage'])
min_CG_percentage=float(probe_specifics['min_CG_percentage'])

#patch backbone together
Backbone=configObject['Backbone_sequence'][0]
rev_compl_univeral_forward_primer=Backbone['Reverse_complement_universal_forward_primer']
common_region=Backbone['Common_region']
univeral_reverse_primer=Backbone['Universal_reverse_primer']
backbone_sequence=rev_compl_univeral_forward_primer+common_region+univeral_reverse_primer

mid_loc=int(target_range)+1 #The middle nucleotide is the target range +1 because the total target is created by adding the target range on both sides.
freq_threshold=float(configObject['SNP_frequency_threshold'])

#--------------------------------------------------------------
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
    delta_Tm=round(abs(Tm_up-Tm_down),1)
    return [Tm_up,Tm_down,delta_Tm]

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

def check_probe_for_hairpin_score(probes,fasta_name,outputname_json,index):
    seq_list=[]
    list_len=len(probes)
    if list_len == 0: # in case there is no probes
        return []
    with open(fasta_name+str(index)+'.fa','w') as output_handle:
        for i,probe in enumerate(probes):
            probename=str(i)
            seq_list.append(SeqRecord.SeqRecord(Seq.Seq(probe),id=probename))
        SeqIO.write(seq_list,output_handle,'fasta')
    os.system('mfeprimer hairpin --in '+fasta_name+str(index)+'.fa -j -o '+outputname_json+str(index)) #Remove/move the output files
    with open(outputname_json+str(index)+'.json') as jsonFile:
        hairpinlist = json.load(jsonFile)
        jsonFile.close()
    os.system('rm '+outputname_json+str(index)+'.json')
    os.system('rm '+outputname_json+str(index))
    os.system('rm '+fasta_name+str(index)+'.fa')
    hairpin_scores=[0 for i in range(list_len)]
    if hairpinlist is None:
        pass
    else:
        for hairpin in hairpinlist:
            hairpin_score=hairpin['Score']
            probe_nr=int(hairpin['Seq']['ID'])
            hairpin_scores[probe_nr]=hairpin_score
    return hairpin_scores

def loci_to_bed_37(locus):#https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt
    chrom=locus.split(':')[0][3:]
    left_index=locus.split(':')[1].split('-')[0]
    right_index=locus.split('-')[1]
    chromnrlist=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
    chromlist=['NC_000001.10','NC_000002.11','NC_000003.11','NC_000004.11','NC_000005.9','NC_000006.11','NC_000007.13','NC_000008.10','NC_000009.11','NC_000010.10','NC_000011.9','NC_000012.11','NC_000013.10','NC_000014.8','NC_000015.9','NC_000016.9','NC_000017.10','NC_000018.9','NC_000019.9','NC_000020.10','NC_000021.8','NC_000022.10','NC_000023.10','NC_000024.9']
    bed_info=[chromlist[chromnrlist.index(chrom)],left_index,right_index]
    return bed_info

def target_bed_to_bed_37(bedpath):
    chromnrlist=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
    chromlist=['NC_000001.10','NC_000002.11','NC_000003.11','NC_000004.11','NC_000005.9','NC_000006.11','NC_000007.13','NC_000008.10','NC_000009.11','NC_000010.10','NC_000011.9','NC_000012.11','NC_000013.10','NC_000014.8','NC_000015.9','NC_000016.9','NC_000017.10','NC_000018.9','NC_000019.9','NC_000020.10','NC_000021.8','NC_000022.10','NC_000023.10','NC_000024.9']
    bed_info_list=[]
    with open(bedpath) as handle:
        reader=csv.reader(handle,delimiter='\t')
        for row in reader:
            chrom=row[0][3:]
            bed_info=[chromlist[chromnrlist.index(chrom)],row[1],row[2],row[3]]
            bed_info_list.append(bed_info)
    new_path=bedpath+'_37'
    with open(new_path, 'w') as output_handle:
        writer=csv.writer(output_handle,delimiter='\t')
        writer.writerows(bed_info_list)
    return new_path

def obtain_SNPs(probe_list,bedpath,freq_threshold):
    #first obtain frequent SNPs in targets
    new_path=target_bed_to_bed_37(bedpath)
    os.system('tabix '+ftp_path_snp_database+' -R '+new_path+' > tmp_output_snp') #put this FTP URL in the config file?
    chr_list=[]
    loc_list=[]
    freq_threshold=1.0-freq_threshold
    with open('tmp_output_snp') as inputfile:
        input_reader=csv.reader(inputfile)
        for row in input_reader:
            if "FREQ" in row[0]:
                freq=float(row[0].split(";")[-1].split(":")[-1])
                if freq < freq_threshold:
                    chrom=row[0].split('\t')[0]
                    locus=int(row[0].split('\t')[1])
                    if chrom in chr_list:
                        loc_list[chr_list.index(chrom)].append(locus)
                    else:
                        chr_list.append(chrom)
                        loc_list.append([locus])
    os.system('rm tmp_output_snp')
    SNP_count=[[0 for probe in probe_arms] for probe_arms in probe_list]
    #first create a bedfile containing each of the arms
    arm1list=[]
    arm2list=[]
    for i,probe_arms in enumerate(probe_list):
        for j,probe in enumerate(probe_arms):
            arm1=loci_to_bed_37(probe[2])
            arm1.append(str(i)+'-'+str(j))
            arm2=loci_to_bed_37(probe[3])
            arm2.append(str(i)+'-'+str(j))
            if arm1[0] in chr_list:
                index_to_check=chr_list.index(arm1[0])
                for locus in loc_list[index_to_check]:
                    if locus in range(int(arm1[1]),int(arm1[2])+1): #check whether the SNP is in the arm
                        SNP_count[i][j]+=1
                    if locus in range(int(arm2[1]),int(arm2[2])+1): #check whether the SNP is in the arm
                        SNP_count[i][j]+=1
    return SNP_count
#----------------------------------------------------------------------------------------------------------

with open(filename) as handle:
    pool=mp.Pool(nr_of_cores)#mp.cpu_count())
    possible_arm_combinations_all_targets=[pool.apply(create_all_possible_arms_both_strands,args=(record,min_target_length,max_target_length,min_arm_length,max_arm_length,min_target_length,max_CG_percentage,cpg_flanks,mid_loc)) for record in SeqIO.parse(handle,"fasta")]
    pool.close()
print('All possible arms obtained')
pool=mp.Pool(nr_of_cores)
tms=[[pool.apply(get_delta_Tm, args=([probe_arms])) for probe_arms in possible_arm_combinations] for possible_arm_combinations in possible_arm_combinations_all_targets]
pool.close()
print('\tTms obtained')
pool=mp.Pool(nr_of_cores)
CpG_conflicts=[[pool.apply(report_CpGs_in_arms, args=([probe_arms])) for probe_arms in possible_arm_combinations] for possible_arm_combinations in possible_arm_combinations_all_targets]
pool.close()
print('\tCpGs in arms counted')
pool=mp.Pool(nr_of_cores)
possible_probes_all_targets=[[pool.apply(add_backbone, args=(probe_arms,backbone_sequence)) for probe_arms in possible_arm_combinations] for possible_arm_combinations in possible_arm_combinations_all_targets]
pool.close()
print('\tBackbone added')
fasta_name='temp_test_fasta'
outputname_json='temp_test_json'
pool=mp.Pool(nr_of_cores)
hairpin_scores=[pool.apply(check_probe_for_hairpin_score,args=(possible_probes,fasta_name,outputname_json,index)) for index,possible_probes in enumerate(possible_probes_all_targets)]
pool.close()
print('\tHairpins detected')

SNP_conflicts=obtain_SNPs(possible_arm_combinations_all_targets,bedfile,freq_threshold)
print('\tSNPs in arms counted')

#------------------------------------------------------ Store all information from the prossible probes into files
#with open(outputdir+'possible_arm_combinations_all_targets.csv', 'w') as file:
#    for item in possible_arm_combinations_all_targets:
#            file.write(",".join(map(str,item)))
#            file.write("\n")
with open(outputdir+'tms.csv', 'w') as file:
    for item in tms:
            file.write(",".join(map(str,item)))
            file.write("\n")
with open(outputdir+'CpG_conflicts.csv', 'w') as file:
    for item in CpG_conflicts:
            file.write(",".join(map(str,item)))
            file.write("\n")
with open(outputdir+'possible_probes_all_targets.csv', 'w') as file:
    for item in possible_probes_all_targets:
            file.write(",".join(map(str,item)))
            file.write("\n")
with open(outputdir+'hairpin_scores.csv', 'w') as file:
    for item in hairpin_scores:
            file.write(",".join(map(str,item)))
            file.write("\n")
with open(outputdir+'SNP_conflicts.csv', 'w') as file:
    for item in SNP_conflicts:
            file.write(",".join(map(str,item)))
            file.write("\n")
