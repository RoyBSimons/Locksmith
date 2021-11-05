#!/usr/bin/python3
#Create panels and choose the best
import os
from argparse import ArgumentParser
import json
from Bio import SeqIO, Seq, SeqRecord
import multiprocessing as mp
from numpy.random import choice
import csv

parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="filename",
                    help="open FILE", metavar="FILE")
parser.add_argument("-o", "--output", dest="output_dir",
                    help="Directory for output", metavar="OUTPUTDIR")
parser.add_argument("-c", "--config_file", dest="config_file",
                    help="Config file with probe specifics", metavar="JSONfile")
parser.add_argument("-m", "--tms", dest="tms",
                    help="csv file containing probe binding temperatures", metavar="TMS")
parser.add_argument("-g", "--cpg", dest="cpg",
                    help="csv file containing amount of cpgs in probe arms", metavar="CPG")
parser.add_argument("-p", "--probes", dest="probes",
                    help="csv file containing the full probes", metavar="PROBES")
parser.add_argument("-h", "--hairpins", dest="hairpins",
                    help="csv file containing hairpins of probes", metavar="HAIRPINS")
parser.add_argument("-n", "--snp", dest="snps",
                    help="csv file containing amount of snps in probe arms", metavar="SNPS")
parser.add_argument("-s", "--score", dest="scores",
                    help="csv file containing scores of probes", metavar="SCORES")
parser.add_argument("-t", "--cores", dest="cores",
                    help="Passed amount of cores to script", metavar="Cores")

args = vars(parser.parse_args())

filename=args["filename"]
outputdir=args["output_dir"]+'/'
config_file=args["config_file"]
tms=args['tms']
cpgs=args['cpg']
probes=args['probes']
hairpins=args['hairpins']
snps=args['snps']
scores=args['scores']
nr_of_cores=int(args["cores"])
with open(config_file) as jsonFile:
    configObject = json.load(jsonFile)
    jsonFile.close()
permutations=int(configObject['Permutations'])

with open(probes, 'r') as handle:
    reader=csv.reader(handle)
    probe_list=[]
    for i,row in enumerate(reader):
        probe_list.append(row[0])
        for probe in row[1:]:
            probe_list[i].append(probe)
with open(scores, 'r') as handle:
    reader=csv.reader(handle)
    scores_list=[]
    for i,row in enumerate(reader):
        scores_list.append(row[0])
        for score in row[1:]:
            scores_list[i].append(score)
#--------------------------------------------------------
def choose_probes_from_scores(probe_scores,possible_arm_combinations,n,counter):
    if possible_arm_combinations==[]:
        return
    else:
        sum_score=sum(probe_scores)
        probability_distribution=[i/sum_score for i in probe_scores]
        probe= choice(possible_arm_combinations,n,p=probability_distribution)
    return probe[counter]

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



#------------------------------------------------------
counter=0
chosen_probes_lists=[]
probes_with_dimers_lists=[]
while counter < permutations:
    # put the below lines in a function that can be called ~100 times to be able to find the best subset of probes.
    pool=mp.Pool(nr_of_cores)
    chosen_probes=[]
    chosen_probes=[pool.apply(choose_probes_from_scores,args=(probe_scores[i],possible_arm_combinations,permutations,counter)) for i,possible_arm_combinations in enumerate(possible_probes_all_targets)] #function that choses probes by the probability which is based on the score
    pool.close()
    print('\tProbes chosen')
    chosen_probes_lists.append(chosen_probes)
    probes_with_dimers=find_dimer_forming_probes(chosen_probes) #function that returns for each chosen probe whether if forms a dimer with one of the other chosen probes
    nr_of_dimer_probes=sum(probes_with_dimers)
    probes_with_dimers_lists.append(nr_of_dimer_probes)
    print('\t'+str(nr_of_dimer_probes)+' probes form a dimer')
    # put the above lines in a function that can be called ~100 times to be able to find the best subset of probes.
    counter+=1
chosen_set_index=probes_with_dimers_lists.index(min(probes_with_dimers_lists)) #Also check the sum of all probe scores?
chosen_set=chosen_probes_lists[chosen_set_index]
print(chosen_set)
seq_list=[]
with open(outputdir+'chosen_probes.fasta','w') as handle:
    for i,probe in enumerate(chosen_set):
        if probe == None:
            pass
        else:
            probename=str(i)
            seq_list.append(SeqRecord.SeqRecord(Seq.Seq(probe),id=probename))
    SeqIO.write(seq_list,handle,'fasta')
#print(probes_with_dimers_lists)
#print(chosen_probes_lists)
#to write
#    get_probe_scores
#    choose_probes_from_scores
#    find_dimer_forming_probes

#test print first probe. First arm combinations from first target site.
#print(tms)
#print(possible_arm_combinations_all_targets[0][0])
#print(CpG_conflicts)
#print(possible_probes_all_targets[0][0])
#print(hairpin_scores)
#print(SNP_conflicts)

