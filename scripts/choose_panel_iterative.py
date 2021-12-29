#!/usr/bin/python3
# Create panels and choose the best
import os
from argparse import ArgumentParser
import json
from Bio import SeqIO, Seq, SeqRecord
import multiprocessing as mp
from numpy.random import choice
import csv
import numpy as np
import pickle


def main():
    parser = get_arg_parser()
    args = vars(parser.parse_args())
    output_name = args["output_name"]
    config_file = args["config_file"]
    tms_file = args['tms']
    cpgs = args['cpg']
    probes = args['probes']
    hairpins = args['hairpins']
    snps = args['snps']
    targets = args['targets']
    probe_arms_file = args['probe_arms_file']
    nr_of_cores = int(args["cores"])
    outputdir = args["output_dir"] + '/'
    conflicting_file = args["conflicting_file"]

    permutations, backbone_length, exclusion_factor, seed = import_config(config_file)

    probe_cpg_id_list, possible_probes_all_targets, tms, cpg_conflicts, snp_conflicts, hairpin_scores, probe_arm_list, \
        arm_upstream_loc_list, arm_downstream_loc_list, probe_id_list, arm_upstream_list, arm_downstream_list, \
        possible_arm_combinations_all_targets = import_probe_parameters(targets, probes, tms_file, cpgs, snps, hairpins,
                                                                        probe_arms_file)

    probe_costs = get_probe_costs_array(tms, cpg_conflicts, snp_conflicts, hairpin_scores)
    print('\tProbe costs computed')

    # Iteratively exclude the most dimer forming probes.
    chosen_probes_lists = []
    probes_with_dimers_lists = []
    conflicting_probe_list = []
    counter = 0
    min_dimers = len(possible_probes_all_targets)

    create_arms_loci_bedfile(outputdir, arm_upstream_loc_list, arm_downstream_loc_list, probe_arm_list)

    while counter < permutations and min_dimers > 0:
        # Choose random probe set
        # put the below lines in a function that can be called ~100 times to be able to find the best subset of probes.
        pool = mp.Pool(nr_of_cores)
        chosen_probes = [
            choose_probes_from_costs(probe_costs[i], possible_arm_combinations, permutations, counter, probe_id_list[i],
                                     seed) for i, possible_arm_combinations in enumerate(
                possible_probes_all_targets)]  
        # function that choses probes by the probability which is based on the cost
        pool.close()
        print('\tRound ' + str(counter) + ' :Probes chosen')

        chosen_probes_lists.append(chosen_probes)
        probes_with_dimers, conflicting_targets, probe_costs, cpgs_of_dimer_forming_probes = \
            increase_costs_dimer_forming_probes_iterative(possible_arm_combinations_all_targets, 
                                                          chosen_probes, probe_costs, probe_cpg_id_list, 
                                                          probe_id_list, arm_upstream_list, arm_downstream_list, 
                                                          arm_upstream_loc_list, arm_downstream_loc_list, 
                                                          exclusion_factor, outputdir, score_cutoff, tm_cutoff)

        nr_of_dimer_probes = len(cpgs_of_dimer_forming_probes)
        probes_with_dimers_lists.append(nr_of_dimer_probes)
        conflicting_probe_list.append(probes_with_dimers)
        print('\t' + str(nr_of_dimer_probes) + ' out of ' + str(len(chosen_probes)) + ' probes form a dimer')
        min_dimers = min(nr_of_dimer_probes, min_dimers)
        output_nested_list(outputdir + 'costs_iter' + str(counter), probe_costs)
        # output_nested_list(outputdir+'chosen_probes'+str(counter),chosen_probes)
        counter += 1

    chosen_set_index = probes_with_dimers_lists.index(
        min(probes_with_dimers_lists))  # Also check the sum of all probe scores?
    chosen_set = chosen_probes_lists[chosen_set_index]

    write_output(targets, output_name, chosen_set, conflicting_file, conflicting_probe_list, min_dimers)


def get_arg_parser():
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
    return parser


def import_config(config_file):
    with open(config_file) as jsonFile:
        config_object = json.load(jsonFile)
        jsonFile.close()
    permutations = int(config_object['permutations'])
    backbone_length = len(config_object["backbone_sequence"][0]["reverse_complement_universal_forward_primer"]) + len(
        config_object["backbone_sequence"][0]["common_region"]) + len(
        config_object["backbone_sequence"][0]["universal_reverse_primer"])
    exclusion_factor = float(config_object["dimer_exclusion_factor"])
    seed = int(config_object['seed'])
    score_cutoff = int(config_object['mfeprimer_dimer_parameters'][0]['score_cutoff'])
    tm_cutoff = int(config_object['mfeprimer_dimer_parameters'][0]['tm_cutoff'])
    return permutations, backbone_length, exclusion_factor, seed, score_cutoff, tm_cutoff


def import_probe_parameters(targets, probes, tms_file, cpgs, snps, hairpins, probe_arms_file):
    with open(targets, 'r') as handle:
        reader = csv.reader(handle, delimiter='\t')
        probe_cpg_id_list = []
        for row in reader:
            probe_cpg_id_list.append(row[3])

    with open(probes, 'r') as handle:
        reader = csv.reader(handle)
        possible_probes_all_targets = []
        for i, row in enumerate(reader):
            possible_probes_all_targets.append([])
            for probe in row:
                possible_probes_all_targets[i].append(probe)

    with open(tms_file, 'rb') as handle:
        tms = pickle.load(handle)
        tms = np.array(tms, dtype=object)

    with open(cpgs, 'r') as handle:
        reader = csv.reader(handle)
        cpg_conflicts = []
        for i, row in enumerate(reader):
            cpg_conflict_list = []
            for probe in row:
                cpg_conflict_list.append(int(probe))
            cpg_conflicts.append(np.array(cpg_conflict_list, dtype=object))
        cpg_conflicts = np.array(cpg_conflicts, dtype=object)

    with open(snps, 'r') as handle:
        reader = csv.reader(handle)
        snp_conflicts = []
        for i, row in enumerate(reader):
            snp_conflict_list = []
            for probe in row:
                snp_conflict_list.append(int(probe))
            snp_conflicts.append(np.array(snp_conflict_list, dtype=object))
        snp_conflicts = np.array(snp_conflicts, dtype=object)

    with open(hairpins, 'r') as handle:
        reader = csv.reader(handle)
        hairpin_scores = []
        for i, row in enumerate(reader):
            hairpin_score_list = []
            for probe in row:
                hairpin_score_list.append(int(probe))
            hairpin_scores.append(np.array(hairpin_score_list, dtype=object))
        hairpin_scores = np.array(hairpin_scores, dtype=object)

    with open(probe_arms_file, 'rb') as handle:
        possible_arm_combinations_all_targets = pickle.load(handle)
        probe_arm_list = []
        arm_upstream_loc_list = []
        arm_downstream_loc_list = []
        probe_id_list = []
        arm_upstream_list = []
        arm_downstream_list = []
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
                probe_id_list[i].append(str(probe[7]) + ':' + str(probe[4]))
                arm_upstream_list[i].append(probe[0])
                arm_downstream_list[i].append(probe[1])
    return probe_cpg_id_list, possible_probes_all_targets, tms, cpg_conflicts, snp_conflicts, hairpin_scores, \
        probe_arm_list, arm_upstream_loc_list, arm_downstream_loc_list, probe_id_list, arm_upstream_list, \
        arm_downstream_list, possible_arm_combinations_all_targets


def get_probe_costs_array(tms, cpg_conflicts, snp_conflicts, hairpin_array):
    probe_cost = np.add(np.add(np.add(np.add(np.multiply(hairpin_array, int(10**10)), cpg_conflicts), snp_conflicts), tms),1)
    # Score is at least 1
    return probe_cost


def output_nested_list(filename, nested_list):
    with open(filename, 'w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerows(nested_list)
    return


def choose_probes_from_costs(probe_costs, possible_arm_combinations, n, counter, probe_id_list, seed):
    np.random.seed(seed)
    if probe_costs.size == 0:
        return
    else:
        probe_scores = 1 / probe_costs.astype(np.float64)
        sum_score = sum(probe_scores)
        if sum_score == 0:  # If all scores are 0, there should still be a probability to choose one of the probes
            probability_distribution = [1 / len(probe_scores)] * len(probe_scores)
        else:
            probability_distribution = np.divide(probe_scores, sum_score)
        probe = choice(possible_arm_combinations, n, p=probability_distribution)
        index = possible_arm_combinations.index(probe[counter])
    return [probe[counter], probe_id_list[index]]


def get_conflict_ranges(conflicting_cpg_list_dimer, conflict_range_dimer, exclusion_factor):
    probe_list = []
    probe_list_count = []
    for probe in conflicting_cpg_list_dimer:  # only look into the first mentioned probes when forming a dimer
        if probe in probe_list:
            probe_list_count[probe_list.index(probe)] += 1
        else:
            probe_list.append(probe)
            probe_list_count.append(1)
    if probe_list_count == []:
        probe_list_count = [0]
    nr_of_times = max(probe_list_count)
    sorted_nr_list = []
    sorted_probe_list = []
    while nr_of_times > 0:
        for i, nr in enumerate(probe_list_count):
            if nr == nr_of_times:
                sorted_probe_list.append(probe_list[i])
                sorted_nr_list.append(nr)
        nr_of_times = nr_of_times - 1
    if exclusion_factor == float(2):
        amount_to_exlude = 1
    else:
        amount_to_exlude = int(int(len(sorted_probe_list) * exclusion_factor))
    if amount_to_exlude == 0:
        most_conflicting_probe_list = sorted_probe_list
    else:
        most_conflicting_probe_list = sorted_probe_list[
                                      0:amount_to_exlude]  # exclude the top portion (top 10% if exclusion factor== 0.1
    most_conflicting_dimer_range_list = []
    for probe in most_conflicting_probe_list:  # report the conflicting locus at first occurrence of dimer
        index = conflicting_cpg_list_dimer.index(probe)
        most_conflicting_dimer_range_list.append(conflict_range_dimer[index])
    return most_conflicting_dimer_range_list


def get_conflict_ranges_s1_and_s2(conflicting_cpg_list_dimer, conflict_range_dimer, conflicting_cpg_list_dimer2,
                                  conflict_range_dimer2, exclusion_factor):
    combined_conflicting_cpg_list_dimer = conflicting_cpg_list_dimer + conflicting_cpg_list_dimer2
    combined_conflict_range_dimer = conflict_range_dimer + conflict_range_dimer2
    combined_most_conflicting_dimer_range_list = get_conflict_ranges(combined_conflicting_cpg_list_dimer,
                                                                     combined_conflict_range_dimer, exclusion_factor)
    return combined_most_conflicting_dimer_range_list


def get_dimer_range_upstream(possible_arm_combinations_all_targets, arm_upstream_loc_list, probe_index, nested_index,
                             dimer_start_num, dimer_end_num, dimer):
    chrom = arm_upstream_loc_list[probe_index][nested_index].split(':')[0]
    start = int(arm_upstream_loc_list[probe_index][nested_index].split(':')[1].split('-')[0].strip("'"))
    end = int(arm_upstream_loc_list[probe_index][nested_index].split(':')[1].split('-')[1].strip("'"))
    template_strand = possible_arm_combinations_all_targets[probe_index][nested_index][-2]
    dimer_start_range = max(start, start + dimer_start_num)
    if template_strand == '+':
        dimer_end_range = min(end, start + dimer_end_num)
    else:
        dimer_end_range = max(end, start + dimer_end_num)
    dimer_range = str(chrom) + ':' + str(dimer_start_range) + '-' + str(dimer_end_range)
    if dimer_start_range > dimer_end_range:
        print('error ' + str([dimer_start_range, dimer_end_range, dimer]))
    return [dimer_range, template_strand]


def get_dimer_range_downstream(possible_arm_combinations_all_targets, arm_downstream_loc_list, probe_index,
                               nested_index, dimer_start_num, dimer_end_num, s_length, dimer):
    chrom = arm_downstream_loc_list[probe_index][nested_index].split(':')[0]
    start = int(arm_downstream_loc_list[probe_index][nested_index].split(':')[1].split('-')[0].strip("'"))
    end = int(arm_downstream_loc_list[probe_index][nested_index].split(':')[1].split('-')[1].strip("'"))
    template_strand = possible_arm_combinations_all_targets[probe_index][nested_index][-2]
    dimer_start_range = max(start, start + dimer_start_num - (s_length - (end - start + 1)))
    if template_strand == '+':
        dimer_end_range = min(end, start + dimer_end_num - (s_length - (end - start + 1)))
    else:
        dimer_end_range = max(end, start + dimer_end_num - (s_length - (end - start + 1)))
    dimer_range = str(chrom) + ':' + str(dimer_start_range) + '-' + str(dimer_end_range)
    if dimer_start_range > dimer_end_range:
        print('error ' + str([dimer_start_range, dimer_end_range, dimer]))
    return [dimer_range, template_strand]


def write_nested_loci_to_bedfile(output_name, loci_list_up, loci_list_down, probe_arm_list):
    to_write_list_up = [
        [locus.split(':')[0], locus.split(':')[1].split('-')[0], locus.split('-')[1], '', '', probe_arm_list[i][j][-2],
         str(i) + ':' + str(j)] for i, row in enumerate(loci_list_up) for j, locus in enumerate(row)]  
    # [chrom,start,end,i:j] e.g. [1,105386,105396,0:304] for chr1:105386-105396 in the locus at loci_list[0][304]
    to_write_list_down = [
        [locus.split(':')[0], locus.split(':')[1].split('-')[0], locus.split('-')[1], '', '', probe_arm_list[i][j][-2],
         str(i) + ':' + str(j)] for i, row in enumerate(loci_list_down) for j, locus in enumerate(row)]  
    # [chrom,start,end,    i:j] e.g. [1,105386,105396,0:304] for chr1:105386-105396 in the locus at loci_list[0][304]
    with open(output_name, 'w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerows(to_write_list_up)
        writer.writerows(to_write_list_down)
    return output_name


def write_loci_to_bedfile(output_name, loci_list):
    to_write_list = [
        [locus[0].split(':')[0], locus[0].split(':')[1].split('-')[0], locus[0].split('-')[1], '', '', locus[1], str(i)]
        for i, locus in enumerate(loci_list)]  
    # [chrom,start,end,i] e.g. [1,105386,105396,0] for chr1:105386-105396 in the locus at loci_list[0]
    with open(output_name, 'w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerows(to_write_list)
    return output_name


def create_arms_loci_bedfile(outputdir, arm_upstream_loc_list, arm_downstream_loc_list, probe_arm_list):
    loci_bedfile_name = outputdir + 'combined.bed'
    write_nested_loci_to_bedfile(loci_bedfile_name, arm_upstream_loc_list, arm_downstream_loc_list, probe_arm_list)
    return


def create_conflicting_indices_list_bedtools(dimer_range_list, dimer_bedfile_name, bedfile_intersect_name, outputdir):
    write_loci_to_bedfile(dimer_bedfile_name, dimer_range_list)
    # Do bedtools intersect here on combined_loci_bedfile_name and dimer_bedfile_name
    os.system('bedtools intersect -wa -s -a ' + outputdir + 'combined.bed -b ' + dimer_bedfile_name + ' -f 7E-9 > ' 
              + bedfile_intersect_name)
    # obtain identifier from bedfile_intersect_name
    new_conflicting_indices_list = []
    with open(bedfile_intersect_name, 'r') as handle:
        reader = csv.reader(handle, delimiter='\t')
        for row in reader:
            identifier = row[-1]
            i = identifier.split(':')[0]
            j = identifier.split(':')[1]
            new_conflicting_indices_list.append([i, j])
    return new_conflicting_indices_list


def increase_costs_dimer_forming_probes_iterative(possible_arm_combinations_all_targets, chosen_probes, probe_costs,
                                                  probe_cpg_id_list, probe_id_list, arm_upstream_list,
                                                  arm_downstream_list, arm_upstream_loc_list, arm_downstream_loc_list,
                                                  exclusion_factor, outputdir, score_cutoff, tm_cutoff):
    # Create a fasta file with the chosen probes in this iteration.
    passed_list = []
    with open('tmp_fasta_chosen_probes.fasta', 'w') as handle:
        for i, probe in enumerate(chosen_probes):
            if probe is None:
                passed_list.append(i)
            else:
                handle.write('>' + probe[1] + '\n')
                handle.write(probe[0] + '\n')
    # test the chosen probes set for dimers
    os.system(
        'mfeprimer dimer -i tmp_fasta_chosen_probes.fasta -j -o tmp_dimers_chosen_probes -s '+str(score_cutoff)+' -t '+str(tm_cutoff))
    # Score cut-off and temperature cut-off are set in configuration file

    # Obtain dimer list from mfeprimer output
    os.system('rm tmp_fasta_chosen_probes.fasta tmp_dimers_chosen_probes')
    dimer_file = 'tmp_dimers_chosen_probes.json'
    with open(dimer_file) as jsonFile:
        dimerlist = json.load(jsonFile)
        jsonFile.close()
    conflicting_cpg_list_dimer = []
    conflicting_cpg_list_dimer2 = []
    conflict_range_dimer = []
    conflict_range_dimer2 = []
    cpgs_of_dimer_forming_probes = set([])
    if dimerlist is None:
        pass
    else:
        for dindex, dimer in enumerate(dimerlist):
            cpg_id = dimer['S1']['ID']
            cpg_id2 = dimer['S2']['ID']
            cpgs_of_dimer_forming_probes.add(cpg_id)
            cpgs_of_dimer_forming_probes.add(cpg_id2)
            dimer_structure = dimer['Aseq']
            s1_length = len(dimer['S1']['Seq'])
            s2_length = len(dimer['S2']['Seq'])
            # Obtain the location of the dimer on the two sequences
            s_1_dangling = dimer['S1Dangling']
            a_seq_dangling = dimer['AseqDangling']
            s_2_dangling = dimer['S2Dangling']
            dimer_start_num = a_seq_dangling - s_1_dangling
            dimer_end_num = dimer_start_num + len(dimer_structure)
            dimer_start_num2 = a_seq_dangling - s_2_dangling
            dimer_end_num2 = dimer_start_num2 + len(dimer_structure)
            conflicting_cpg_list_dimer.append(cpg_id)
            conflicting_cpg_list_dimer2.append(cpg_id2)
            probe_index = probe_cpg_id_list.index(cpg_id.split(':')[0])
            nested_index = probe_id_list[probe_index].index(cpg_id)
            probe_index2 = probe_cpg_id_list.index(cpg_id2.split(':')[0])
            nested_index2 = probe_id_list[probe_index2].index(cpg_id2)
            # obtain the genomic locations of the dimer forming regions for S1
            # check whether the S1 and S2 sequences are 5' to 3' or the other way around.)
            if arm_upstream_list[probe_index][nested_index].strip("'") in dimer['S1']['Seq'][:len(
                    arm_upstream_list[probe_index][nested_index].strip("'"))]:
                if dimer_start_num < len(
                        arm_upstream_list[probe_index][nested_index].strip("'")):  # dimer is in upstream arm
                    dimer_range = get_dimer_range_upstream(possible_arm_combinations_all_targets, arm_upstream_loc_list,
                                                           probe_index, nested_index, dimer_start_num, dimer_end_num,
                                                           dimer)
                elif dimer_end_num >= s1_length - len(arm_downstream_list[probe_index][nested_index].strip(
                        "'")):  # dimer is in downstream arm.
                    # If dimer end has no overlap with the downsteam or upstream arm is should be in the backbone.
                    dimer_range = get_dimer_range_downstream(possible_arm_combinations_all_targets,
                                                             arm_downstream_loc_list, probe_index, nested_index,
                                                             dimer_start_num, dimer_end_num, s1_length, dimer)
                else:  # dimer forming region is in backbone,
                    # only report the other probe that forms the dimer together with this probe.
                    pass
            else:  # In this situation the S1 sequence is 3'-5', therefore we need to correct the indices.
                dimer_start_num, dimer_end_num = s1_length - dimer_end_num, s1_length - dimer_start_num
                if dimer_start_num < len(
                        arm_upstream_list[probe_index][nested_index].strip("'")):  # dimer is in upstream arm
                    dimer_range = get_dimer_range_upstream(possible_arm_combinations_all_targets, arm_upstream_loc_list,
                                                           probe_index, nested_index, dimer_start_num, dimer_end_num,
                                                           dimer)
                elif dimer_end_num >= s1_length - len(arm_downstream_list[probe_index][nested_index].strip(
                        "'")):  # dimer is in downstream arm
                    # If dimer end has no overlap with the downsteam or upstream arm is should be in the backbone.
                    dimer_range = get_dimer_range_downstream(possible_arm_combinations_all_targets,
                                                             arm_downstream_loc_list, probe_index, nested_index,
                                                             dimer_start_num, dimer_end_num, s1_length, dimer)
                else:  # dimer forming region is in backbone,
                    # only report the other probe that forms the dimer together with this probe.
                    pass
            if arm_upstream_list[probe_index2][nested_index2].strip("'") in dimer['S2']['Seq'][:len(
                    arm_upstream_list[probe_index2][nested_index2].strip("'"))]:
                if dimer_start_num2 < len(
                        arm_upstream_list[probe_index2][nested_index2].strip("'")):  # dimer is in upstream arm
                    dimer_range2 = get_dimer_range_upstream(possible_arm_combinations_all_targets,
                                                            arm_upstream_loc_list, probe_index2, nested_index2,
                                                            dimer_start_num2, dimer_end_num2, dimer)
                elif dimer_end_num2 >= s2_length - len(arm_downstream_list[probe_index2][nested_index2].strip(
                        "'")):  # dimer is in downstream arm.
                    # If dimer end has no overlap with the downsteam or upstream arm is should be in the backbone.
                    dimer_range2 = get_dimer_range_downstream(possible_arm_combinations_all_targets,
                                                              arm_downstream_loc_list, probe_index2, nested_index2,
                                                              dimer_start_num2, dimer_end_num2, s2_length, dimer)
                else:  # dimer forming region is in backbone,
                    # only report the other probe that forms the dimer together with this probe.
                    pass
            else:  # In this situation the S2 sequence is 3'-5', therefore we need to correct the indices.
                dimer_start_num2, dimer_end_num2 = s2_length - dimer_end_num2, s2_length - dimer_start_num2
                if dimer_start_num2 < len(
                        arm_upstream_list[probe_index2][nested_index2].strip("'")):  # dimer is in upstream arm
                    dimer_range2 = get_dimer_range_upstream(possible_arm_combinations_all_targets,
                                                            arm_upstream_loc_list, probe_index2, nested_index2,
                                                            dimer_start_num2, dimer_end_num2, dimer)
                elif dimer_end_num2 >= s2_length - len(arm_downstream_list[probe_index2][nested_index2].strip(
                        "'")):  # dimer is in downstream arm.
                    # If dimer end has no overlap with the downsteam or upstream arm is should be in the backbone.
                    dimer_range2 = get_dimer_range_downstream(possible_arm_combinations_all_targets,
                                                              arm_downstream_loc_list, probe_index2, nested_index2,
                                                              dimer_start_num2, dimer_end_num2, s2_length, dimer)
                else:  # dimer forming region is in backbone,
                    # only report the other probe that forms the dimer together with this probe.
                    pass
            conflict_range_dimer.append(dimer_range)
            conflict_range_dimer2.append(dimer_range2)

    # Only report the probe that formed the most amount of dimers
    most_conflicting_dimer_range_list = get_conflict_ranges_s1_and_s2(conflicting_cpg_list_dimer, conflict_range_dimer,
                                                                      conflicting_cpg_list_dimer2,
                                                                      conflict_range_dimer2,
                                                                      exclusion_factor)
    # obtain conflict ranges from dimer forming probes

    # find all probes that have common regions with the most_conflicting_dimer_range and report as conflicting probes
    all_conflicting_probe_list = []
    conflicting_targets = set([])
    print('Dimer ranges obtained')
    dimer_bedfile_name = outputdir + 'conflicting_dimer.bed'
    bedfile_intersect_name = outputdir + 'conflicting_loci.bed'
    new_conflicting_indices_list = create_conflicting_indices_list_bedtools(most_conflicting_dimer_range_list,
                                                                            dimer_bedfile_name, bedfile_intersect_name,
                                                                            outputdir)
    for indexes in new_conflicting_indices_list:
        i = int(indexes[0])
        j = int(indexes[1])
        all_conflicting_probe_list.append(probe_id_list[i][j])
        conflicting_targets.add(probe_id_list[i][0].split(':')[0])
        probe_costs[i][j] = probe_costs[i][j] + float(10 ** 10)  # Adjust cost of conflicting probe
    print(str(len(conflicting_targets)) + ' conflicts found')
    return all_conflicting_probe_list, conflicting_targets, probe_costs, cpgs_of_dimer_forming_probes


def write_output(targets, output_name, chosen_set, conflicting_file, conflicting_probe_list, min_dimers):
    seq_list = []
    print('minimal number of dimers is ' + str(min_dimers))
    with open(targets, 'r') as handle:
        reader = csv.reader(handle, delimiter='\t')
        description_list = []
        for row in reader:
            loc=int((int(row[1])+int(row[2])-1)/2)  #Obtain middle value from target bed file. Middle value is the location of the CpG.
            description_list.append(str(row[0] + ':' + str(loc)  + '\t' + row[-1]))

    with open(output_name, 'w') as handle:
        for i, probe in enumerate(chosen_set):
            if probe is None:
                pass
            else:
                probename = str(i)
                probe_description = description_list[i]
                seq_list.append(SeqRecord.SeqRecord(Seq.Seq(probe[0]), id=probename, description=probe_description))
                # probe[1]))
        SeqIO.write(seq_list, handle, 'fasta')
    with open(conflicting_file, 'w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        for row in conflicting_probe_list:
            writer.writerow(row)


if __name__ == '__main__':
    main()
