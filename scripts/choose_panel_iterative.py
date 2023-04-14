#!/usr/bin/python3
# Create panels and choose the best
import os
from argparse import ArgumentParser
import json
from Bio import SeqIO, Seq, SeqRecord
from Bio.SeqUtils import GC
import multiprocessing as mp
from numpy.random import choice
import csv
import numpy as np
import pickle
import statistics


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
    panel_output_file = args['panel_output_file']

    permutations, backbone_sequence, backbone_length, exclusion_factor, seed, score_cutoff, tm_cutoff, \
            scoring_weights, max_delta_tm_panel = import_config(config_file)  # Import parameters from configuration file.

    probe_cpg_id_list, possible_probes_all_targets, tms, cpg_conflicts, snp_conflicts, hairpin_scores, probe_arm_list, \
        arm_upstream_loc_list, arm_downstream_loc_list, probe_id_list, arm_upstream_list, arm_downstream_list, \
        possible_arm_combinations_all_targets, tms_up, tms_down = import_probe_parameters(targets, probes, tms_file, 
                                                                                            cpgs, snps, hairpins,
                                                                                            probe_arms_file)  # Import output create probe rule.

    probe_costs = get_probe_costs_array(tms, cpg_conflicts, snp_conflicts, hairpin_scores, scoring_weights)  # Obtain probe costs from delta Tm, hairpin formation, CpG and SNP conflicts.
    print('\tProbe costs computed')

    # Iteratively exclude the most dimer forming probes.
    chosen_probes_lists = []
    probes_with_dimers_lists = []
    conflicting_probe_list = []
    tm_amount_out_of_range_list = []
    counter = 0
    min_dimers = len(possible_probes_all_targets)
    tm_amount_out_of_range = 1

    # Create a bedfile for all possible upstream and downstream arms.
    # This file is used for the create_conflicting_indices_list_bedtools function inside the increase_costs_dimer_forming_probes_iterative function.
    loci_bedfile_name = outputdir + 'possible_probe_info/combined.bed'
    write_nested_loci_to_bedfile(loci_bedfile_name, arm_upstream_loc_list, arm_downstream_loc_list, probe_arm_list)

    #Loop to choose probes and rescore dimer-forming probes untill the maximum amount of iterations is encountered or there are no dimer-forming probes in the panel.
    while counter < permutations and (min_dimers > 0 or tm_amount_out_of_range > 0):
        # Choose probe set based on probability which is set by cost of a probe.
        chosen_probes = [choose_probes_from_costs(probe_costs[i], possible_arm_combinations,
                                                    permutations, counter, probe_id_list[i],
                                                    seed)
                        for i, possible_arm_combinations in enumerate(possible_probes_all_targets)]

        print('\tRound ' + str(counter) + ' :Probes chosen')

        chosen_probes_lists.append(chosen_probes) #Keep panel of probes for later.

        # Increase cost of probes that are chosen as conflicting probes.
        # The dimer exclusion factor determines the percentage of dimer-forming probes that are reported as conflicting.
        probes_with_dimers, conflicting_targets, probe_costs, cpgs_of_dimer_forming_probes = \
            increase_costs_dimer_forming_probes_iterative(possible_arm_combinations_all_targets, 
                                                          chosen_probes, probe_costs, probe_cpg_id_list, 
                                                          probe_id_list, arm_upstream_list, arm_downstream_list, 
                                                          arm_upstream_loc_list, arm_downstream_loc_list, 
                                                          exclusion_factor, outputdir, score_cutoff, tm_cutoff, nr_of_cores)

        # Increase cost of all probes (also not-chosen probes) that do not fall into the panel Tm range.
        # Tm is out of bounds when it is outside of: median Tm + or - half of the Tm range threshold.
        probe_costs, tm_amount_out_of_range = increase_cost_probes_with_out_of_bounds_tm(chosen_probes,
                probe_id_list, probe_costs, tms_down, tms_up, scoring_weights, max_delta_tm_panel)
        tm_amount_out_of_range_list.append(tm_amount_out_of_range)
        nr_of_dimer_probes = len(cpgs_of_dimer_forming_probes)
        probes_with_dimers_lists.append(nr_of_dimer_probes)
        conflicting_probe_list.append(probes_with_dimers)
        print('\t\t' + str(nr_of_dimer_probes) + ' out of ' + str(len(chosen_probes)) + ' probes form a dimer')
        print('\t\t' + str(tm_amount_out_of_range) + ' out of ' + str(len(chosen_probes)*2)+ ' annealing sites are outside the Tm range')
        min_dimers = min(nr_of_dimer_probes, min_dimers)  # Store the minimum number of dimers, to choose this panel after all iterations.
        counter += 1

    chosen_set_index_list = [ i for i,val in enumerate(probes_with_dimers_lists) if val == min_dimers ]  # Find index of panel with the least amount of dimers.
    if len(chosen_set_index_list) == 1:
        chosen_set_index = chosen_set_index_list[0]
        min_tm_chosen_panel = tm_amount_out_of_range_list[chosen_set_index]
    else:  # There are multiple panels with the same minimum amount of dimers. Keep the one with the least amount of probes out of the Tm range.
         min_tm_low_dimer_panels = [ tm_amount_out_of_range_list[index] for index in chosen_set_index_list]
         min_tm_chosen_panel = min(min_tm_low_dimer_panels)
         chosen_set_index = chosen_set_index_list[min_tm_low_dimer_panels.index(min_tm_chosen_panel)]
    chosen_set = chosen_probes_lists[chosen_set_index]  # Get panel with the least amount of dimers.
    print('\nStatistics of probe panel created by Locksmith:')
    print('\tMinimal number of dimers is ' + str(min_dimers))
    print('\tAmount of annealing sites outside of Tm range is ' + str(min_tm_chosen_panel))

    dimer_scores = get_dimer_scores(chosen_set, score_cutoff, tm_cutoff, targets, nr_of_cores)

    write_output(targets, output_name, chosen_set, conflicting_file, conflicting_probe_list, min_dimers,
                panel_output_file, probe_arm_list, backbone_sequence, tms, cpg_conflicts, snp_conflicts,
                hairpin_scores, probe_id_list, dimer_scores, tms_up, tms_down)  # Write the output to a file.
    os.system('rm ' + outputdir + '/possible_probe_info/combined.bed')

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
    parser.add_argument("-e", "--panel_output_file", dest="panel_output_file",
                        help="Comma separated file with information of all probes in the panel", metavar="PANEL")
    return parser


def import_config(config_file):
    # Import parameters set in configuration file
    with open(config_file) as jsonFile:
        config_object = json.load(jsonFile)
        jsonFile.close()
    permutations = int(config_object['permutations'])
    backbone_sequence = config_object["backbone_sequence"][0]["reverse_complement_universal_forward_primer"] \
        + config_object["backbone_sequence"][0]["common_region"] \
        + config_object["backbone_sequence"][0]["universal_reverse_primer"]
    backbone_length = len(backbone_sequence)
    exclusion_factor = float(config_object["dimer_exclusion_factor"])
    seed = int(config_object['seed'])
    score_cutoff = int(config_object['mfeprimer_dimer_parameters'][0]['score_cutoff'])
    tm_cutoff = int(config_object['mfeprimer_dimer_parameters'][0]['tm_cutoff'])
    scoring_weights = config_object['scoring_weights'][0]
    max_delta_tm_panel = int(config_object["probe_specifics"][0]["max_delta_tm_panel"])
    return permutations, backbone_sequence, backbone_length, exclusion_factor, seed, score_cutoff, tm_cutoff, scoring_weights, max_delta_tm_panel


def loadall(f):
    while True:
        try:
            yield pickle.load(f)
        except EOFError:
            break

def import_probe_parameters(targets, probes, tms_file, cpgs, snps, hairpins, probe_arms_file):
    # Import the probe parameters from each of the separate files created by the Create_probes rule
    # All 2D array are a object dtype as the amount of probes per CpG differ.

    # Get CpG id from bedfile containing all targets
    with open(targets, 'r') as handle:
        reader = csv.reader(handle, delimiter='\t')
        probe_cpg_id_list = []
        for row in reader:
            probe_cpg_id_list.append(row[3])

    # 2D array containing the padlock probes. Each row consist of an array of strings containing the full
    # padlock probe sequences.
    with open(probes, 'r') as handle:
        reader = csv.reader(handle)
        possible_probes_all_targets = []
        for i, row in enumerate(reader):
            possible_probes_all_targets.append([])
            for probe in row:
                possible_probes_all_targets[i].append(probe)

    # Obtain 2D array containing the difference in Tm between the upstream and downstream arm of the padlock probe.
    # Each row consist of the delta-tms for all possible padlock probes to target one CpG.
    with open(tms_file, 'rb') as handle:
        tms = loadall(handle)
        tms = np.array([np.array(tms_array,dtype=object) for tms_array in tms],dtype=object)
    tms_up_file = tms_file[:-7]+'_up'+tms_file[-7:]
    with open(tms_up_file, 'rb') as handle:
        tms_up = loadall(handle)
        tms_up = np.array([np.array(tms_array,dtype=object) for tms_array in tms_up],dtype=object)
    tms_down_file = tms_file[:-7]+'_down'+tms_file[-7:]
    with open(tms_down_file, 'rb') as handle:
        tms_down = loadall(handle)
        tms_down = np.array([np.array(tms_array,dtype=object) for tms_array in tms_down],dtype=object)
    # Obtain a 2D array containing the amount of CpGs in the arms of the padlock probe.
    # Each row consist of the amount of CpGs in the arms for all possible padlock probes to target one CpG.
    with open(cpgs, 'r') as handle:
        reader = csv.reader(handle)
        cpg_conflicts = []
        for i, row in enumerate(reader):
            cpg_conflict_list = []
            for probe in row:
                cpg_conflict_list.append(int(probe))
            cpg_conflicts.append(np.array(cpg_conflict_list, dtype=object))
        cpg_conflicts = np.array(cpg_conflicts, dtype=object)

    # Obtain a 2D array containing the amount of frequent SNPs in the arms of the padlock probe.
    # Each row consist of the amount of SNPs in the arms for all possible padlock probes to target one CpG.
    with open(snps, 'r') as handle:
        reader = csv.reader(handle)
        snp_conflicts = []
        for i, row in enumerate(reader):
            snp_conflict_list = []
            for probe in row:
                snp_conflict_list.append(int(probe))
            snp_conflicts.append(np.array(snp_conflict_list, dtype=object))
        snp_conflicts = np.array(snp_conflicts, dtype=object)

    # Obtain a 2D array containing the hairpin information of all created padlock probes.
    # Each row consist of an array of integers; a 0 when no hairpin is formed, a 1 when a hairpin is formed.
    with open(hairpins, 'r') as handle:
        reader = csv.reader(handle)
        hairpin_scores = []
        for i, row in enumerate(reader):
            hairpin_score_list = []
            for probe in row:
                hairpin_score_list.append(int(probe))
            hairpin_scores.append(np.array(hairpin_score_list, dtype=object))
        hairpin_scores = np.array(hairpin_scores, dtype=object)

    # Obtain the 2D array of possible arm combinations for each target CpG in the Fasta file.
    # Each probe is a list containing the following information:
        # 0 [Sequence Upstream arm
        # 1 Sequence Downstream arm
        # 2 genomic location upstream arm
        # 3 genomic location downstream arm
        # 4 nr of probe for this target CpG
        # 5 target length
        # 6 standedness
        # 7 CpG Id]
    with open(probe_arms_file, 'rb') as handle:
        possible_arm_combinations_all_targets = loadall(handle)
        possible_arm_combinations_all_targets = np.array([np.array(possible_arm_combinations_one_target,dtype=object) for possible_arm_combinations_one_target in possible_arm_combinations_all_targets],dtype=object)
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
                probe_arm_list[i].append(probe)                                 # All info
                arm_upstream_loc_list[i].append(probe[2])                       # Genomic location of upstream arm sequence.
                arm_downstream_loc_list[i].append(probe[3])                     # Genomic location of downstream arm sequence.
                probe_id_list[i].append(str(probe[7]) + ':' + str(probe[4]))    # CpG identifier + number of created probe
                arm_upstream_list[i].append(probe[0])                           # Upstream arm sequence
                arm_downstream_list[i].append(probe[1])                         # Downstream arm sequence

    return probe_cpg_id_list, possible_probes_all_targets, tms, cpg_conflicts, snp_conflicts, hairpin_scores, \
        probe_arm_list, arm_upstream_loc_list, arm_downstream_loc_list, probe_id_list, arm_upstream_list, \
        arm_downstream_list, possible_arm_combinations_all_targets, tms_up, tms_down


def get_probe_costs_array(tms, cpg_conflicts, snp_conflicts, hairpin_array, scoring_weights):
    # Calculate probe cost from delta tm, hairpin formation and CpG and SNP conflicts.
    # Cost = (hairpin * 10^10)+ CpG conflicts + SNP conflicts + delta_tm
    weight_hairpin = int(scoring_weights['hairpin'])
    weight_cpg = int(scoring_weights['cpg'])
    weight_snp = int(scoring_weights['snp'])
    weight_tm = int(scoring_weights['tm'])
    probe_cost = np.add(np.add(np.add(np.add(np.multiply(hairpin_array, weight_hairpin),np.multiply(cpg_conflicts,weight_cpg)),np.multiply(snp_conflicts,weight_snp)), np.multiply(tms,weight_tm)),1)
    # Score is at least 1
    return probe_cost


def choose_probes_from_costs(probe_costs, possible_arm_combinations, n, counter, probe_id_list, seed):
    # Choose probes based on the current cost.
    # Cost changes each iteration
    np.random.seed(seed)  # Set seed from configuration file
    if probe_costs.size == 0:
        return
    else:
        probe_scores = 1 / probe_costs.astype(np.float64)
        sum_score = sum(probe_scores)
        if sum_score == 0:  # If all scores are 0, there should still be a probability to choose one of the probes
            probability_distribution = [1 / len(probe_scores)] * len(probe_scores)  # Create list of probabilities to choose probe
        else:
            probability_distribution = np.divide(probe_scores, sum_score)
        probe = choice(possible_arm_combinations, n, p=probability_distribution)  # Choose probe based on probability list
        index = possible_arm_combinations.index(probe[counter])
    return [probe[counter], probe_id_list[index]]


def get_conflict_ranges(conflicting_cpg_list_dimer,conflict_range_dimer, exclusion_factor):
    probe_list = []
    probe_list_count = []

    # Count how often a probe forms a dimer in the panel.
    for probe in conflicting_cpg_list_dimer:
        if probe in probe_list:
            probe_list_count[probe_list.index(probe)] += 1
        else:
            probe_list.append(probe)
            probe_list_count.append(1)
    if probe_list_count == []: #If there are no conflicts add a 0.
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
    # Combine conflict ranges from the S1 and S2 probes from mfeprimer dimer
    combined_conflicting_cpg_list_dimer = conflicting_cpg_list_dimer + conflicting_cpg_list_dimer2
    combined_conflict_range_dimer = conflict_range_dimer + conflict_range_dimer2

    # Keep top percentage of combined conflict dimers, based on exclusion factor.
    combined_most_conflicting_dimer_range_list = get_conflict_ranges(combined_conflicting_cpg_list_dimer,
                                                                     combined_conflict_range_dimer, exclusion_factor)
    return combined_most_conflicting_dimer_range_list


def get_dimer_range_upstream(possible_arm_combinations_all_targets, arm_upstream_loc_list, probe_index, nested_index,
                             dimer_start_num, dimer_end_num, dimer):
    # Obtain dimer range and template strand from the dimer forming region in the upstream arm of one precific probe
    # This translates a conflicting probe into a conflicting region to form a probe arm from
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
    # Obtain dimer range and template strand from the dimer forming region in the downstream arm of one precific probe
    # This translates a conflicting probe into a conflicting region to form a probe arm from
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


def create_conflicting_indices_list_bedtools(dimer_range_list, dimer_bedfile_name, bedfile_intersect_name, outputdir):
    write_loci_to_bedfile(dimer_bedfile_name, dimer_range_list)
    # Do bedtools intersect here on combined_loci_bedfile_name and dimer_bedfile_name
    os.system('bedtools intersect -wa -s -a ' + outputdir + 'possible_probe_info/combined.bed -b ' + dimer_bedfile_name + ' -f 7E-9 > ' 
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
                                                  exclusion_factor, outputdir, score_cutoff, tm_cutoff, nr_of_cores):
    # Create a fasta file with the chosen probes in this iteration.
    passed_list = []
    with open('tmp_fasta_chosen_probes.fasta', 'w') as handle:
        for i, probe in enumerate(chosen_probes):
            if probe is None:
                passed_list.append(i)
            else:
                handle.write('>' + probe[1] + '\n')
                handle.write(probe[0] + '\n')

    # Test the chosen probes set for dimers
    os.system(
        'mfeprimer dimer -i tmp_fasta_chosen_probes.fasta -j -o tmp_dimers_chosen_probes -s '+str(score_cutoff)+' -t '+str(tm_cutoff) + ' -c ' + str(nr_of_cores))
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

    # Obtain dimer ranges from dimer list of mfeprimer dimer output.
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
            probe_index = probe_cpg_id_list.index(cpg_id.split(':')[0])
            nested_index = probe_id_list[probe_index].index(cpg_id)
            probe_index2 = probe_cpg_id_list.index(cpg_id2.split(':')[0])
            nested_index2 = probe_id_list[probe_index2].index(cpg_id2)
            # obtain the genomic locations of the dimer forming regions for S1
            # check whether the S1 sequence is 5' to 3' or the other way around.
            if arm_upstream_list[probe_index][nested_index].strip("'") in dimer['S1']['Seq'][:len(
                    arm_upstream_list[probe_index][nested_index].strip("'"))]:
                # Check whether the dimer in in the upstream or downstream arm of the S1 probe
                if dimer_start_num < len(
                        arm_upstream_list[probe_index][nested_index].strip("'")):  # dimer is in upstream arm
                    dimer_range = get_dimer_range_upstream(possible_arm_combinations_all_targets, arm_upstream_loc_list,
                                                           probe_index, nested_index, dimer_start_num, dimer_end_num,
                                                           dimer)
                    conflict_range_dimer.append(dimer_range)
                    conflicting_cpg_list_dimer.append(cpg_id)
                elif dimer_end_num >= s1_length - len(arm_downstream_list[probe_index][nested_index].strip(
                        "'")):  # dimer is in downstream arm.
                    # If dimer end has no overlap with the downsteam or upstream arm is should be in the backbone.
                    dimer_range = get_dimer_range_downstream(possible_arm_combinations_all_targets,
                                                             arm_downstream_loc_list, probe_index, nested_index,
                                                             dimer_start_num, dimer_end_num, s1_length, dimer)
                    conflict_range_dimer.append(dimer_range)
                    conflicting_cpg_list_dimer.append(cpg_id)
                else:  # dimer forming region is in backbone,
                    # only report the other probe that forms the dimer together with this probe.
                    pass
            else:  # In this situation the S1 sequence is 3'-5', therefore we need to correct the indices.
                dimer_start_num, dimer_end_num = s1_length - dimer_end_num, s1_length - dimer_start_num
                # Check whether the dimer in in the upstream or downstream arm of the S1 probe
                if dimer_start_num < len(
                        arm_upstream_list[probe_index][nested_index].strip("'")):  # dimer is in upstream arm
                    dimer_range = get_dimer_range_upstream(possible_arm_combinations_all_targets, arm_upstream_loc_list,
                                                           probe_index, nested_index, dimer_start_num, dimer_end_num,
                                                           dimer)
                    conflict_range_dimer.append(dimer_range)
                    conflicting_cpg_list_dimer.append(cpg_id)
                elif dimer_end_num >= s1_length - len(arm_downstream_list[probe_index][nested_index].strip(
                        "'")):  # dimer is in downstream arm
                    # If dimer end has no overlap with the downsteam or upstream arm is should be in the backbone.
                    dimer_range = get_dimer_range_downstream(possible_arm_combinations_all_targets,
                                                             arm_downstream_loc_list, probe_index, nested_index,
                                                             dimer_start_num, dimer_end_num, s1_length, dimer)
                    conflict_range_dimer.append(dimer_range)
                    conflicting_cpg_list_dimer.append(cpg_id)
                else:  # dimer forming region is in backbone,
                    # only report the other probe that forms the dimer together with this probe.
                    pass

            # obtain the genomic locations of the dimer forming regions for S2
            # check whether the S2 sequence is 5' to 3' or the other way around.
            if arm_upstream_list[probe_index2][nested_index2].strip("'") in dimer['S2']['Seq'][:len(
                    arm_upstream_list[probe_index2][nested_index2].strip("'"))]:
                # Check whether the dimer in in the upstream or downstream arm of the S2 probe
                if dimer_start_num2 < len(
                        arm_upstream_list[probe_index2][nested_index2].strip("'")):  # dimer is in upstream arm
                    dimer_range2 = get_dimer_range_upstream(possible_arm_combinations_all_targets,
                                                            arm_upstream_loc_list, probe_index2, nested_index2,
                                                            dimer_start_num2, dimer_end_num2, dimer)
                    conflict_range_dimer2.append(dimer_range2)
                    conflicting_cpg_list_dimer2.append(cpg_id2)
                elif dimer_end_num2 >= s2_length - len(arm_downstream_list[probe_index2][nested_index2].strip(
                        "'")):  # dimer is in downstream arm.
                    # If dimer end has no overlap with the downsteam or upstream arm is should be in the backbone.
                    dimer_range2 = get_dimer_range_downstream(possible_arm_combinations_all_targets,
                                                              arm_downstream_loc_list, probe_index2, nested_index2,
                                                              dimer_start_num2, dimer_end_num2, s2_length, dimer)
                    conflict_range_dimer2.append(dimer_range2)
                    conflicting_cpg_list_dimer2.append(cpg_id2)
                else:  # dimer forming region is in backbone,
                    # only report the other probe that forms the dimer together with this probe.
                    pass
            else:  # In this situation the S2 sequence is 3'-5', therefore we need to correct the indices.
                dimer_start_num2, dimer_end_num2 = s2_length - dimer_end_num2, s2_length - dimer_start_num2
                # Check whether the dimer in in the upstream or downstream arm of the S2 probe
                if dimer_start_num2 < len(
                        arm_upstream_list[probe_index2][nested_index2].strip("'")):  # dimer is in upstream arm
                    dimer_range2 = get_dimer_range_upstream(possible_arm_combinations_all_targets,
                                                            arm_upstream_loc_list, probe_index2, nested_index2,
                                                            dimer_start_num2, dimer_end_num2, dimer)
                    conflict_range_dimer2.append(dimer_range2)
                    conflicting_cpg_list_dimer2.append(cpg_id2)
                elif dimer_end_num2 >= s2_length - len(arm_downstream_list[probe_index2][nested_index2].strip(
                        "'")):  # dimer is in downstream arm.
                    # If dimer end has no overlap with the downsteam or upstream arm is should be in the backbone.
                    dimer_range2 = get_dimer_range_downstream(possible_arm_combinations_all_targets,
                                                              arm_downstream_loc_list, probe_index2, nested_index2,
                                                              dimer_start_num2, dimer_end_num2, s2_length, dimer)
                    conflict_range_dimer2.append(dimer_range2)
                    conflicting_cpg_list_dimer2.append(cpg_id2)
                else:  # dimer forming region is in backbone,
                    # only report the other probe that forms the dimer together with this probe.
                    pass

    # Only report the probes that formed the most amount of dimers.
    # Exclusion factor determines the amount of conflicting probes that is selected. 1 = all, 0.1 = top 10%, etc.
    most_conflicting_dimer_range_list = get_conflict_ranges_s1_and_s2(conflicting_cpg_list_dimer, conflict_range_dimer,
                                                                      conflicting_cpg_list_dimer2,
                                                                      conflict_range_dimer2,
                                                                      exclusion_factor)

    # Find all probes that have common regions with the most_conflicting_dimer_range and report as conflicting probes
    all_conflicting_probe_list = []
    conflicting_targets = set([])
    print('\t\t'+'Dimer ranges obtained')
    dimer_bedfile_name = outputdir + 'conflict_info/conflicting_dimer.bed'
    bedfile_intersect_name = outputdir + 'conflict_info/conflicting_loci.bed'
    new_conflicting_indices_list = create_conflicting_indices_list_bedtools(most_conflicting_dimer_range_list,
                                                                            dimer_bedfile_name, bedfile_intersect_name,
                                                                            outputdir)
    # Increase the cost of probes are conflicting and are more likely to form dimers.
    for indexes in new_conflicting_indices_list:
        i = int(indexes[0])
        j = int(indexes[1])
        all_conflicting_probe_list.append(probe_id_list[i][j])
        conflicting_targets.add(probe_id_list[i][0].split(':')[0])
        probe_costs[i][j] = probe_costs[i][j] + float(10 ** 10)  # Increase cost of conflicting probe
    return all_conflicting_probe_list, conflicting_targets, probe_costs, cpgs_of_dimer_forming_probes

def increase_cost_probes_with_out_of_bounds_tm(chosen_probes, probe_id_list, probe_costs, tms_down, tms_up, scoring_weights, max_delta_tm_panel):
    # Get indices of chosen_probes
    indices=[]
    for i,probe in enumerate(chosen_probes):
        if probe == None:
            pass
        else:
            probe_id = probe[1]
            index = probe_id_list[i].index(probe_id)
            indices.append([i,index])
    # Find median Tm of chosen panel
    tm_up_list = [tms_up[index[0]][index[1]] for index in indices]
    tm_down_list = [tms_down[index[0]][index[1]] for index in indices]
    tm_list = tm_up_list + tm_down_list
    median_tm = statistics.median(tm_list)
    upperbound_tm = median_tm + 0.5 * max_delta_tm_panel  # Set upper Tm limit.
    lowerbound_tm = median_tm - 0.5 * max_delta_tm_panel  # Set lower Tm limit.
    weight_tm = int(scoring_weights['tm'])
    
    # Find amount of chosen probes that have a Tm outside of the Tm range around the median Tm.
    tm_amount_out_of_range = sum([1 for tm in tm_list if tm > upperbound_tm]) + sum([1 for tm in tm_list if tm  < lowerbound_tm])

    # Find all probes that have a Tm outside of the Tm range around the median Tm.
    # Add the following amount to the score: scoring weight Tm * difference between probe Tm and median Tm.
    increase_1 = np.array([np.where(row < upperbound_tm, row , np.multiply(np.abs(np.subtract(row,median_tm)),weight_tm)) for row in tms_down],dtype = object)
    increase_2 = np.array([np.where(row > lowerbound_tm, row , np.multiply(np.abs(np.subtract(row,median_tm)),weight_tm)) for row in tms_down],dtype = object)
    increase_3 = np.array([np.where(row < upperbound_tm, row , np.multiply(np.abs(np.subtract(row,median_tm)),weight_tm)) for row in tms_up],dtype = object)
    increase_4 = np.array([np.where(row > lowerbound_tm, row , np.multiply(np.abs(np.subtract(row,median_tm)),weight_tm)) for row in tms_up],dtype = object)
    total_increase = np.add(np.add(increase_1,increase_2),np.add(increase_3,increase_4))
    new_probe_costs = np.add(probe_costs, total_increase)
    return new_probe_costs, tm_amount_out_of_range

def get_dimer_scores(chosen_set, score_cutoff, tm_cutoff, targets, nr_of_cores):
    # Create a fasta file with the chosen probes in this iteration.
    passed_list = []
    with open('tmp_fasta_chosen_probes.fasta', 'w') as handle:
        for i, probe in enumerate(chosen_set):
            if probe is None:
                passed_list.append(i)
            else:
                handle.write('>' + probe[1] + '\n')
                handle.write(probe[0] + '\n')

    # Test the chosen probes set for dimers
    os.system(
        'mfeprimer dimer -i tmp_fasta_chosen_probes.fasta -j -o tmp_dimers_chosen_probes -s '+str(score_cutoff)+' -t '+str(tm_cutoff) + ' -c ' + str(nr_of_cores))
    # Score cut-off and temperature cut-off are set in configuration file

    # Obtain dimer list from mfeprimer output
    os.system('rm tmp_fasta_chosen_probes.fasta tmp_dimers_chosen_probes')
    dimer_file = 'tmp_dimers_chosen_probes.json'
    with open(dimer_file) as jsonFile:
        dimerlist = json.load(jsonFile)
        jsonFile.close()

    dimer_scores = [0]*len(chosen_set)
    genomic_loci_list=[]
    with open(targets,'r') as handle:
        reader = csv.reader(handle, delimiter = '\t')
        for row in reader:
            genomic_loci=row[-1]
            genomic_loci_list.append(genomic_loci)


    # Obtain dimer scores from dimer list of mfeprimer dimer output.
    if dimerlist is None:
        pass
    else:
        for dindex, dimer in enumerate(dimerlist):
            cpg_loc1 = dimer['S1']['ID'].split(':')[0]
            cpg_loc2 = dimer['S2']['ID'].split(':')[0]
            tm = round(float(dimer['Tm']),2)
            for genomic_loci in [cpg_loc1,cpg_loc2]:
                index = genomic_loci_list.index(genomic_loci)
                dimer_scores[index] = max(dimer_scores[index],tm)
    return dimer_scores


def write_output(targets, output_name, chosen_set, conflicting_file, conflicting_probe_list, min_dimers, panel_output_file, probe_arm_list, backbone_sequence, tms, cpg_conflicts, snp_conflicts, hairpin_scores, probe_id_list, dimer_scores, tms_up, tms_down):
    # Write the output of the Choose_probes rule.
        # 1. Fasta file with chosen panel
        # 2. Tab delimiter file with conflicting probes per iteration

    # Create discriptions for chosen probes, based on the targets.
    seq_list = []
    with open(targets, 'r') as handle:
        reader = csv.reader(handle, delimiter='\t')
        description_list = []
        for row in reader:
            loc=int((int(row[1])+int(row[2])-1)/2)  #Obtain middle value from target bed file. Middle value is the location of the CpG.
            description_list.append(str(row[0] + ':' + str(loc)  + '\t' + row[-1]))  # Genomic locus \t nr of chosen probe for CpG

    # Choose panel with the least amount of dimers.
    # Write probes of chosen panel to a fasta file.
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

    # Return conflicting probes per iteration to a file.
    with open(conflicting_file, 'w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        for row in conflicting_probe_list:
            writer.writerow(row)

    
    # Return tab separated file with information about probe panel
    with open(panel_output_file, 'w') as handle:
        writer = csv.writer(handle, delimiter = ',')
        header = ['Genomic_locus', 'CpG_id', 'Probe_sequence', 'Upstream_arm_sequence', 'Backbone_sequence', 'Downstream_arm_sequence', 
                    'Genomic_locus_upstream_arm', 'Genomic_locus_downstream_arm', 'Target_length', 'Target_strand', 'Delta_Tm',
                    'Tm_upstream_arm', 'Tm_downstream_arm', 'CpGs_in_arms', 'SNPs_in_arms', 'Hairpin_score', 'Dimer_score', 'upstream_cg_percentage', 'downstream_cg_percentage']
        writer.writerow(header)
        for i,probe in enumerate(chosen_set):
            if probe == None:
                locus = description_list[i].split('\t')[0]  # 0 Genomic locus
                cpg_id = description_list[i].split('\t')[1]  # 1 CpG id
                print('It was not possible to design a probe with the current parameters for target nr '+str(cpg_id))
                row = [locus, cpg_id]
                writer.writerow(row)
            else:
                probe_id=int(probe[1].split(':')[-1])
                locus = description_list[i].split('\t')[0]  # 0 Genomic locus
                cpg_id = description_list[i].split('\t')[1]  # 1 CpG id
                probe_sequence = probe[0]  # 2 sequence
                upstream_sequence = probe_arm_list[i][probe_id][0]  # 3 Upstream arm sequence
                backbone_sequence = backbone_sequence  # 4 Backbone sequence
                downstream_sequence = probe_arm_list[i][probe_id][1]  # 5 Downstream arm sequence
                upstream_locus =  probe_arm_list[i][probe_id][2]  # 6 Genomic location upstream arm
                downstream_locus = probe_arm_list[i][probe_id][3]  # 7 Genomic location downstream arm
                target_length = probe_arm_list[i][probe_id][5]  # 8 Target length
                strandedness = probe_arm_list[i][probe_id][6]  # 9 Which strand is the target strand
                delta_tm = tms[i][probe_id]  # 10 Difference in Tm between probe arms
                tm_up = tms_up[i][probe_id]  # 11 Tm upstream arm
                tm_down = tms_down[i][probe_id]  # 12 Tm downstream arm
                cpg_conflict = cpg_conflicts[i][probe_id]  # 13 Amount of CpGs in arms
                snp_conflict = snp_conflicts[i][probe_id]  # 14 Amount of frequent SNPs in arms
                hairpin_score = hairpin_scores[i][probe_id]  # 15 Hairpin score of probe
                dimer_score = dimer_scores[i]  # 16 Highest dimer score of probe
                upstream_cg_percentage = GC(upstream_sequence)
                downstream_cg_percentage = GC(downstream_sequence)
                row=[locus, cpg_id, probe_sequence, upstream_sequence, backbone_sequence, downstream_sequence, upstream_locus,
                        downstream_locus, target_length, strandedness, delta_tm, tm_up, tm_down, cpg_conflict, snp_conflict,
                        hairpin_score, dimer_score, upstream_cg_percentage, downstream_cg_percentage]
                writer.writerow(row)

if __name__ == '__main__':
    main()
