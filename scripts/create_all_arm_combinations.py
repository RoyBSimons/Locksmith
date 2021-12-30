#!/usr/bin/python3
# Create all possible Padlock probes
import os
from argparse import ArgumentParser
import json
from Bio import SeqIO, Seq, SeqRecord
from Bio.SeqUtils import GC
import multiprocessing as mp
import csv
import numpy as np
import pickle
import time


def main():
    parser = get_arg_parser()  # Parse input to variables
    args = vars(parser.parse_args())
    filename = args["filename"]
    bed_file = args["bed_file"]
    output_dir = args["output_dir"] + '/'
    config_file = args["config_file"]
    nr_of_cores = int(args["cores"])
    acc_nr_to_chrom_nr_file = args["acc_nr_to_chrom_nr_file"]

    ftp_path_snp_database, probe_specifics, min_arm_length, max_arm_length, min_target_length, max_target_length, \
        target_range, cpg_flanks, max_cg_percentage, min_cg_percentage, backbone_sequence, mid_loc, freq_threshold, \
        score_cutoff, tm_cutoff = import_config(config_file)  # Import parameters set in configuration file

    start = time.time()
    cpg_id_list = create_cpg_id_list(bed_file)  # Obtain CpG list from target bed_file
    with open(filename) as handle:  # Create a list of possible arm combinations for each target CpG in the Fasta file.
        with mp.Pool(nr_of_cores) as pool_handle:  # Use the maximum amount of cores from the configuration file
            possible_arm_combinations_all_targets = [pool_handle.apply(create_all_possible_arms_both_strands, 
                                                                       args=(record, min_target_length, 
                                                                             max_target_length, min_arm_length, 
                                                                             max_arm_length, min_target_length, 
                                                                             max_cg_percentage, cpg_flanks, mid_loc, 
                                                                             cpg_id_list[i])) for i, record in 
                                                     enumerate(SeqIO.parse(handle, "fasta"))]  
            # Create all possible arms taking into account the parameters set in the configuration file
    end = time.time()
    print(round(end - start, 2))  # print elapsed time to log file
    print('All possible arms obtained')  # print progress to log file

    start = time.time()
    tms, tms_up, tms_down = get_delta_tm_array(possible_arm_combinations_all_targets)  
    # Obtain a 2D array containing the difference in Tm between the upstream and downstream arm of the padlock probe. 
    # Each row consist of the delta-tms for all possible padlock probes to target one CpG.
    end = time.time()
    print(round(end - start, 2))  # print elapsed time to log file
    print('\ttms obtained')  # print progress to log file

    start = time.time()
    cpg_conflicts = report_cpgs_in_arms(possible_arm_combinations_all_targets)  
    # Obtain a 2D array containing the amount of CpGs in the arms of the padlock probe. 
    # Each row consist of the amount of CpGs in the arms for all possible padlock probes to target one CpG.
    end = time.time()
    print(round(end - start, 2))  # print elapsed time to log file
    print('\tCpGs in arms counted')  # print progress to log file

    start = time.time()
    possible_probes_all_targets = add_backbone_array(possible_arm_combinations_all_targets, backbone_sequence)  
    # Obtain a 2D array containing the padlock probe. Each row consist of an array of strings containing the full 
    # padlock probe sequences by combining the arms together with the backbone sequence from the configuration file.
    end = time.time()
    print(round(end - start, 2))  # print elapsed time to log file
    print('\tbackbone added')  # print progress to log file

    start = time.time()
    fasta_name = 'temp_test_fasta'
    output_name_json = 'temp_test_json'
    with mp.Pool(nr_of_cores) as pool_handle:
        hairpin_scores = [pool_handle.apply(check_probe_for_hairpin_score, 
                                            args=(possible_probes, fasta_name, output_name_json, index, score_cutoff, tm_cutoff))
                          for index, possible_probes in enumerate(possible_probes_all_targets)]  
        # Obtain a 2D array containing the hairpin information of all created padlock probes. 
        # Each row consist of an array of integers; a 0 when no hairpin is formed, a 1 when a hairpin is formed.
    end = time.time()
    print(round(end - start, 2))  # print elapsed time to log file
    print('\tHairpins detected')  # print progress to log file

    start = time.time()
    snp_conflicts = obtain_snps(possible_arm_combinations_all_targets, bed_file, freq_threshold, ftp_path_snp_database,
                                acc_nr_to_chrom_nr_file) 
    # Obtain a 2D array containing the amount of frequent SNPs in the arms of the padlock probe. 
    # Each row consist of the amount of SNPs in the arms for all possible padlock probes to target one CpG. 
    # Frequency threshold is set in configuration file
    end = time.time()
    print(round(end - start, 2))  # print elapsed time to log file
    print('\tSNPs in arms counted')  # print progress to log file

    start = time.time()
    create_output_files(output_dir, tms, cpg_conflicts, possible_probes_all_targets, hairpin_scores, snp_conflicts,
                        possible_arm_combinations_all_targets, tms_up, tms_down)  
    # Write all created padlock probes and their parameters to files.
    end = time.time()
    print(round(end - start, 2))  # print elapsed time to log file
    print('\tPadlock probes and their parameters are written to their output files')  # print progress to log file


def get_arg_parser():
    # parse all files
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", dest="filename",
                        help="open FILE", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output_dir",
                        help="Directory for output", metavar="output_dir")
    parser.add_argument("-c", "--config_file", dest="config_file",
                        help="Config file with probe specifics", metavar="JSONfile")
    parser.add_argument("-b", "--bed", dest="bed_file",
                        help="Bed file containing the targets", metavar="BED")
    parser.add_argument("-t", "--cores", dest="cores",
                        help="Passed amount of cores to script", metavar="Cores")
    parser.add_argument("-a", "--acc_2_chrom", dest="acc_nr_to_chrom_nr_file",
                        help="Tab delimiter file containing the chromosome nr and the accompanied accesion nr "
                             "used in the ftp database used for SNP detection", metavar="ACC2CHROM")
    return parser


def import_config(config_file):  # import config file
    with open(config_file) as jsonFile:
        config_object = json.load(jsonFile)
        jsonFile.close()
    ftp_path_snp_database = config_object['ftp_path_snp_database']
    probe_specifics = config_object['probe_specifics'][0]
    min_arm_length = probe_specifics['min_arm_length']
    max_arm_length = probe_specifics['max_arm_length']
    min_target_length = probe_specifics['min_target_length']
    max_target_length = probe_specifics['max_target_length']
    target_range = config_object['target_range']
    cpg_flanks = probe_specifics['cpg_flanks']
    max_cg_percentage = float(probe_specifics['max_cg_percentage'])
    min_cg_percentage = float(probe_specifics['min_cg_percentage'])

    # patch backbone together
    backbone = config_object['backbone_sequence'][0]
    rev_compl_universal_forward_primer = backbone['reverse_complement_universal_forward_primer']
    common_region = backbone['common_region']
    universal_reverse_primer = backbone['universal_reverse_primer']
    backbone_sequence = rev_compl_universal_forward_primer + common_region + universal_reverse_primer

    mid_loc = int(target_range) + 1  
    # The middle nucleotide is the target range +1: the total target is created by adding the target range on both ends.
    freq_threshold = float(config_object['snp_frequency_threshold'])

    score_cutoff = int(config_object['mfeprimer_hairpin_parameters'][0]['score_cutoff'])
    tm_cutoff = int(config_object['mfeprimer_hairpin_parameters'][0]['tm_cutoff'])

    return ftp_path_snp_database, probe_specifics, min_arm_length, max_arm_length, min_target_length, \
        max_target_length, target_range, cpg_flanks, max_cg_percentage, min_cg_percentage, backbone_sequence, mid_loc, \
        freq_threshold, score_cutoff, tm_cutoff


def create_cpg_id_list(bed_file):
    cpg_id_list = []  # initialize list of target CpGs
    with open(bed_file, 'r') as handle:
        reader = csv.reader(handle, delimiter='\t')
        for row in reader:
            cpg_id_list.append(row[3])
    return cpg_id_list


def create_all_possible_arms_both_strands(record, min_target_length, max_target_length, min_arm_length, max_arm_length,
                                          min_cg_percentage, max_cg_percentage, cpg_flanks, mid_loc, cpg_id):
    # Obtain the list of possible arm combinations for one record (target CpG in the Fasta file).
    # Each nested value is a list containing the following information:
        # 0 [Sequence Upstream arm
        # 1 Sequence Downstream arm
        # 2 genomic location upstream arm
        # 3 genomic location downstream arm
        # 4 nr of probe for this target CpG
        # 5 target length
        # 6 standedness
        # 7 CpG Id]

    probe_list = []
    i = 0
    rev_record = record.reverse_complement()
    record_length = len(record.seq)
    for target_length in range(min_target_length,
                               max_target_length + 1):  # loop over range of target lengths, including the maximum
        for arm_length_downstream in range(min_arm_length, max_arm_length + 1):
            start_loc_downstream = mid_loc - (target_length - cpg_flanks)
            end_loc_downstream = mid_loc - cpg_flanks
            for start_loc in range(start_loc_downstream, end_loc_downstream):
                start_loc_upstream = start_loc + target_length
                new_d_arm = record.seq[start_loc - arm_length_downstream:start_loc]
                if GC(new_d_arm) > max_cg_percentage or GC(new_d_arm) < min_cg_percentage:
                    break
                else:
                    pass
                downstream_id = record.id.split(':')[0] + ":" + str(
                    int(record.id.split(':')[1].split("-")[0]) + start_loc - arm_length_downstream) + "-" + str(
                    int(record.id.split(':')[1].split("-")[0]) + start_loc - 1)
                for arm_length_upstream in range(min_arm_length, max_arm_length + 1):
                    new_u_arm = record.seq[start_loc_upstream:start_loc_upstream + arm_length_upstream]
                    if GC(new_u_arm) > max_cg_percentage or GC(new_d_arm) < min_cg_percentage:
                        break
                    else:
                        pass
                    upstream_id = record.id.split(':')[0] + ":" + str(
                        int(record.id.split(':')[1].split("-")[0]) + start_loc_upstream) + "-" + str(
                        int(record.id.split(':')[1].split("-")[0]) + start_loc_upstream + arm_length_upstream - 1)
                    probe = [str(new_u_arm), str(new_d_arm), upstream_id, downstream_id, i, target_length, '+', cpg_id]
                    i += 1
                    probe_list.append(probe)
    for target_length in range(min_target_length,
                               max_target_length + 1):  # loop over range of target lengths, including the maximum
        for arm_length_downstream in range(min_arm_length, max_arm_length + 1):
            start_loc_downstream = mid_loc - (target_length - cpg_flanks)
            end_loc_downstream = mid_loc - cpg_flanks
            for start_loc in range(start_loc_downstream, end_loc_downstream):
                start_loc_upstream = start_loc + target_length
                new_d_arm_rev = rev_record.seq[start_loc - arm_length_downstream:start_loc]
                if GC(new_d_arm_rev) > max_cg_percentage or GC(new_d_arm_rev) < min_cg_percentage:
                    break
                else:
                    pass
                for arm_length_upstream in range(min_arm_length, max_arm_length + 1):
                    new_u_arm_rev = rev_record.seq[start_loc_upstream:start_loc_upstream + arm_length_upstream]
                    if GC(new_u_arm_rev) > max_cg_percentage or GC(new_u_arm_rev) < min_cg_percentage:
                        break
                    else:
                        pass
                    upstream_id_rev = record.id.split(':')[0] + ":" + str(int(record.id.split(':')[1].split("-")[0]) + (
                                record_length - start_loc_upstream - arm_length_upstream)) + "-" + str(
                        int(record.id.split(':')[1].split("-")[0]) + (record_length - start_loc_upstream - 1))
                    downstream_id_rev = record.id.split(':')[0] + ":" + str(
                        int(record.id.split(':')[1].split("-")[0]) + (record_length - start_loc)) + "-" + str(
                        int(record.id.split(':')[1].split("-")[0]) + (
                                    record_length - start_loc + arm_length_downstream - 1))
                    probe = [str(new_u_arm_rev), str(new_d_arm_rev), upstream_id_rev, downstream_id_rev, i,
                             target_length, '-', cpg_id]
                    i += 1
                    probe_list.append(probe)
    return probe_list


def get_delta_tm_array(probe_arms_array):
    # Calculate the delta Tm between the arms of each probe.
    A = [[string[0].count('A') for string in row] for row in probe_arms_array]
    T = [[string[0].count('T') for string in row] for row in probe_arms_array]
    G = [[string[0].count('G') for string in row] for row in probe_arms_array]
    C = [[string[0].count('C') for string in row] for row in probe_arms_array]
    tm_up = np.array([np.array(np.round(np.add(64.9, np.divide(np.multiply(41, np.add(G[i], np.subtract(C[i], 16.4))),
                                    np.add(np.add(A[i], T[i]), np.add(C[i], G[i])))),1)) for i in
             range(len(A))])  # Tm=64.9+41*(G+C-16.4)/(A+T+C+G)
    A = [[string[1].count('A') for string in row] for row in probe_arms_array]
    T = [[string[1].count('T') for string in row] for row in probe_arms_array]
    G = [[string[1].count('G') for string in row] for row in probe_arms_array]
    C = [[string[1].count('C') for string in row] for row in probe_arms_array]
    tm_down = np.array([np.array(np.round(np.add(64.9, np.divide(np.multiply(41, np.add(G[i], np.subtract(C[i], 16.4))),
                                      np.add(np.add(A[i], T[i]), np.add(C[i], G[i])))),1)) for i in
               range(len(A))])  # Tm=64.9+41*(G+C-16.4)/(A+T+C+G)
    tms = np.array([np.array([np.round(np.absolute(np.subtract(tm_down[i][j], tm_up[i][j])), 1) for j in range(len(A[i]))]) for i in
           range(len(A))])
    return tms, tm_up, tm_down


def report_cpgs_in_arms(probe_arms_array):
    # Count the amount of CpGs in the arms of each probe.
    upstream_count = [[string[0].count('CG') for string in row] for row in probe_arms_array]
    downstream_count = [[string[1].count('CG') for string in row] for row in probe_arms_array]
    counts_array = [np.add(upstream_count[i], downstream_count[i]) for i in range(len(upstream_count))]
    return counts_array


def add_backbone_array(probe_arms_array, backbone_sequence):
    # Construct probe sequence (5'-3') from arms and backbone
    probe_array = [[string[0] + backbone_sequence + string[1] for string in row] for row in probe_arms_array]
    return probe_array


def check_probe_for_hairpin_score(probes, fasta_name, output_name_json, index, score_cutoff, tm_cutoff):
    # Report score of hairpin when a probe forms one; else report a zero.
    seq_list = []
    list_len = len(probes)
    if list_len == 0:  # in case there is no probes
        return []
    with open(fasta_name + str(index) + '.fa', 'w') as output_handle:
        for i, probe in enumerate(probes):
            probe_name = str(i)
            seq_list.append(SeqRecord.SeqRecord(Seq.Seq(probe), id=probe_name))
        SeqIO.write(seq_list, output_handle, 'fasta')
    os.system('mfeprimer hairpin --in ' + fasta_name + str(index) + '.fa -j -o ' + output_name_json + str(
        index) + ' -s '+str(score_cutoff)+' -t '+str(tm_cutoff))  # Remove/move the output files;
    # Report hairpins which go over the set thresholds for score and Tm
    with open(output_name_json + str(index) + '.json') as jsonFile:
        hairpin_list = json.load(jsonFile)
        jsonFile.close()
    os.system('rm ' + output_name_json + str(index) + '.json')
    os.system('rm ' + output_name_json + str(index))
    os.system('rm ' + fasta_name + str(index) + '.fa')
    hairpin_scores = [0]*list_len
    if hairpin_list is None:
        pass
    else:
        for hairpin in hairpin_list:
            hairpin_score = hairpin['Score']
            probe_nr = int(hairpin['Seq']['ID'])
            hairpin_scores[probe_nr] = hairpin_score
    return hairpin_scores


def obtain_acc_nr_and_chrom_nr_lists(acc_nr_to_chrom_nr_file):
    # This conversion is needed for couting the frequent SNPs in the probe arms.
    chrom_nr_list = []
    acc_nr_list = []
    with open(acc_nr_to_chrom_nr_file, 'r') as handle:
        reader = csv.reader(handle, delimiter='\t')
        next(reader)  # Read the header
        for row in reader:
            chrom_nr_list.append(row[0])
            acc_nr_list.append(row[1])
    return chrom_nr_list, acc_nr_list


def locus_with_accession_nr_bed_info_with_chromosome_nr(locus, chrom_nr_list, acc_nr_list):  
    # Tab seperated file with header.
    # #Chromosome	Accession.version
    chrom = locus.split(':')[0][3:]
    left_index = locus.split(':')[1].split('-')[0]
    right_index = locus.split('-')[1]
    bed_info = [acc_nr_list[chrom_nr_list.index(chrom)], left_index, right_index]
    return bed_info


def rewrite_bed_file_from_accession_nr_to_chromosome_nr(bed_path, acc_nr_to_chrom_nr_file):
    # This conversion is needed for couting the frequent SNPs in the probe arms.
    chrom_nr_list, acc_nr_list = obtain_acc_nr_and_chrom_nr_lists(acc_nr_to_chrom_nr_file)
    bed_info_list = []
    with open(bed_path) as handle:
        reader = csv.reader(handle, delimiter='\t')
        for row in reader:
            chrom = row[0][3:]
            bed_info = [acc_nr_list[chrom_nr_list.index(chrom)], row[1], row[2], row[3]]
            bed_info_list.append(bed_info)
    new_path = bed_path + '_chrom_nrs'
    with open(new_path, 'w') as output_handle:
        writer = csv.writer(output_handle, delimiter='\t')
        writer.writerows(bed_info_list)
    return new_path


def obtain_snps(probe_list, bed_path, freq_threshold, ftp_path_snp_database, acc_nr_to_chrom_nr_file):
    # first obtain frequent SNPs in targets
    new_path = rewrite_bed_file_from_accession_nr_to_chromosome_nr(bed_path, acc_nr_to_chrom_nr_file)
    os.system('tabix ' + ftp_path_snp_database + ' -R ' + new_path + ' > tmp_output_snp')
    chr_list = []
    loc_list = []
    freq_threshold = 1.0 - freq_threshold
    with open('tmp_output_snp') as input_file:
        input_reader = csv.reader(input_file)
        for row in input_reader:
            if "FREQ" in row[0]:
                freq = float(row[0].split(";")[-1].split(":")[-1])
                if freq < freq_threshold:
                    chrom = row[0].split('\t')[0]
                    locus = int(row[0].split('\t')[1])
                    if chrom in chr_list:
                        loc_list[chr_list.index(chrom)].append(locus)
                    else:
                        chr_list.append(chrom)
                        loc_list.append([locus])
    os.system('rm tmp_output_snp')

    #Find SNPs in probe arms, and return found amount
    chrom_nr_list, acc_nr_list = obtain_acc_nr_and_chrom_nr_lists(acc_nr_to_chrom_nr_file)
    snp_count=[]
    for i, probe_arms in enumerate(probe_list):
        snp_count.append([0]*len(probe_arms))
        for j, probe in enumerate(probe_arms):
            arm1 = locus_with_accession_nr_bed_info_with_chromosome_nr(probe[2], chrom_nr_list, acc_nr_list)
            arm1.append(str(i) + '-' + str(j))
            arm2 = locus_with_accession_nr_bed_info_with_chromosome_nr(probe[3], chrom_nr_list, acc_nr_list)
            arm2.append(str(i) + '-' + str(j))
            if arm1[0] in chr_list:
                index_to_check = chr_list.index(arm1[0])
                for locus in loc_list[index_to_check]:
                    if locus in range(int(arm1[1]), int(arm1[2]) + 1):  # check whether the SNP is in the arm
                        snp_count[i][j] += 1
                    if locus in range(int(arm2[1]), int(arm2[2]) + 1):  # check whether the SNP is in the arm
                        snp_count[i][j] += 1
    return snp_count


# Store all information from the possible probes into files
def create_output_files(output_dir, tms, cpg_conflicts, possible_probes_all_targets, hairpin_scores, snp_conflicts,
                        possible_arm_combinations_all_targets, tms_up, tms_down):
    start = time.time()
    with open(output_dir + 'tms.csv', 'wb') as file:
        pickle.dump(tms, file)
    with open(output_dir + 'cpg_conflicts.csv', 'w') as file:
        for item in cpg_conflicts:
            file.write(",".join(map(str, item)))
            file.write("\n")
    with open(output_dir + 'possible_probes_all_targets.csv', 'w') as file:
        for item in possible_probes_all_targets:
            file.write(",".join(map(str, item)))
            file.write("\n")
    with open(output_dir + 'hairpin_scores.csv', 'w') as file:
        for item in hairpin_scores:
            file.write(",".join(map(str, item)))
            file.write("\n")
    with open(output_dir + 'snp_conflicts.csv', 'w') as file:
        for item in snp_conflicts:
            file.write(",".join(map(str, item)))
            file.write("\n")
    with open(output_dir + 'probe_arms.csv', 'wb') as file:
        pickle.dump(possible_arm_combinations_all_targets, file)
    with open(output_dir + 'tms_down.csv', 'wb') as file:
        pickle.dump(tms_up, file)
    with open(output_dir + 'tms_up.csv', 'wb') as file:
        pickle.dump(tms_down, file)
    end = time.time()
    print(end - start)
    return

if __name__ == '__main__':
    main()    
