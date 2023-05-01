#!/usr/bin/python3
# Check probes in panel for specificity
import os
from argparse import ArgumentParser
import json
import csv
from itertools import combinations

def main():
    parser = get_arg_parser()
    args = vars(parser.parse_args())
    config_file = args["config_file"]
    panel_path_fasta = args["panel_path_fasta"]
    output_path = args["output_path"]
    conversion, ref_genomes_list, tile_size, min_match, min_score, outputdir = import_config(config_file)
    aspecific_probe_list = check_specificity(conversion, ref_genomes_list, panel_path_fasta, outputdir, tile_size, min_match, min_score)
    with open(output_path, 'w') as handle:
        writer = csv.writer(handle, delimiter = '\t')
        writer.writerows(aspecific_probe_list)


def get_arg_parser():
    parser = ArgumentParser()
    parser.add_argument("-c", "--config_file", dest="config_file",
                        help="Config file with probe specifics", metavar="JSONfile")
    parser.add_argument("-g", "--reference_genome", dest="reference_genome",
                        help="Path to reference genome", metavar="REFERENCE_GENOME")
    parser.add_argument("-p", "--panel_fasta", dest="panel_path_fasta",
                        help="Path to FASTA file of panel", metavar="FASTA_PANEL")
    parser.add_argument("-o", "--panel_output", dest="output_path",
                        help="Path for output panel file", metavar="OUTPUT_PANEL")
    return parser

def import_config(config_file):
    # Import parameters set in configuration file
    with open(config_file) as jsonFile:
        config_object = json.load(jsonFile)
        jsonFile.close()
    conversion = config_object['probe_specifics'][0]['conversion']
    outputdir = config_object['output_directory']+'/'
    tile_size = config_object['specificity_specifics'][0]['tile_size']
    min_match = config_object['specificity_specifics'][0]['min_match']
    min_score = config_object['specificity_specifics'][0]['min_score']
    ref_genomes_list = []
    if conversion == 'bisulfite':
        PATH_to_reference_bisulfite_genome_folder = config_object['PATH_to_reference_bisulfite_genome_folder']
        ref_genomes_list.append(PATH_to_reference_bisulfite_genome_folder + 'CT_conversion/genome_mfa.CT_conversion.fa')
        ref_genomes_list.append(PATH_to_reference_bisulfite_genome_folder + 'GA_conversion/genome_mfa.GA_conversion.fa')
    else:
        path_to_reference_genome = config_object['PATH_to_reference_genome_fasta']
        ref_genomes_list.append(path_to_reference_genome)
    return conversion, ref_genomes_list, tile_size, min_match, min_score, outputdir

def check_specificity(conversion, ref_genomes_list, panel_path_fasta, outputdir, tile_size, min_match, min_score):
    # If we design a panel for bisulfite converted DNA, blast agains 2 versions of the DNA.
    # C-to-T converted and G-to-A converted.
    # Since we design probes with prefereably no CpGs int he arms, we do not have to take into account whether CpGs are methylated or not.
    for i,genome in enumerate(ref_genomes_list):
        os.system('blat ' + genome + ' ' + panel_path_fasta + ' ' + outputdir  + 'specificity_output_' + str(i) + '.blast8 -tileSize=' + str(tile_size) + ' -out=blast8 -minMatch=' + str(min_match) + ' -minScore=' + str(min_score))
    hit_list = []
    
    # Read all hits
    hits_per_probe = {}
    for i,genome in enumerate(ref_genomes_list):
        with open(outputdir + 'specificity_output_' + str(i) + '.blast8') as handle:
            reader = csv.reader(handle, delimiter = '\t')
            for row in reader:
                hit_list.append(row)
                probe_ref_gen_combination = row[0] + '.' + row[1]
                if probe_ref_gen_combination in hits_per_probe.keys():
                    pass
                else:
                    hits_per_probe[probe_ref_gen_combination] = set()
    for hit in hit_list:
        probe_ref_gen_combination = hit[0] + '.' + hit[1]
        locus = int(hit[8])
        hits_per_probe[probe_ref_gen_combination].add(locus)
    
    aspecific_probe_list = []
    for key in hits_per_probe.keys():
        split_key = key.split('.')
        probe_nr = split_key[0]
        chr_nr = split_key[1].split('_')[0]
        conversion_type = split_key[1].split('_')[1]
        locus_list = list(hits_per_probe[key])
        differences_list = [(abs(x-y),(x,y)) for (x), (y) in combinations(locus_list, 2)]
        alignment_differences = sorted(filter(lambda x: x[0] <=500, differences_list))
#        print(alignment_differences)
        amount_of_alignments = len(alignment_differences)
        if amount_of_alignments > 2 :
            print('Probe ' + str(key) + ' has ' + str(amount_of_alignments) + ' alignments')

            # Report the genomic locations that this probe aligns to.
            aspecific_probe_list.append([probe_nr,chr_nr,conversion_type,value[1]] for value in alignment_differences)
             
        else:
            pass
        # we expect to have at least one value in the differences list to be low. as this will be the locus where the probe actually binds.
        # continue form here to report on the amount of loci in the genome where the probe might bind in 2 spots in close proximity.
        # The strategy now includes the situation where one arm binds in two places on the genome, which will not even result in a product at the capture step... Think abou thtis

    return aspecific_probe_list

if __name__ == '__main__':
    main()
