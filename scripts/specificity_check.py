#!/usr/bin/python3
# Check probes in panel for specificity
import os
from argparse import ArgumentParser
import json
import csv
import itertools

def main():
    parser = get_arg_parser()
    args = vars(parser.parse_args())
    config_file = args["config_file"]
    panel_path_csv = args["panel_path_csv"]
    output_path = args["output_path"]
    conversion, ref_genomes_list, tile_size, min_match, min_score, target_region, outputdir, nr_cores = import_config(config_file)
    panel_arm_path_fasta = create_probe_arm_fasta(outputdir, panel_path_csv)
    aspecific_probe_list = check_specificity(conversion, ref_genomes_list, panel_arm_path_fasta, outputdir, tile_size, min_match, min_score, target_region, nr_cores)
    with open(output_path, 'w') as handle:
        writer = csv.writer(handle, delimiter = ',')
        writer.writerows(aspecific_probe_list)


def get_arg_parser():
    parser = ArgumentParser()
    parser.add_argument("-c", "--config_file", dest="config_file",
                        help="Config file with probe specifics", metavar="JSONfile")
    parser.add_argument("-g", "--reference_genome", dest="reference_genome",
                        help="Path to reference genome", metavar="REFERENCE_GENOME")
    parser.add_argument("-p", "--panel_csv", dest="panel_path_csv",
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
    target_region = config_object['specificity_specifics'][0]['target_region']
    nr_cores = config_object['max_threads']
    ref_genomes_list = []
    if conversion == 'bisulfite':
        PATH_to_reference_bisulfite_genome_folder = config_object['PATH_to_reference_bisulfite_genome_folder']
        ref_genomes_list = os.listdir(PATH_to_reference_bisulfite_genome_folder + '/CT_conversion')
        ref_genomes_list_CT = [PATH_to_reference_bisulfite_genome_folder + '/CT_conversion/'  + ref_genomes for ref_genomes in ref_genomes_list if ref_genomes[-3:] == '.fa']
        ref_genomes_list = os.listdir(PATH_to_reference_bisulfite_genome_folder + '/GA_conversion')
        ref_genomes_list_GA = [PATH_to_reference_bisulfite_genome_folder + '/GA_conversion/'  + ref_genomes for ref_genomes in ref_genomes_list if ref_genomes[-3:] == '.fa']
        ref_genomes_list = [*ref_genomes_list_CT,*ref_genomes_list_GA] 
    else:
        path_to_reference_genome = config_object['PATH_to_reference_genome_fasta']
        ref_genomes_list.append(path_to_reference_genome)
    return conversion, ref_genomes_list, tile_size, min_match, min_score, target_region, outputdir, nr_cores

def create_probe_arm_fasta(outputdir, panel_path_csv):
    arms = []
    with open(panel_path_csv,'r') as file:
        reader = csv.reader(file)
        next(reader,None)  #header
        for i,row in enumerate(reader):
            if len(row) > 2:
                u_arm = row[3]
                d_arm = row[5]
                probe_nr = i
                arms.append([str(probe_nr),u_arm,d_arm])
    panel_arm_path_fasta = outputdir + 'probe_arms.fasta'
    with open(panel_arm_path_fasta,'w') as outfile:
        writer = csv.writer(outfile,delimiter = '\n')
        for row in arms:
            writer.writerow(['>'+row[0]+'_u'])
            writer.writerow([row[1]])
            writer.writerow(['>'+row[0]+'_d'])
            writer.writerow([row[2]])
    return panel_arm_path_fasta

def check_specificity(conversion, ref_genomes_list, panel_arm_path_fasta, outputdir, tile_size, min_match, min_score, target_region, nr_cores):
    # If we design a panel for bisulfite converted DNA, blast agains 2 versions of the DNA.
    # C-to-T converted and G-to-A converted.
    # Since we design probes with prefereably no CpGs int he arms, we do not have to take into account whether CpGs are methylated or not.
    aspecific_probe_list = []
    print(ref_genomes_list)
    print(panel_arm_path_fasta)
    for i, genome in enumerate(ref_genomes_list):
        os.system('blastn -db ' + genome + ' -query ' + panel_arm_path_fasta + ' -out ' + outputdir  + 'specificity_output_' + str(i) + '.blast8 ' + '-outfmt 6 -word_size ' + str(tile_size) + ' -min_raw_gapped_score ' + str(min_score) + ' -num_threads ' + str(nr_cores))
    # Report hits of the two arms in one probe in close proximity of eachother for each of the strands of the (converted) DNA.
    # These hits can be the locations in which a product can be formed at probe capture.
    for i,genome in enumerate(ref_genomes_list):
        hits_per_probe_u = {}
        hits_per_probe_d = {}
        hit_list = []
        origin_fasta_filename = genome.split('.')[0]
        with open(outputdir + 'specificity_output_' + str(i) + '.blast8') as handle:
            reader = csv.reader(handle, delimiter = '\t')
            for row in reader:
                hit_list.append(row)
                probe_ref_gen_combination = row[0] + '.' + row[1]
                arm_type = row[0][-1]
                if arm_type == 'u':
                    if probe_ref_gen_combination in hits_per_probe_u.keys():
                        pass
                    else:
                        hits_per_probe_u[probe_ref_gen_combination] = set()
                else:
                    if probe_ref_gen_combination in hits_per_probe_d.keys():
                        pass
                    else:
                        hits_per_probe_d[probe_ref_gen_combination] = set()
        for hit in hit_list:
            arm_type = hit[0][-1] #u denotes upstream arm, d denotes downstream arm
            probe_ref_gen_combination = hit[0] + '.' + hit[1]
            locus = (int(hit[8]),int(hit[9]))
            if arm_type == 'u':
                hits_per_probe_u[probe_ref_gen_combination].add(locus)
            else:
                hits_per_probe_d[probe_ref_gen_combination].add(locus)
            
        for key in list(hits_per_probe_u.keys()):
            split_key = key.split('.')
            probe_nr = split_key[0].split('_')[0]
            chr_nr = split_key[1] #.split('_')[0]
            d_key = probe_nr + '_' + 'd' + '.' + chr_nr
            if d_key in list(hits_per_probe_d.keys()):
                locus_list_u = list(hits_per_probe_u[key])
                locus_list_d = list(hits_per_probe_d[d_key])
                differences_list = [(abs(x2-y1),(x2,y1)) for (x1,x2), (y1,y2) in itertools.product(locus_list_u,locus_list_d)]
                alignment_differences = sorted(filter(lambda x: x[0] <= target_region, differences_list))
                amount_of_alignments = len(alignment_differences)
                if amount_of_alignments > 0 :
                    print('Probe ' + str(key) + ' has ' + str(amount_of_alignments) + ' alignments')
                    # Report the genomic locations that this probe aligns to.
                    #aspecific_probe_list.append([probe_nr,chr_nr,conversion_type,value[1]] for value in alignment_differences)
                    for value in alignment_differences:
                        aspecific_probe_list.append([probe_nr,chr_nr,origin_fasta_filename,value[1][0],value[1][1]]) 
                else:
                    pass
            # we expect to have at least one value in the differences list to be low. as this will be the locus where the probe actually binds.
            # continue form here to report on the amount of loci in the genome where the probe might bind in 2 spots in close proximity.
            # The strategy now includes the situation where one arm binds in two places on the genome, which will not even result in a product at the capture step... Think abou thtis
            else:
                pass

    return aspecific_probe_list

if __name__ == '__main__':
    main()
