#!/usr/bin/python3
# Check probes in panel for specificity
import os
from argparse import ArgumentParser
import json
import csv
import itertools
import collections

def main():
    parser = get_arg_parser()
    args = vars(parser.parse_args())
    config_file = args["config_file"]
    panel_path_csv = args["panel_path_csv"]
    output_path = args["output_path"]
    conversion, ref_genomes_list, tile_size, min_match, min_score, target_region, outputdir, nr_cores = import_config(config_file)
    panel_arm_path_fasta = create_probe_arm_fasta(outputdir, panel_path_csv)
    aspecific_probe_list = check_specificity(conversion, ref_genomes_list, panel_arm_path_fasta, outputdir, tile_size, min_match, min_score, target_region, nr_cores)
    #write hits to file
    with open(output_path, 'w') as handle:
        writer = csv.writer(handle, delimiter = ',')
        writer.writerows(aspecific_probe_list)
    #append hits to chosen panel output file
    append_nr_of_hits_to_panel(panel_path_csv,aspecific_probe_list)


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
    # Since we design probes with prefereably no CpGs int he arms, we do not have to take into account whether CpGs are methylated or not. For the probes that do have a CpG or SNP in the arms, this could result in a single mismatch when aligning to the reference genome. These probes should still report back a hit since we look for a lower score than perfect match.
    aspecific_probe_list = []
    print(ref_genomes_list)
    print(panel_arm_path_fasta)
    for i, genome in enumerate(ref_genomes_list):
        # Define alignment: bottom strand aslignemnt at CT genome and top strand alignment at GA genome as these are the strand which should correspond to the arms of the padlock probes.
        if 'GA' in genome:
            strandedness = 'plus'
            print('GA genome strandedness is ' + str(strandedness))
        elif 'CT' in genome:
            strandedness = 'minus'
            print('CT genome strandedness is ' + str(strandedness))
        else:
            strandedness = 'both'
            print('genome strandedness is ' + str(strandedness))
        os.system('blastn -db ' + genome + ' -query ' + panel_arm_path_fasta + ' -out ' + outputdir  + 'specificity_output_' + str(i) + '.blast8 ' + '-outfmt 6 -word_size ' + str(tile_size) + ' -min_raw_gapped_score ' + str(min_score) + ' -num_threads ' + str(nr_cores) + ' -strand ' + strandedness + ' -dust no -soft_masking false')
    # Report hits of the two arms in one probe in close proximity of eachother for each of the strands of the (converted) DNA.
    # These hits can be the locations in which a product can be formed at probe capture.
    # Since we collate the hits per reference genome, we already restrict for arms binding on the same strand.
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
                # Report the hits for which the two arms are in a distance closer than the target_region parameter
                # To implement: Directionality check
                    #Could this combination of aligned arms result in a circularized product after capture?
                if 'GA' in genome:
                    differences_generator = ((abs(x2-y1),(x2,y1)) for (x1,x2), (y1,y2) in itertools.product(locus_list_u,locus_list_d) if y1 < x2)
                    # only keep hits at which y1 is lower than x2. --> The downstream arm binds downstream of the upstream arm and the probe can be elongated and circularized.
                elif 'CT' in genome:
                    # only keep hits at which x2 is lower than y1. --> The downstream arm binds downstream of the upstream arm and the probe can be elongated and circularized.
                    differences_generator = ((abs(x2-y1),(x2,y1)) for (x1,x2), (y1,y2) in itertools.product(locus_list_u,locus_list_d) if y1 > x2)

                alignment_differences = sorted(filter(lambda x: x[0] <= target_region, differences_generator))
                amount_of_alignments = len(alignment_differences)

                # To implement: Strand check
                    # --> Check whether I should blast against a ref genome in which CpG assumed to be converted/unconverted (4 ref genomes instead of 2)
                    # we are now assuming that all CpGs are unmethylated. C-to-T. 

                # To implement: Location check
                    # Check whether the location of the CpG is in between the two arms of the hit


    
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

def append_nr_of_hits_to_panel(panel_path_csv,aspecific_probe_list):
    # Function to add amount of hits and inclusion of correct hit to the final chosen panel file.
    # Create a dictionary counter with the amount of hits per probe.
    nr_of_probes_with_more_than_one_hit = 0
    spec_hit_list = [item[0] for item in aspecific_probe_list]
    probe_spec_dict = collections.Counter(spec_hit_list)
    with open(panel_path_csv,'r') as file:
        reader = csv.reader(file)
        header = next(reader,None)  # Skip header
        new_panel_path_csv = panel_path_csv.split('.')[0] + '_specificity' + '.csv'
        with open(new_panel_path_csv, 'w') as file_out:
            writer = csv.writer(file_out)
            writer.writerow(header)
            for i,row in enumerate(reader):
                if len(row) <= 2:
                    # No probe designed for this target
                    pass
                else:
                    # Report hits per probe
                    if str(i) in probe_spec_dict:
                        count = probe_spec_dict[str(i)]
                        row[-3] = count
                        if count > 1:
                            nr_of_probes_with_more_than_one_hit += 1
                    else:
                        row[-3] = 0

                    hits_for_target_i_list = [item for item in aspecific_probe_list if int(item[0]) == i]
                    hit_boolean = 0
                    if row[9] == '+':
                        start = int(row[7].split(':')[1].split('-')[0])+1
                        end = int(row[6].split(':')[1].split('-')[1])+1
                    elif row[9] == '-':
                        end = int(row[7].split(':')[1].split('-')[1])+1
                        start = int(row[6].split(':')[1].split('-')[0])+1
                    # Check whether the intended target is included in the hits.
                    for hit in hits_for_target_i_list:
                        hit_u = hit[3]
                        hit_d = hit[4]
                        if hit_u >= start and hit_u <= end and hit_d >= start and hit_d <= end:
                            hit_boolean = 1
                            break
                    # Return a 1 if the intented target is in included in the hits.
                    row[-2] = hit_boolean
                writer.writerow(row)
    #Add specificity parameters to the final chosen panel file.
    os.system('rm ' + panel_path_csv)
    os.system('mv ' + new_panel_path_csv + ' ' + panel_path_csv)
    print('There are ' + str(nr_of_probes_with_more_than_one_hit) + ' probes which target more than one location.' )
    return


if __name__ == '__main__':
    main()
