#!/usr/bin/python3
# Create bed file with chr headers
import os
from argparse import ArgumentParser
import csv


def main():
    parser = get_arg_parser()  # Parse input to variables
    args = vars(parser.parse_args())
    acc_nr_to_chrom_nr_file = args["acc_2_chrom"]
    bed_path = args["bed_path"]
    output_file = args["output_file"]
    rewrite_bed_file_from_accession_nr_to_chromosome_nr(bed_path, acc_nr_to_chrom_nr_file, output_file)
    print('happened')

def  obtain_acc_nr_and_chrom_nr_lists(acc_nr_to_chrom_nr_file):
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


def rewrite_bed_file_from_accession_nr_to_chromosome_nr(bed_path, acc_nr_to_chrom_nr_file, output_file):

    # This conversion is needed for couting the frequent SNPs in the probe arms.
    chrom_nr_list, acc_nr_list = obtain_acc_nr_and_chrom_nr_lists(acc_nr_to_chrom_nr_file)
    bed_info_list = []
    with open(bed_path) as handle:
        reader = csv.reader(handle, delimiter='\t')
        for row in reader:
            chrom = row[0][3:]
            bed_info = [acc_nr_list[chrom_nr_list.index(chrom)], row[1], row[2], row[3]]
            bed_info_list.append(bed_info)
    new_path = output_file
    print(new_path)
    with open(new_path, 'w') as output_handle:
        writer = csv.writer(output_handle, delimiter='\t')
        writer.writerows(bed_info_list)
    return new_path


def get_arg_parser():
    # parse all files
    parser = ArgumentParser()
    parser.add_argument("-b", "--bed", dest="bed_path",
                        help="Bed file containing the targets", metavar="BED")
    parser.add_argument("-a", "--acc_2_chrom", dest="acc_2_chrom",
                        help="Tab delimiter file containing the chromosome nr and the accompanied accesion nr "
                             "used in the ftp database used for SNP detection", metavar="ACC2CHROM")
    parser.add_argument("-o", "--output", dest="output_file",
                        help="Path of chr version of bed file", metavar="output_file")
    return parser
if __name__ == '__main__':
    main()
