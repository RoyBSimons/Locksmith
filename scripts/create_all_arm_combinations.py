#!/usr/bin/python3
#extract high frequent SNPs from target VCF file.
import os
from argparse import ArgumentParser
import json

parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="filename",
                    help="open FILE", metavar="FILE")
parser.add_argument("-o", "--output", dest="outputname",
                    help="write report to FILE", metavar="OUTPUTFILE")
parser.add_argument("-c", "--config_file", dest="config_file",
                    help="Config file with probe specifics", metavar="JSONfile")

args = vars(parser.parse_args())

filename=args["filename"]
outputname=args["outputname"]
config_file=args["config_file"]
with open(config_file) as jsonFile:
    configObject = json.load(jsonFile)
    jsonFile.close()


import csv
probe_specifics=configObject['probe_specifics'][0]
min_arm_length=probe_specifics['min_arm_length']
max_arm_length=probe_specifics['max_arm_length']
min_target_length=probe_specifics['min_target_length']
max_target_length=probe_specifics['max_target_length']


#with open(filename) as inputfile:

