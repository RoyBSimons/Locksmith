#!/usr/bin/python3
#creates a bed file from cglist
#python ../scripts/cg_illumina_to_hg38_bed.py -c cpglist_horvath.txt -i ~/Documents/Illumina_450K/HM450.hg38.manifest.tsv -o target_list_test.bed -p ~/opt/Locksmith/data/hg38/* -e .bed
import os
from argparse import ArgumentParser
import glob
import csv
parser = ArgumentParser()
parser.add_argument("-c", "--cgid_list", dest="cglistfile",
                    help="open FILE", metavar="FILE1")
parser.add_argument("-i", "--illumina_cgid_list", dest="illumina_cgid_listfile",
                    help="open FILE", metavar="FILE2")
parser.add_argument("-o", "--output", dest="outputname",
                    help="write report to FILE", metavar="OUTPUTFILE")
parser.add_argument("-p", dest='folder_path', nargs='+', help='Path of a file or a folder of files.')
parser.add_argument('-e', '--extension', default='', help='File extension to filter by.')
args = parser.parse_args()
# Parse paths
full_paths = [os.path.join(os.getcwd(), path) for path in args.folder_path]
files = set()

for path in full_paths:
    if os.path.isfile(path):
        fileName, fileExt = os.path.splitext(path)
        if args.extension == '' or args.extension == fileExt:
            files.add(path)
    else:
        pass
        #        if (args.recursive):
#            full_paths += glob.glob(path + '/*')
chr_list_cpg=[]
loc_list_cpg=[]
cpg_list_cpg=[]
for f in files:
    with open(f) as handle:
        reader=csv.reader(handle,delimiter='\t')
        for row in reader:
            chr_list_cpg.append(row[0])
            loc_list_cpg.append(row[1])
            cpg_list_cpg.append(row[3])
print(full_paths)            
cglistfile=args.cglistfile
illumina_cgid_listfile=args.illumina_cgid_listfile
outputname=args.outputname
import csv
cglist=[]

with open(cglistfile) as inputfile1:
    input_reader1=csv.reader(inputfile1)
    for row in input_reader1:
        cglist.append(row)

chr_list=[]
loc_list=[]
cpg_list=[]
strand_list=[]
with open(illumina_cgid_listfile) as inputfile:
    input_reader=csv.reader(inputfile,delimiter="\t")
    for row in input_reader:
        if [row[4]] in cglist:
            chr_list.append(str(row[0]))
            loc_list.append(int(row[1])+1)
            strand_list.append(str(row[3]))
            i_chr=[i for i,x in enumerate(chr_list_cpg) if x == row[0]]
            i_loc=[i for i,x in enumerate(loc_list_cpg) if x == str(int(row[1])+1)]
            index=list(set(i_chr).intersection(i_loc))
            print(index)
#            print(i_chr)
#            print(i_loc)
            print(row[0])
            print(chr_list_cpg[0])
            cpg_list.append(cpg_list_cpg[index[0]]) #place where loc and chr are both correct
        else:
            pass
with open(outputname,'w') as handle:
    writer=csv.writer(handle,delimiter='\t')
    for i,loc in enumerate(loc_list):
        writer.writerow([chr_list[i],str(loc),str(loc+1),cpg_list[i],strand_list[i]])

    #os.system("grep chr"+str(chr_list[i])+"'\t'"+str(loc)+"'\t' ./data/CpG_bed_nomenclature/CpG_hg19_chr"+str(chr_list[i])+".bed >> "+outputname+" -h")
