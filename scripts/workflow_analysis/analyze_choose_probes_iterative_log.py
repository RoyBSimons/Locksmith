#!/usr/bin/python3
#analyzing log file of choose probes iterative rule
#python ~/opt/Locksmith/scripts/workflow_analysis/analyze_choose_probes_iterative_log.py -i /media/disk3/roy/output_Locksmith/21_dec_horvath* -o dimers_per_iteration.png
import csv
from argparse import ArgumentParser
import matplotlib
import matplotlib.pyplot as plt
import os

def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input",nargs='+', dest="input_name",
                        help="input log files of choose probes iterative rule", metavar="INPUTNAME")
    parser.add_argument("-o", "--output", dest="output_name",
                        help="Fasta file which will contain the chosen probes as output", metavar="OUTPUTNAME")
    
    args=vars(parser.parse_args())
    input_files=args['input_name']
    output_file=args['output_name']
    conflicts_per_file=[]
    labels=[]
    
    for file in input_files:
        conflict_list=[]
        label=file.split('/')[5].split('_')[3]
        labels.append(label)
        with open(file,'r') as handle:
            reader=handle.readlines()
            for row in reader:
                if 'out of ' in row:
                    nr_of_conflicts=int(row.rstrip().split(' out of ')[0])
                    conflict_list.append(nr_of_conflicts)
        conflicts_per_file.append(conflict_list)
    
    colors=['r','g','b','c','m','y','k','w']
    font = {'weight' : 'bold',
            'size'   : 20}
    plt.rc('font',**font)
    fig, ax = plt.subplots()
    for i,conflict_list in enumerate(conflicts_per_file):
        y=conflict_list
        x=list(range(0,len(y)))
        y2=create_min_list(y)
        print(min(y))
        plt.plot(x,y,color=colors[i],alpha=0.5,linewidth=3)
        plt.plot(x,y2,color=colors[i])
    plt.legend(labels)
    leg = ax.get_legend()
    for i in range(0,len(conflicts_per_file)):
        leg.legendHandles[i].set_color(colors[i])
    plt.xlabel('Iterations')
    plt.ylabel('Nr of conflicts')
    fig.set_size_inches(15, 8)
    fig.savefig(output_file)
    os.system('display '+output_file)


def create_min_list(y):
    min_list=[]
    minimum=y[0]
    for val in y:
        new_min=min(minimum,val)
        min_list.append(new_min)
        minimum=new_min
    return min_list

if __name__ == '__main__':
    main()

