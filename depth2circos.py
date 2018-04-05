import os, sys
import re
import argparse


# Prepare chromosome conversio ndictionary
chrDict = dict() 
chromosomeList = [str(i) for i in range(1,23)] + ['X', 'Y', 'MT']
for i in chromosomeList:
    chrDict.update({i: 'hs'+i})

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='output from samtools depth command')
    args = parser.parse_args()
    return args.input

def output_path(inputfile):
    return inputfile + '.circosData'

def parse_samtools_depth(inputfile, outputfile):
    """takes samtools depth output file as input and parse them into 
    appropriately formatted data file for Circos"""
    with open(outputfile, 'w') as f:
        # first line 
        f.write('# chr1 start1 end1 chr2 start2 end2 [options]\n')
        with open(inputfile, 'r') as g:
            for line in g:
                chr, start, depth, *args = line.strip().split()
                if chr in chrDict.keys(): # ignore MT, GLxxx contigs
                    circosChr = chrDict[chr]
                    f.write(f'{circosChr}\t{start}\t{start}\t{depth}\n')

                    
def main():
    # take input file 
    input = argument_parser()
    # prepare ouptut file path
    output = output_path(input)
    # parse input cnvkit file and output circos data file
    parse_samtools_depth(input, output)



if __name__=='__main__':
    main()

