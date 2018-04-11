import os, sys
import re
import argparse
import vcf
import subprocess
import shlex
import tabix

# Prepare chromosome conversio ndictionary
chrDict = dict() 
chromosomeList = [str(i) for i in range(1,23)] + ['X', 'Y']
for i in chromosomeList:
    chrDict.update({i: 'hs'+i})

svtype2color = dict({'BND': 'black', 'DEL': 'yellow', 'DUP': 'blue', 'INV': 'orange', 'INS': 'green', 'TRA': 'black'})

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='input smoothened CNV file')
    parser.add_argument('-o', '--output_dir', required=False, default=os.getcwd(), help='folder to write output')
    args = parser.parse_args()
    return args.input, args.output_dir


def output_path(inputfile, output_dir):
    return os.path.join(output_dir, os.path.basename(inputfile) + '.circosData')


def parse_smoothened_cnv_file(inputfile, output_dir):
    """takes structural variation vcf file as input and parse them into 
    appropriately formatted data file for Circos"""
    outputfile = output_path(inputfile, output_dir)

    with open(outputfile, 'w') as g:
        with open(inputfile, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    chromosome, position, _, _, cnv = line.strip().split()
                    if cnv != 'NA':
                        
                        if chromosome in chrDict.keys():
                            hschr = chrDict[chromosome]
                            information = f'{hschr}\t{position}\t{position}\t{cnv}\n'
                            g.write(information)



    
    return 0                    
        
                    
def main():
    # take input file 
    input, output_dir = argument_parser()
    parse_smoothened_cnv_file(input, output_dir)



if __name__=='__main__':
    main()

