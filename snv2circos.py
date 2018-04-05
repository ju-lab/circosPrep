import os, sys
import re
import argparse
import vcf
import subprocess
import shlex
#import tabix
import math
# Prepare chromosome conversio ndictionary
chrDict = dict() 
chromosomeList = [str(i) for i in range(1,23)] + ['X', 'Y']
for i in chromosomeList:
    chrDict.update({i: 'hs'+i})



complementSNV_change = dict({'G>C': 'C>G', 
                            'G>T': 'C>A', 
                            'G>A': 'C>T', 
                            'A>T': 'T>A', 
                            'A>G': 'T>C', 
                            'A>C': 'T>G'})
SNV_change_color = dict({'C>G': 'black', 'C>A': 'blue', 'C>T': 'red', 'T>A': 'magenta', 'T>C': 'yellow', 'T>G': 'green'})


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='SNV file in vcf or vcf.gz file')
    parser.add_argument('-c', '--caller', required=True, choices=['mutect', 'strelka', 'freebayes'], help='SNV caller used')
    parser.add_argument('-y', '--yaxis', required=True, choices=['vaf', 'rainfall'], help='value of Y-axis to represent')
    args = parser.parse_args()
    return args.input, args.caller, args.yaxis

def output_path(inputfile, yaxis):
    return inputfile +'.' + yaxis +  '.circosData'

def copy_vcfgz(vcfgzfile):
    '''if input is a compressed vcf.gz file, create a temporary copy and unzip it'''

    tempfile = re.sub(string=vcfgzfile, pattern=r'.vcf.gz$', repl='.temp.vcf.gz')
    
    copycmd = f'cp {vcfgzfile} {tempfile}'
    subprocess.call(shlex.split(copycmd))
    unzipcmd = f'gunzip {tempfile}'
    subprocess.call(shlex.split(unzipcmd))
    uncompressed_copy = re.sub(string=tempfile, pattern=r'.gz$', repl='')
    print(uncompressed_copy)
    return uncompressed_copy

def parse_snv_vaf_vcf(inputfile, outputfile, caller):
    """takes structural variation vcf file as input and parse them into 
    appropriately formatted data file for Circos"""
    if inputfile.endswith('.gz'):
        vcffile = copy_vcfgz(inputfile)

    else:
        vcffile = inputfile
    with open(vcffile, 'r') as g:
        with open(outputfile, 'w') as f:
             for line in g:
                distance_filter = False # distance filter, if SV too small don't show
                if not line.startswith('#'): # skip vcf header
                    
                    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, TUMOR, NORMAL = line.strip().split()
                    if CHROM in chrDict.keys(): # only consider major chromosome contigs
                        #print(TUMOR.split(';'))
                        vaf = TUMOR.split(':')[4] # FA column for mutect

                        snvCircosChrom = chrDict[CHROM]
                         
                        basechange = f'{REF}>{ALT}'
                        if REF in ['G', 'A']:
                            basechange = complementSNV_change[basechange]
                        
                        color = SNV_change_color[basechange]


                        f.write(f'{snvCircosChrom}\t{POS}\t{POS}\t{vaf}\tcolor={color}\n')


def parse_snv_rainfall_vcf(inputfile, outputfile):
    """prepare for circos rainfall plot from a given vcf file"""
    if inputfile.endswith('.gz'):
        vcffile = copy_vcfgz(inputfile)

    else:
        vcffile = inputfile

    with open(vcffile, 'r') as g:
        prev_location = 0
        prev_chr = '0'
        with open(outputfile, 'w') as f:
            for line in g:
                if not line.startswith('#'):
                    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, TUMOR, NORMAL = line.strip().split() 
                    if CHROM in chrDict.keys():
                        snvCircosChrom = chrDict[CHROM]
                        basechange = f'{REF}>{ALT}'
                        if REF in ['G', 'A']:
                            basechange = complementSNV_change[basechange]

                        color = SNV_change_color[basechange]

                        if prev_chr == CHROM:
                            
                            distance = int(POS) - int(prev_location)
                            prev_location = int(POS)
                            prev_chr = CHROM        
                        else:
                            distance = 1
                            prev_location = int(POS)
                            prev_chr = CHROM
                        log_distance = math.log10(distance)

                        f.write(f'{snvCircosChrom}\t{POS}\t{POS}\t{log_distance}\tcolor={color}\n')
                    else:
                        pass # if on extra contigs other than the major ones then skip  

    # remove temporary copy if vcf.gz file was the input
    if inputfile.endswith('.vcf.gz'):
        subprocess.call(shlex.split('rm -rf ' + vcffile))

    return 0

def main():
    # take input file 
    input, caller, yaxis = argument_parser()
    # prepare ouptut file path
    output = output_path(input, yaxis)
    # parse input cnvkit file and output circos data file
    if yaxis == 'vaf':
        parse_snv_vaf_vcf(input, output, caller)
    elif yaxis == 'rainfall':
        parse_snv_rainfall_vcf(input, output)


    print(f'Output Circos Data file is written: {output}')


if __name__=='__main__':
    main()

