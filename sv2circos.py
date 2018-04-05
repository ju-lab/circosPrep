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
    parser.add_argument('-i', '--input', required=True, help='SV file in vcf or vcf.gz file')
    parser.add_argument('-c', '--caller', required=True, choices=['delly', 'manta'], help='SV caller used')
    parser.add_argument('-d', '--distance_threshold', required=False, default=100000, type=int, help='Distance threshold within same chromosome to display')
    args = parser.parse_args()
    return args.input, args.caller, args.distance_threshold

def output_path(inputfile):
    return inputfile + '.circosData'

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

def parse_sv_vcf(inputfile, outputfile, caller, distance_threshold):
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
                    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, *args = line.strip().split()
                    svtype = re.sub(string=re.search(string=INFO, pattern=r'SVTYPE=[A-Z]+').group(0), pattern='SVTYPE=', repl='')
                    color = svtype2color[svtype]

                    if CHROM in chrDict.keys(): # only consider major chromosome contigs
                        if caller=='delly':
                            bp2Chr = re.sub(string=re.search(string=INFO, pattern=r'CHR2=[A-Za-z0-9]*').group(0), pattern='CHR2=', repl='')
                            bp2POS = re.sub(string=re.search(string=INFO, pattern=r'END=[0-9]+').group(0), pattern='END=', repl='')
                    
                        elif caller=='manta':
                            if re.search(string=INFO, pattern=r'SVTYPE=BND'): # if sv type is bnd
                                bnd_pos= re.search(string=ALT, pattern=r'[a-zA-Z]*[0-9]*:[0-9]+').group(0)
                                bp2Chr, bp2POS = bnd_pos.split(':')
                            else: # if variant is on the same chromosome (i.e. not BND)
                                
                                svlen = int(re.sub(string=re.search(string=INFO, pattern=r'SVLEN=-*[0-9]+').group(0), pattern=r'SVLEN=', repl=''))
                                bp2Chr=CHROM
                                bp2POS = int(POS) + svlen                        

                        bp1CircosChrom = chrDict[CHROM]
                        bp2CircosChrom = chrDict[bp2Chr] 
                        
                        # if Breakpoints on different chromosome, then display
                        # if on same chromosome, don't show unless it is greater than threshold
                        if CHROM == bp2Chr:
                            if abs(int(POS) - int(bp2POS)) > distance_threshold:
                                distance_filter = True
                        else:
                            distance_filter = True

                        if distance_filter == True:
                            f.write(f'{bp1CircosChrom}\t{POS}\t{POS}\t{bp2CircosChrom}\t{bp2POS}\t{bp2POS}\tcolor={color}\n')

    
    # remove temporary copy if vcf.gz file was the input
    if inputfile.endswith('.vcf.gz'):
        subprocess.call(shlex.split('rm -rf ' + vcffile))



                    
        
                    
def main():
    # take input file 
    input, caller, distance_threshold = argument_parser()
    # prepare ouptut file path
    output = output_path(input)
    # parse input cnvkit file and output circos data file
    parse_sv_vcf(input, output, caller, distance_threshold)

    print(f'Output Circos Data file is written: {output}')


if __name__=='__main__':
    main()

