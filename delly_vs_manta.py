"""
Combined Structural variation callset from both DELLY and MANTA

"""

import argparse
import cyvcf2
import re
import time
import datetime
import os
import shutil
import sys
import vcf


# Prepare chromosome conversion dictionary for circos
chrDict = dict()
chromosomeList = [str(i) for i in range(1,23)] + ['X', 'Y']
for i in chromosomeList:
    chrDict.update({i: 'hs'+i})

# make sure this matches with `sv2circos.py` for consistency
svtype2color = dict({'BND': 'black', 'DEL': 'yellow', 'DUP': 'blue', 'INV': 'orange', 'INS': 'green', 'TRA': 'black'})

class Position():
    ''' python class for handling genomic positions
    0-based
    '''
    def __init__(self, chromosome, start, end, is_bp=None, clipped_reads=None):
        self.chromosome = chromosome
        self.start = start
        self.end = end

    def __repr__(self):
        return str(self.chromosome) + ":" + str(self.start) + '-' + str(self.end)

    def __str__(self):
        return str(self.chromosome) + ":" + str(self.start) + '-' + str(self.end)

    
    @classmethod
    def fromstring(cls, position_string):
        chromosome = position_string.split(':')[0]
        start = int(position_string.split(':')[1].split('-')[0])
        end = int(position_string.split(':')[1].split('-')[1])
        return Position(chromosome, start, end)

    @staticmethod
    def overlap(position1, position2):
        '''true if position1 and position2 has more than 1 overlapping base'''
        try:
            if isinstance(position1, Position) and isinstance(position2, Position):
                if position1.chromosome == position2.chromosome:
                    if min(position1.end, position2.end) > max (position1.start, position2.start):
                        return True
                    else:
                        return False
                else:
                    return False  # cannot compare if two positions are in different chromosome
            else:
                return None # has to be Posiiton class.
        except:
            Exception
    
    def extend(self, direction, basepairs):
        """extends objects in by specified base pairs, either upstream, downstream, or both"""
        if direction=="up":
            return Position(self.chromosome, max(0, self.start-basepairs), end)
        elif direction=="down":
            return Position(self.chromosome, self.start, self.end + basepairs)
        elif direction=="both":
            return Position(self.chromosome, max(0, self.start - basepairs), self.end + basepairs)
        else:
            print('direction has to be either up, down, or both')
            raise ValueError


def argument_parser():
    """parses argument passed on from command line"""
    parser = argparse.ArgumentParser(
        description='Find fusion evidence from WGS results')

    parser.add_argument('-d', '--distance_threshold', required=True, type=int)
    parser.add_argument('-o', '--outputDIR', required=False, default=os.getcwd())
    parser.add_argument('--delly', required=True, help='Delly SV vcf file')
    parser.add_argument('--manta', required=True, help='Manta SV vcf file')
    parser.add_argument('-s', '--sampleName', required=True, help='Sample Name to be used for output .circosData file')
    args = vars(parser.parse_args())

    distance_threshold = args['distance_threshold']
    outputDIR = args['outputDIR']
    delly = args['delly']
    manta = args['manta']
    sampleName = args['sampleName']
    return outputDIR, distance_threshold, delly, manta, sampleName


def vcf2SVPosition(vcf_file):
    ''' create a generator of Position object of break points from a given vcf file
    added Manta/Lumpy support Dec 21 2017, and commented out old function above
    March 3 2018 edit. Created Dictionary to call each SV types separately. 
    '''
    bp_dict = dict({'BND':[], 'DUP':[], 'INS':[], 'DEL':[], 'INV':[]})

    for variant in cyvcf2.VCF(vcf_file):
        # if variant.FILTER == None:
        variant_type = variant.INFO.get('SVTYPE')
        if variant_type  == "BND":
            # this one can be used for any SV VCF with BND type
            bnd_pos= re.search(string=variant.ALT[0], pattern=r'[a-zA-Z]*[0-9]*:[0-9]+').group(0)
            bnd_chrom, bnd_pos = bnd_pos.split(':')
            bp2String = f'{bnd_chrom}:{bnd_pos}-{int(bnd_pos) + 1}'
        elif variant_type == "TRA":
            # this one is specific for Delly v0.7.6 annotation
            bp2String = f"{variant.INFO.get('CHR2')}:{variant.INFO.get('END')}-{variant.INFO.get('END') + 1}"
        else:
            bp2String = f'{variant.CHROM}:{variant.INFO.get("END")}-{variant.INFO.get("END") + 1}'

        bp1 = Position.fromstring(f'{variant.CHROM}:{variant.POS}-{variant.POS + 1}')
        bp2 = Position.fromstring(bp2String)

        bp_dict[variant_type].append((bp1, bp2))
        
    return bp_dict

def vcf2SVPosition_select(vcf_file, select):
    ''' create a generator of Position object of break points from a given vcf file
    added Manta/Lumpy support Dec 21 2017, and commented out old function above

    from a vcf file with multiple structural variant types, select a single variant type for result list
    select can be either BND, TRA, INV, DEL, DUP
    Updated March 3 to use vcf2SVPosition
    '''
            
    return vcf2SVPosition(vcf_file)[select]
 


# this dictionary is a map that converts bp support overlap -> number of supports
# each edge has weight of 1, 2, 4, 8 and sum of these weights are converted into 
# support level to see if there is unique two overlap  
# used in 'support_breakpoint' function      
support_dict = ({
    0:0,
    1:1,
    2:1,
    3:1,
    4:1,
    5:1,
    6:2,
    7:2,
    8:1,
    9:2,
    10:1,
    11:1,
    12:1, 
    13:2,
    14:2,
    15:2
})


def support_breakpoint(upstream_bp1, downstream_bp1, upstream_bp2, downstream_bp2, distance_threshold):
    '''takes in the two breakpoints of fusion and sv records to compare
    compares the two bed files each (4 total) and see if there's support for both breakpoints from SV, or just one, or none'''
    bp_support_status = dict({1: 0, 2:0, 4:0, 8:0})
    edge_number = 1 
    for bp1 in [upstream_bp1, downstream_bp1]:       
        for bp2 in [upstream_bp2, downstream_bp2]:
            # check whether s_bp is within 'dist' base of f_bp
            result = Position.overlap(bp1.extend('both', distance_threshold), bp2) 
            if result:
                bp_support_status[edge_number] = 1
            edge_number *= 2 # multiply edge number by 2 to get binary representation

    
    return support_dict[sum([key for key, value in bp_support_status.items() if value == 1])]          # this is the number of break point support
    

def compare_breakpoints(vcf1, vcf2, distance_threshold, circosDataFile):
    vcf1_positions = vcf2SVPosition(vcf1)
    vcf2_positions = vcf2SVPosition(vcf2)
    with open(circosDataFile, 'w') as f:
        for svtype in ['BND', 'INV', 'DUP', 'DEL', 'INS']:
            print(svtype)
            color = svtype2color[svtype]
            sv_list = []
            for vcf1_sv in vcf1_positions[svtype]:
                upstream_vcf1_bp, downstream_vcf1_bp = vcf1_sv
                status = False

                for vcf2_sv in vcf2_positions[svtype]:
                    upstream_vcf2_bp, downstream_vcf2_bp = vcf2_sv
                    support = (support_breakpoint(upstream_vcf1_bp, downstream_vcf1_bp, upstream_vcf2_bp, downstream_vcf2_bp, distance_threshold))
                    if support ==2 :
                        status = True
            
                if status == True:
                    print(vcf1_sv)
                    bp1, bp2 = vcf1_sv
                    bp1circosChrom = chrDict[bp1.chromosome]
                    bp2circosChrom = chrDict[bp2.chromosome]
                    f.write(f"{bp1circosChrom}\t{bp1.start}\t{bp1.start}\t{bp2circosChrom}\t{bp2.start}\t{bp2.start}\tcolor={color}\n")
    

    return 0

def main():
    outputDIR, distance_threshold, delly, manta, sampleName = argument_parser()

    circosDataFile = os.path.join(outputDIR, sampleName + '.circosData')
    compare_breakpoints(delly, manta, 1000, circosDataFile)
    print(f'Output Circos Data file is written: {circosDataFile}')
    return 0

if __name__=='__main__':
    main()

