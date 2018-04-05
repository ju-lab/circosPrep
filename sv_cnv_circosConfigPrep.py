"""
Combined Structural variation callset from both DELLY and MANTA

"""

import argparse
import re
import os
import subprocess
import shlex


def argument_parser():
    """parses argument passed on from command line"""
    parser = argparse.ArgumentParser(
        description='Writes a configuration file for Circos Plot from .circosData inputs')


    parser.add_argument('--cnv', required=True, help='CNV circosData')
    parser.add_argument('--sv', required=True, help='SV circosData')
    parser.add_argument('--sampleName', required=True, help='Sample name prefix for output configuration file')
    parser.add_argument('-o', '--outputDIR', required=False, default=os.getcwd())
    parser.add_argument('-t', '--template', required=True, help='Template configuration file')

    args = vars(parser.parse_args())

    outputDIR = args['outputDIR']
    sampleName = args['sampleName']
    sv = args['sv']
    cnv = args['cnv']
    template = args['template']
    return sv, cnv, outputDIR, sampleName, template


def create_configuration(template, sv, cnv, circosConfig, circosPlot):
    """read in template file, 
    substitue <SVCIRCOS>, <CNVCIRCOS>, <OUTPUTFILE> to respective file names"""
    with open(template, 'r') as f:
        with open(circosConfig, 'w') as g:
            for line in f:
                if re.search(string=line, pattern=r'<OUTPUTFILE>') != None:
                    replaced_line = re.sub(string=line, pattern=r'<OUTPUTFILE>', repl=circosPlot)
                elif re.search(string=line, pattern=r'<SVCIRCOS>') != None:
                    replaced_line = re.sub(string=line, pattern=r'<SVCIRCOS>', repl=sv)
                elif re.search(string=line, pattern=r'<CNVCIRCOS>') != None:
                    replaced_line = re.sub(string=line, pattern=r'<CNVCIRCOS>', repl=cnv)
                else:
                    replaced_line = line 
            
                g.write(replaced_line)



def main():
    sv, cnv, outputDIR, sampleName, template = argument_parser()

    circosConfig= os.path.join(outputDIR, sampleName + '.conf')
    circosPlot = os.path.join(outputDIR, sampleName + '.circosPlot.png')
    create_configuration(template, sv, cnv, circosConfig, circosPlot)
    print(f'Output Circos Configuration file is written: {circosConfig}')
    return 0

if __name__=='__main__':
    main()

