import os
import subprocess
import re
import shlex

#line1_bed = '/home/users/cjyoon/Projects/myeloma/analysis/trafic/GRCh37_L1HS.strand.sorted.bed.gz'


with open('/home/users/cjyoon/Projects/myeloma/sample_info/tumor_match_abs_all.txt', 'r') as f:
    for line in f:
        sampleID, tumorBam, normalBam = line.strip().split()
        tumorID = re.sub(string=re.sub(string=tumorBam, pattern=r'/home/users/cjyoon/Projects/myeloma/bam/', repl=''), pattern=r'.sorted.md.indel.br.bam', repl='')
        normalID = re.sub(string=re.sub(string=normalBam, pattern=r'/home/users/cjyoon/Projects/myeloma/bam/', repl=''), pattern=r'.sorted.md.indel.br.bam', repl='')
        
        sv_circos = f'/home/users/cjyoon/Projects/myeloma/analysis/sv_intersect/{sampleID}_sv_filtered.circosData'
        cnv_circos = f'/home/users/cjyoon/Projects/myeloma/bam/{sampleID}.sorted.md.indel.br.bam.mpileup.100kbcov.absCN.circosData'
        snv_rainfall = f'/home/users/cjyoon/Projects/myeloma/analysis/snv_filter/{sampleID}_mns.vcf.rainfall.circosData'
        snv_vaf = f'/home/users/cjyoon/Projects/myeloma/analysis/snv_filter/{sampleID}_mns.vcf.vaf.circosData'

        
        cmd = f'python /home/users/cjyoon/scripts/circosPrep/circosConfigPrep.py -t /home/users/cjyoon/Projects/myeloma/analysis/circos/circos_full_template.conf -o . --sampleName {sampleID}_v2 --sv {sv_circos} --cnv {cnv_circos} --rainfall {snv_rainfall} --vaf {snv_vaf}'
        print(cmd)
        subprocess.call(shlex.split(cmd))


