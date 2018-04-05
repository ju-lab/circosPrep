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
        manta_vcf = f'/home/users/cjyoon/Projects/myeloma/analysis/manta/{sampleID}_{tumorID}_{normalID}/results/variants/somaticSV.passonly.vcf.gz'
        delly_vcf = f'/home/users/cjyoon/Projects/myeloma/analysis/delly/{tumorID}_{normalID}_delly.somatic.passonly.vcf.gz'
#        if os.path.isfile(mutect_vcf) and os.path.isfile(manta_vcf) and os.path.isfile(delly_vcf):
        cmd = f'python /home/users/cjyoon/scripts/delly_vs_manta/delly_vs_manta.py -o . --manta {manta_vcf} --delly {delly_vcf} -s {sampleID}_sv_filtered -d 1000'
        print(cmd)
        subprocess.call(shlex.split(cmd))


