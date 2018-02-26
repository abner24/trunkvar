# ##################################################
# Title : segment2trunk.py
# Author : Lim Ming Sun (Abner)
# Date : Feb 23 2018 09:33:41
# Description : A script to generate the trunk file 
# ##################################################

import os
import io
import gzip
import numpy as np
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--seg', help = 'sequenza segments file')
parser.add_argument('--vcf', help = 'vcf file from MuTect')
parser.add_argument('-p', help='Mean ploidy value obtained from sequenza')
args=vars(parser.parse_args())

seg_file = args['seg']
vcf_file = args['vcf']
ploidy = int(args['p'])

#if  sys.argv[2] != None:
#    seg_file=sys.argv[2]
#else:
#    print('continuing without a segment_output text file from sequenza')
#if sys.argv[3] !=None:
#    ploidy=int(sys.argv[3])
#else:
#    print('please input a ploidy number')

out_dir=os.path.dirname(seg_file)
out_file=os.path.join(out_dir,'output.trunk_var')

def A_allele(ploidy):
    if ploidy%2==0:
        x=ploidy/2
    else:
        x=np.floor(ploidy/2)+1
    return x

def snv(ref,alt):
	#returns the form of the allele
	mutation_matrix={'N':['N','N','N'],
				  'A':['G','C','T'],
				  'G':['A','C','T'],
				  'C':['T','A','G'],
				  'T':['C','A','G']}
	for form,values in enumerate(mutation_matrix[ref]):
		if values==alt:
			return form

#chrom, start and end are lists
#_1 are values to be appended to the above lists
def add_entry(chrom,start,end,chrom_1,start_1,end_1):
	global chrom start end
	chrom.append(chrom_1)
	start.append(start_1)
	end.append(end_1)
				
def main():
    chrom=[]
    start=[]
    end=[]
    form=[]
    hap=[]
    A_al=A_allele(ploidy)

    with io.TextIOWrapper(io.BufferedReader(gzip.open(vcf_file,'rb'))) as vcf_read:
        for line in vcf_read:
            elements=line.strip().split()
            if line[0]!='#':
                chrom.append(elements[0])
                start.append(int(elements[1])-1)
                end.append(int(elements[1]))
			 form.append(snv(elements[3],elements[4]))
                if (elements[10].split(':')[0])=="0/1":
                    hap.append('0')
                else:
                    hap.append(elements[10].split(':')[0])
    vcf_read.close()


    with open(seg_file,'r') as segment:
        for line in segment:
            elements=line.split('\t')
            if "chromosome" in elements[0]:
            	continue
            else:
				chrom_1=elements[0]
				start_1=int(elements[1]-1)
				end_1=int(elements[2]-1)

                #chrom.append(elements[0])
                #start.append(int(elements[1])-1)
                #end.append(int(elements[2])-1)
                #10 is A 11 is B
                hap_0 = int(elements[10])
                hap_1 = int(elements[11])
                if hap0+hap1==2:
				    if (hap0&&hap1)==1:
						continue
				    elif hap0<1:
						add_entry(chrom,start,end,chrom_1,start_1,end_1)
								
                    form.append("0")
                    hap.append('0')
                elif (hap_0+hap_1) > ploidy:
                    state.append('+'+str((hap_0+hap_1)-ploidy))
                    if hap_0>hap_1:
                        hap.append('0')
                    else:
                        hap.append('1')
                else:
                    form.append("-1")
                    if hap_1>hap_0:
                        hap.append('0')
                    else:
                        hap.append('1')
    segment.close()

    with open(out_file,'w') as out:
        for i in range(len(chrom)):
                if i==0:
                    out.write('#chr'+'\t' 'hap'+'\t' +'start'+'\t'+ 'end'+'\t'+ 'var'+'\n')
                else:
                    out.write(chrom[i].strip('\"')+'\t'+hap[i]+'\t'+str(start[i])+'\t'+str(end[i])+'\t'+form[i]+'\n')
    out.close()

main()
