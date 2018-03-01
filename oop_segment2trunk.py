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
parser.add_argument('-o', help='output directory')
parser.add_argument('-s', help='(y/n) or yes/no option to include sex chromosomes')
args=vars(parser.parse_args())

seg_file = args['seg']
vcf_file = args['vcf']
sex_chrom = args['s']
out_dir = args['o']

out_file=os.path.join(out_dir,'trunkfile.txt')

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

def amplification(count):
	return '+'+str(count-1)

def deletion():
	return '-1'

def check_if_in(position,s_arr,e_arr):
	for i in range(len(s_arr)):
		if int(s_arr[i])<=(position)<=int(e_arr[i]):
			return True

#chrom, start and end are lists
#_1 are values to be appended to the above lists
def main():
	class allele:
		def __init__(self):
			self.chrom = []
			self.start=[]
			self.end=[]
			self.form=[]
			self.hap=[]

		def add_entry(self,chrom_1,start_1,end_1):
			self.chrom.append(chrom_1)
			self.start.append(start_1)
			self.end.append(end_1)
			#hap.append(hap_1) dont know if i want to add this to the function or not

	a = allele()
	#A_al=A_allele(ploidy)

	with open(seg_file,'r') as segment:
		for line in segment:
			elements=line.split('\t')
			if "chromosome" in elements[0]:
				continue
			else:
				chrom_1=elements[0].strip('\"')
				if sex_chrom.lower() == 'n':
					if chrom_1.lower() == 'x' or chrom_1.lower()=='y':
						continue
				start_1=int(elements[1])-1
				end_1=int(elements[2])-1


				hap_seg = [int(elements[10]),int(elements[11])]
				for counter,i in enumerate(hap_seg):
					if i > 1:
						a.add_entry(chrom_1,start_1,end_1)
						a.hap.append(str(counter))
						a.form.append(amplification(i))
					if i<1:
						a.add_entry(chrom_1,start_1,end_1)
						a.hap.append(str(counter))
						a.form.append(deletion())
	segment.close()
	print('segments done')

	segment_length = len(a.chrom)
	with io.TextIOWrapper(io.BufferedReader(gzip.open(vcf_file,'rb'))) as vcf_read:
	#with open(vcf_file,'r') as vcf_read:
		for line in vcf_read:
			# print('file loaded')
			elements=line.strip().split()
			# print('elements splitted')
			if line[0]!='#':
			# print('first line skipped')
			# print(len(elements))
			# print(elements[1])
			# print(len(a.start))
			# print(len(a.end))
				if check_if_in(int(elements[1]),a.start,a.end)==True:
					continue
				chrom_1=elements[0]
				if sex_chrom.lower() == 'n':
					if chrom_1.lower() == 'x' or chrom_1.lower()=='y' or chrom_1.lower()=='mt':
						continue

				a.chrom.append(elements[0])
				a.start.append(int(elements[1])-1)
				a.end.append(int(elements[1]))
				a.form.append(snv(elements[3],elements[4]))
				a.hap.append(np.random.randint(2))
	vcf_read.close()




	with open(out_file,'w') as out:
		for i in range(len(a.chrom)):
			if i==0:
				out.write('#chr'+'\t' 'hap'+'\t' +'start'+'\t'+ 'end'+'\t'+ 'var'+'\n')
			else:
				out.write(str(a.chrom[i].strip('\"'))+'\t'+str(a.hap[i])+'\t'+str(a.start[i])+'\t'+str(a.end[i])+'\t'+str(a.form[i])+'\n')
	out.close()

main()
