#! /usr/bin/env python3
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
parser.add_argument('--vcf', help = 'vcf file from MuTect [optional]')
parser.add_argument('--stm', help = 'tsv file from Estimate Clonality')
parser.add_argument('-o', help='output directory')
parser.add_argument('-s', help='(y/n) or yes/no option to include sex chromosomes')
args=vars(parser.parse_args())

seg_file = args['seg']
vcf_file = args['vcf']
est_file = args['stm']
sex_chrom = args['s']
out_dir = args['o']

out_file=os.path.join(out_dir,'alltrunkfile.txt')

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
#since in sequqnza files are already sorted according to starting length no need to check the final condition where a bigger segment totally overlaps with the current segment
## need to write a function to return the copy number which the mutation falls under

def return_cn(s_arr,e_arr,hap_arr,chrom_arr,pos,chrom,hap,form):
    for i in range(len(s_arr)):
        if chrom_arr[i].split('\"') == chrom and hap_arr[i]==hap and s_arr[i]<=int(pos)<e_arr[i]:
            return i

def check_if_in(pos,s_arr,e_arr,lens,hap,h_arr,chrom,c_arr):
    for i in range(lens):
        if chrom == c_arr[i].strip('\"') and int(hap)==int(h_arr[i]) and (int(s_arr[i])<= int(pos)<int(e_arr[i])):
            return True,i

def main():
    class allele:
        def __init__(self):
            self.chrom = []
            self.start=[]
            self.end=[]
            self.form=[]
            self.hap=[]
            self.bearer=[]

        def add_entry(self,chrom_1,start_1,end_1):
            self.chrom.append(chrom_1)
            self.start.append(start_1)
            self.end.append(end_1)
            #hap.append(hap_1) dont know if i want to add this to the function or not

        def bearer_append(self,hap_1,bearer_1):
            self.hap.append(hap_1)
            self.bearer.append(bearer_1)

        def is_overlap(self,cur_pos,cur_hap,cur_chrom):
            for i in range(len(self.chrom)):
                if self.chrom[i].strip('\"') == cur_chrom and int(self.hap[i]) == int(cur_hap) and (int(self.start[i])<=cur_pos<int(self.end[i])):
                    return True, i

        def is_deletion(self,cur_pos,cur_hap,cur_chrom):
            condition,index = self.is_overlap(cur_pos,cur_hap,cur_chrom)
            if condition == True and self.form[index] == '-1':
                return True

        def append_index(self,chrm,start,end,form,hap,bearer):
            self.chrom.append(chrm)
            self.start.append(start)
            self.end.append(end)
            self.form.append(form)
            self.hap.append(hap)
            self.bearer.append(bearer)

    #initialising object instance
    seg = allele()

    if seg_file:
        with open(seg_file,'r') as segment:
            for line in segment:
                    rand_hap=np.random.randint(2)
                    elements=line.split('\t')
                    segment_length = len(seg.chrom)
                    if "chromosome" in elements[0]:
                            continue
                    elif int(elements[10]) > 6:
                        continue
                    else:
                            chrom_1=elements[0].strip('\"')
                            if sex_chrom.lower() == 'n':
                                    if chrom_1.lower() == 'x' or chrom_1.lower()=='y':
                                            continue
                            start_1=int(elements[1])
                            end_1=int(elements[2])
                            hap_seg = [int(elements[10]),int(elements[11])]
                            rand_hap=np.random.randint(2)
                            if hap_seg[1] > 1:
                                if rand_hap==1:
                                    seg.add_entry(chrom_1,start_1,end_1)
                                    seg.hap.append(0)
                                    seg.form.append(amplification(hap_seg[0]))
                                else:
                                    seg.add_entry(chrom_1,start_1,end_1)
                                    seg.hap.append(1)
                                    seg.form.append(amplification(hap_seg[1]))
                            elif hap_seg[0] < 1:
                                if rand_hap==1:
                                    seg.add_entry(chrom_1,start_1,end_1)
                                    seg.hap.append(0)
                                    seg.form.append(deletion())
                                else:
                                    seg.add_entry(chrom_1,start_1,end_1)
                                    seg.hap.append(1)
                                    seg.form.append(deletion())
                            else:
                                    for counter,i in enumerate(hap_seg):
                                            if i>1:
                                                    seg.add_entry(chrom_1,start_1,end_1)
                                                    if rand_hap==1:
                                                            seg.hap.append(0)
                                                    if rand_hap==0:
                                                            seg.hap.append(1)
                                                    seg.form.append(amplification(i))
                                            if i<1:
                                                    seg.add_entry(chrom_1,start_1,end_1)
                                                    seg.hap.append(rand_hap)
                                                    seg.form.append(deletion())
    segment.close()
    print('segments done',len(seg.chrom))

    #fn_open = gzip.open if vcf_file.endswith('.gz') else open
    if vcf_file:
        if vcf_file.endswith('.gz'):
            vcf_read = io.TextIOWrapper(io.BufferedReader(gzip.open(vcf_file,'rb')))
        else:
            vcf_read = open(vcf_file,'r')
        for line in vcf_read:
            rand_hap=np.random.randint(2)
            elements=line.strip().split()
            if line[0]!='#':
            #do not append to vcf file if variant overlaps with segment
                chrom_1=elements[0]
                if chrom_1.lower()=='mt':
                    continue
                #skips sex chromosomes based on input parameter
                if sex_chrom.lower() == 'n':
                    if chrom_1.lower() == 'x' or chrom_1.lower()=='y':
                        continue
                    if check_if_in(int(elements[1]),seg.start,seg.end,segment_length,rand_hap,seg.hap,chrom_1,seg.chrom):
                        continue
                    seg.chrom.append(elements[0])
                    seg.start.append(int(elements[1]))
                    seg.end.append(int(elements[1])+1)
                    seg.form.append(snv(elements[3],elements[4]))
                    seg.hap.append(rand_hap)
        vcf_read.close()

    stm = allele()
    if est_file:
        with open(est_file,'r') as tsv_read:
            for line in tsv_read:
                if 'patient' in line:
                    continue
                else:
                    elements = line.split()
                    ## 11 is the multiplier and determines early and late mutation
                    if float(elements[11]) > 1 and elements[21] == 'early':
                        sam,chm,pos,ref = elements[2].split(':')
                        alt = elements[4]
                        cn = [int(elements[9]),int(elements[10])] #[minor,major]
                        parental = [0,1]

                        for hap in parental:
                            stm.chrom.append(chm)
                            stm.start.append(pos)
                            stm.end.append(int(pos)+1)
                            stm.hap.append(np.random.randint(2))
                            stm.form.append(snv(ref,alt))
                            stm.bearer.append(range(cn[hap]))
        tsv_read.close()

    stm_filter = allele()
    i_index = -1
    for i in range(len(seg.chrom)):
        for j in range(len(stm.chrom)):
            if seg.chrom[i] == stm.chrom[j]:
                if int(seg.start[i]) <= int(stm.start[j]) < int(seg.end[i]):
                    if max(stm.bearer[j],default=99) <= int(seg.form[i]):
                        if i_index < i:
                            stm_filter.append_index(stm.chrom[j],stm.start[j],stm.end[j],stm.form[j],seg.hap[i],stm.bearer[j])
                            i_index = i
                            continue
                else:
                    if check_if_in(stm.start[j],stm_filter.start,stm_filter.end,len(stm_filter.chrom),stm.hap[j],stm_filter.hap,stm.chrom[j],stm_filter.chrom):
                        continue
                    if check_if_in(stm.start[j],seg.start,seg.end,len(seg.chrom),stm.hap[j],seg.hap,stm.chrom[j],seg.chrom):
                        continue
                    stm_filter.append_index(stm.chrom[j],stm.start[j],stm.end[j],stm.form[j],stm.hap[j],None)
                    break


    print('svn + cnv: ',len(seg.chrom))
    print('clonal variants: ',len(stm.chrom))
    print('filtered mutation: ',len(stm_filter.chrom))

    with open(out_file,'w') as out:
        '''uncomment to write all variants to out_file'''
        out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format('#chr','hap','start','end','var','bearer'))
        for i in range(len(seg.chrom)):
            if (int(seg.form[i])>6):
                continue
            out.write(str(seg.chrom[i].strip('\"'))+'\t'+str(seg.hap[i])+'\t'+str(seg.start[i])+'\t'+str(seg.end[i])+'\t'+str(seg.form[i])+'\n')

        for i in range(len(stm_filter.chrom)):
            if len(stm.bearer[i]) == 0 or stm.bearer[i] == None:
                out.write('{}\t{}\t{}\t{}\t{}\n'.format(stm_filter.chrom[i],stm_filter.hap[i],stm_filter.start[i],stm_filter.end[i],stm_filter.form[i]))
                continue
            if check_if_in(stm_filter.start[i],seg.start,seg.end,len(seg.start),stm_filter.hap[i],seg.hap,stm_filter.chrom[i],seg.chrom):
                out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(stm_filter.chrom[i],stm_filter.hap[i],stm_filter.start[i],stm_filter.end[i],stm_filter.form[i],','.join(repr(e) for e in stm.bearer[i])))
    out.close()


if __name__=="__main__":
    main()
