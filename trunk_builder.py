import pandas as pd
import numpy as np
import yaml
import argparse

seg_file = "/home/abner/aquila/ITH/data/estclonality/hep198/sequenza/WHT260/output_segments.txt"
est_file = "/home/abner/aquila/ITH/data/estclonality/results/WHT260_earlylate.tsv"
data = "/home/abner/aquila/ITH/src/csite_scripts/config.yaml"

seg = pd.read_table(seg_file,header = 'infer')
est = pd.read_table(est_file, header = 'infer')

seg = seg.iloc[:,[0,1,2,10,11]]
est = est.iloc[:,[2,3,4,9,10,11,21]]

## finding mutations which are clonal in est
est = est[est['timing']=='early']

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
def opp_hap(hap):
    if hap == 1:
        return 0
    else:
        return 1

def cn_var(cn):
    if cn > 1:
        return '+'+str(cn-1)
    elif cn == 0:
        return '-1'

est['form'] = est.apply(lambda row: snv(row['Reference_Base'],row['Alternate_Base']), axis = 1)
est['sample'],est['chromosome'], est['start'], est['ref'] = zip(*est['mutation_id'].map(lambda x: x.split(':')))
est = est.drop(['sample','ref','Reference_Base','Alternate_Base','mut.multi','timing','mutation_id'], axis=1)
est['start'] = est['start'].astype(int)

#to find out which SNV are bearing segments
segment_end = seg['end.pos'].values
segment_start = seg['start.pos'].values
mutation_loci = est['start'].values
est_chr = est['chromosome'].values
seg_chr = seg['chromosome'].values
i, j = np.where((mutation_loci[:, None] >= segment_start) & (mutation_loci[:, None] < segment_end) & (est_chr[:,None]==seg_chr))

# non SNV bearing segments
seg_free = pd.DataFrame(seg.values[[index for index in range(len(seg)) if index not in j]], columns=seg.columns)
seg_free['major']=np.random.randint(2,size=len(seg_free))
seg_free['minor']=seg_free.major.apply(lambda x: opp_hap(x))
seg_free.columns
seg_free = seg_free.melt(id_vars=['chromosome','start.pos','end.pos','A','B'])
seg_free.head()

# SNV bearing segments
overlaps = pd.DataFrame(np.column_stack([est.values[i],seg.values[j,1:5]]),columns=np.append(est.columns.get_values(),seg.columns.get_values()[1:5]))
overlaps['end'] = overlaps.start.apply(lambda x: x+1)
overlaps['major']=np.random.randint(2,size=len(overlaps))
overlaps['minor']=overlaps.major.apply(lambda x: opp_hap(x))
overlaps = overlaps.melt(id_vars = overlaps.columns[0:10])

## To prepare DataFrame for output
'''
header for output file:
#chr hap start end var focal_cp
'''
## preparing segments for output
# preparing segments without bearer | filtering normal plodic segments
seg_free = seg_free[((seg_free['variable']=='major') & (seg_free['A'] != 1)) | ((seg_free['variable']=='minor') & (seg_free['B'] != 1))]
seg_free.columns
seg_free_major = seg_free[seg_free['variable']=='major'][['chromosome','value','start.pos','end.pos']]
seg_free_major['var'] = seg_free['A'].apply(lambda x : cn_var(x))
seg_free_minor = seg_free[seg_free['variable']=='minor'][['chromosome','value','start.pos','end.pos']]
seg_free_minor['var'] = seg_free['B'].apply(lambda x : cn_var(x))
out_seg_free = pd.concat([seg_free_major,seg_free_minor],ignore_index =True)
out_seg_free.columns = ['chromosome', 'value' ,'start','end','var']

# preparing mutations output (segments)
mutation_seg = overlaps[((overlaps['variable']=='major') & (overlaps['A'] != 1)) | ((overlaps['variable']=='minor') & (overlaps['B'] != 1))]
mutations_major_seg = overlaps[overlaps['variable']=='major'][['chromosome','value','start','end']]
mutations_major_seg['var'] = overlaps['A'].apply(lambda x: cn_var(x))
mutations_minor_seg = overlaps[overlaps['variable']=='major'][['chromosome','value','start','end']]
mutations_minor_seg['var'] = overlaps['B'].apply(lambda x : cn_var(x))
out_mutation_seg = pd.concat([mutations_major_seg,mutations_minor_seg], ignore_index = True)

# preparing mutation output (variants)
'''cherry picked all bearing mutations to have cn > 1 '''
mutations_major = overlaps[(overlaps['variable']=='major') & (overlaps['major_cn']>0)][['chromosome','value','start','end','form']]
mutations_major['focal_cp'] = overlaps[(overlaps['variable']=='major') & (overlaps['major_cn']>1)]['major_cn'].apply(lambda x: ','.join(repr(focal) for focal in range(x)))
mutations_minor = overlaps[(overlaps['variable']=='minor') & (overlaps['minor_cn']>0)][['chromosome','value','start','end','form']]
mutations_minor['focal_cp'] = overlaps[(overlaps['variable']=='minor') & (overlaps['minor_cn']>1)]['minor_cn'].apply(lambda x: ','.join(repr(focal) for focal in range(x)))
out_mutations = pd.concat([mutations_major,mutations_minor], ignore_index = True)
out_mutations.columns = ['chromosome', 'value' ,'start','end','var','focal_cp']

trunk_variants = pd.concat([out_seg_free,out_mutation_seg,out_mutations], ignore_index = True)[out_mutations.columns]
# filter out sex chromosomes
trunk_variants = trunk_variants[trunk_variants['chromosome']=='x' | trunk_variants['chromosome'] =='y']
trunk_variants[trunk_variants['end']-trunk_variants['start'] > 1].end.count()


config = yaml.load(open(data,'r'))
missing_haps = []
chrs = [str(i) for i in range(22)]
for chromosome in config['chromosomes']:
    for key, value in chromosome.items():
        chrs.remove(key)
        for nested_keys, nested_value in value.items():
            if nested_keys == "parental" and len(set(nested_value) < 2:
                    missing_haps.append([key,set(nested_value)])

for chromosomes in chr:
        trunk_variants = trunk_variants[trunk_variants['chromosome'] != chromosomes]
for values in missing_haps:
        trunk_variants = trunk_variants[(trunk_variants['chromosome'] != values[0]) & (trunk_variants['value'] != values[1])]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(describe='Script to generate trunk format for CSiTE')
    parser.add_argument('--stm', help = 'clonal mutations from estimate clonality.')
    parser.add_argument('--seg', help = 'CNV segments from sequenza')
    parser.add_arguments('--config', help = 'the same yaml format used in CSiTE.')
    args = parser.parse_args()

    ### Current status is that im 90% done. 
    ## 1) to do convert codes above to function and excute from here
    ## 2) set default == None for data.yaml in the event all chromosomes are simulated
    ## 3) perform sanity checks
