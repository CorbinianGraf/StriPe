#!/usr/bin/env python
# coding: utf-8

# In[31]:


###
import argparse

import numpy as np
import h5py
import matplotlib.pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
import bcolz
import pandas as pd
import allel
from collections import Counter
import sys
import getopt

### 
from pybedtools import BedTool
import pybedtools

### 
import primer3

### 
from pydna.dseqrecord import Dseqrecord
from pydna.readers import read
#from pydna.amplify import pcr
from pydna.amplify import pcr
from pydna.primer import Primer
from Bio.SeqRecord import SeqRecord
#np.seterr(divide='ignore', invalid='ignore')

import Bio.Blast.NCBIWWW as bb
from Bio.Blast.Applications import NcbiblastxCommandline as bbx



######   default parameters

###   General

inputfile = '/Users/corbinian/Documents/projects/primer_desgin/scikit_allel/data/step5_germline_highcom.vcf'
outputfile = '/Users/corbinian/Documents/projects/primer_desgin/scikit_allel/data/Primer3_output.csv'
referencefile = '/Users/corbinian/Documents/projects/primer_desgin/scikit_allel/data/Ahypochondriacus_459_v2.0.fa'
delsequencefile = '/Users/corbinian/Documents/projects/primer_desgin/scikit_allel/data/Ahyp_ref_seqs.fa'
primer_fas = '/Users/corbinian/Documents/projects/primer_desgin/scikit_allel/data/primer_desgin_primer_out.fasta'
blastinputfile = '/Users/corbinian/Documents/projects/primer_desgin/scikit_allel/data/primer_desgin_primer_out.fasta'
blastDBfile = '/Users/corbinian/Documents/projects/primer_desgin/scikit_allel/data/Ahypochondriacus_459_v2.0'
blastoutformat = 10
blastoutfile = '/Users/corbinian/Documents/projects/primer_desgin/scikit_allel/data/primer_desgin_blast_out.csv'
primer3out = '/Users/corbinian/Documents/projects/primer_desgin/scikit_allel/data/primer3_out.csv'

# how to read multiple separate elements for same arg? -spop1 PI1 PI2 PI3 => ['PI1', 'PI2' ...]
subpopulation1 = ['PI642741', 'PI490518', 'PI608019'] #, 'PI490612']  # caudatus
subpopulation2 = ['PI511745', 'PI490466', 'PI652426'] #, 'PI669836']  # quitensis
#subpopulation2 = ['PI643058', 'PI576481', 'PI511717', 'PI511714'] # cruentus


###   Primer_info

# products need to be run on the same gel = cannot difffer in size to much
# for now: product size range = 100-1000bp;  > 1.4-1.6% Agarose
product_size_max = 2000
primer_regions = 200

primer_size_min = 19
primer_size_opt = 20
primer_size_max = 21
primer_GC_min = 20
primer_GC_max = 80
primer_TM_min = 57
primer_TM_opt = 60
primer_TM_max = 63
# product_size_min = define later by each individuals actual deletion size


###    Filters

# additional columns/infos from the vcf file to keep (e.g. "Quality") next to the necessary ones (e.g. 'CHROM', 'POS', 'END', 'SVTYPE')
info_keep_opt = []

Fst_min = 1
del_size_min = 20
del_size_max = product_size_max - 2*primer_regions



#####    parse arguments


parser = argparse.ArgumentParser()

parser.add_argument('--requirments', help='list of packages needed to run')
parser.add_argument('--input', required=True, help='path to input vcf file')
parser.add_argument('--out', required=True, help='desired output file and path')
parser.add_argument('--ref', required=True, help='path to reference genome')
parser.add_argument('--dseq', help='desired output fasta file and path for the sequences of the structual variances')
parser.add_argument('--pfas', help='desired output fasta file and path for the sequences of the potential primers')
parser.add_argument('--binput')
parser.add_argument('--bout')
parser.add_argument('--bDB', required=True)
parser.add_argument('--boutf')
parser.add_argument('--p3out')
parser.add_argument('--sub1', required=True)   #### Textfile read
parser.add_argument('--sub2', required=True)
parser.add_argument('--prod_max', default=2000, type=int)
parser.add_argument('--p_regions', type=int)
parser.add_argument('--ps_min', type=int)
parser.add_argument('--ps_max', type=int)
parser.add_argument('--ps_opt', type=int)
parser.add_argument('--pGC_min', type=int)
parser.add_argument('--pGC_max', type=int)
parser.add_argument('--pTM_min', type=float)
parser.add_argument('--pTM_max', type=float)
parser.add_argument('--pTM_opt', type=float)
parser.add_argument('--opt_info')
parser.add_argument('--f_fst', type=int)
parser.add_argument('--f_dels_min', type=int)

args = parser.parse_args()

# files and strings
if args.input:
    inputfile = args.input
if args.out:
    outputfile = args.out
if args.ref:
    referencefile = args.ref
if args.dseq:
    delsequencefile = args.dseq
if args.pfas:
    primer_fas = args.pfas
if args.binput:
    blastinputfile = args.binput
if args.bout:
    blastoutfile = args.bout
if args.bDB:
    blastDBfile = args.bDB
if args.boutf:
    blastoutformat = args.boutf
if args.p3out:
    primer3out = args.p3out
    
if args.sub1:
    subpopulation1 = args.sub1
if args.sub2:
    subpopulation2 = args.sub2

# Primer parameters   help = default values
if args.prod_max:
    product_size_max = args.prod_max
if args.p_regions:
    primer_regions = args.p_regions
if args.ps_min:
    primer_size_min = args.ps_min
if args.ps_max:
    primer_size_max = args.ps_max
if args.ps_opt:
    primer_size_opt = args.ps_opt
if args.pGC_min:
    primer_GC_min = args.pGC_min
if args.pGC_max:
    primer_GC_max = args.pGC_max
if args.pTM_min:
    primer_TM_min = args.pTM_min
if args.pTM_max:
    primer_TM_opt = args.pTM_max
if args.pTM_opt:
    primer_TM_max = args.pTM_opt

if args.opt_info:
    info_keep_opt = args.opt_info

# filters
if args.f_fst:
    Fst_min = args.f_fst
if args.f_dels_min:
    del_size_min = args.f_dels_min
del_size_max = product_size_max - 2*primer_regions



###   reading in populations

my_samples = subpopulation1 + subpopulation2

# later needed to reffere to each sample by position in genome array (calc Fst)
subpop1 = np.arange(len(subpopulation1))
subpop2 = np.arange(len(subpopulation1), len(my_samples))



###   data extraction

# needed for Genotype data
callset = allel.read_vcf(inputfile, fields='*', samples=my_samples, log=sys.stdout)
callset.keys()


##  what categories do exist in vcf file / what is their name
df = allel.vcf_to_dataframe(inputfile, fields='*')

# if X == True:  print(df)


#####   filtering

print('filtering vcf file \n')

###   filter out unimportant vcf info (dependent on what you want to do later)

# set info_keep to user specifications, when given as ARGV
info_keep = ['ID', 'CHROM', 'POS', 'END', 'SVTYPE', 'QUAL', 'PRECISE', 'REF']

info_keep = info_keep + info_keep_opt

if 'info_filter' in vars():
    info_keep = info_filter
    
df = df.filter(info_keep, axis=1)


###    filter for Deletions only

# define data entry to be filtered by
del_true = '(SVTYPE == "DEL") | (SVTYPE == "DUP")'

# determine rows with same entry (True / False)
variant_selection = df.eval(del_true)[:]

print(np.count_nonzero(df.eval('SVTYPE == "DEL"')[:]), 'deletions')
print(np.count_nonzero(df.eval('SVTYPE == "DUP"')[:]), 'duplications')
print(np.count_nonzero(variant_selection), 'total structual variances')

# extract all rows with DEL entries (variant_selection=True) 
variants_pass = df[variant_selection]
variants_pass



###   genotype filtering

# get all genotypes
gt = allel.GenotypeArray(callset['calldata/GT'])

# extract genotypes for all DEL variants
gt_del = gt.subset(variant_selection)

# find all homozygoth Genotypes (True / False)
homo = gt_del.is_hom()

# extract all variants that are homozygoth for all samples
homo = np.all(homo, axis=1)
gt_homo = gt_del.subset(homo)

#print(np.count_nonzero(gt_homo), 'structual variances homozygous for all individuals')

# get all variant rows for gt_homo
variants_pass = variants_pass[homo]



###  find all variants at which the subpop/ sample Genotypes are segregating

ac1 = gt_homo.count_alleles(subpop=subpop1)
ac2 = gt_homo.count_alleles(subpop=subpop2)
acu = allel.AlleleCountsArray(ac1 + ac2)


flt = acu.is_segregating()
gt_seg = gt_homo[flt]

print(np.count_nonzero(flt), 'positions are homoyzgous and segregating')

# get all variant data for segregating positions (flt=True)
variants_pass = variants_pass[flt]
variants_pass



###   calc hudson Fst
print('calculating Hudson Fst for each position')
ac1 = gt_seg.count_alleles(subpop=subpop1)
ac2 = gt_seg.count_alleles(subpop=subpop2)

num, den = allel.hudson_fst(ac1, ac2)
num_fix = np.nan_to_num(num)
den_fix = np.nan_to_num(den)

fst = num_fix / den_fix



###    append Fst to variants

# Fst array to dataframe
fst_df = pd.DataFrame(fst, columns=['Fst_hudson'])

# set index of Fst to index of variants_pass
fst_df.index = variants_pass.index

# join
finished_selection = variants_pass.join(fst_df)


##### add Genotypes for each individual (PI6456547)



###  Filter for Fst

#Fst_min = 1
finished_selection = finished_selection[finished_selection['Fst_hudson'] >= Fst_min]

print(np.count_nonzero(finished_selection), 'positions are segregating between populations')



###    filter for deletion length

# determine and add deletion length for each variant
finished_selection['DEL_LEN'] = finished_selection['END'] - finished_selection['POS']

# sort by deletion size
finished_selection = finished_selection.sort_values('DEL_LEN')
finished_selection


##  apply deletion size filters

finished_selection = finished_selection[(finished_selection.SVTYPE == 'DEL') & (finished_selection.DEL_LEN < del_size_max) | (finished_selection.SVTYPE == 'DUP') & (finished_selection.DEL_LEN < (del_size_max//2))]

#invalid_del_size = finished_selection.DEL_LEN < del_size_max
#finished_selection = finished_selection[invalid_del_size]

invalid_del_size = finished_selection.DEL_LEN > del_size_min
finished_selection = finished_selection[invalid_del_size]

print(np.count_nonzero(finished_selection), 'positions passed PCR product restrictions')



#####    Primer search

print('starting Primer search')

###   get deletion position

del_pos = finished_selection[['CHROM', 'POS', 'END', 'SVTYPE', 'DEL_LEN']]

del_pos['POS'] = del_pos.POS - primer_regions
del_pos['END'] = del_pos.END + primer_regions

del_pos_del = del_pos[del_pos.SVTYPE == 'DEL']
del_pos_dup = del_pos[del_pos.SVTYPE == 'DUP']

del_pos_dup['POS'] = del_pos_dup.POS - del_pos_dup.DEL_LEN
del_pos_dup['END'] = del_pos_dup.END - del_pos_dup.DEL_LEN

del_pos = pd.concat([del_pos_del, del_pos_dup])
del_pos = del_pos.drop(['SVTYPE', 'DEL_LEN'], axis=1)



###   get fasta sequence of deletion from reference genome

bedtools_input = del_pos.to_string(index=False, header=False)

bedtool_element = pybedtools.BedTool(bedtools_input, from_string=True)
fasta = pybedtools.example_filename(referencefile)
fastaout = delsequencefile

# save sequence to fasta file
bedtool_element = bedtool_element.sequence(fi=fasta, fo=fastaout)



##  extract sequences from fasta file

try:
    templates = open(delsequencefile)
    templates_array = templates.readlines()
    
finally:
    templates.close()
    
templates_array = [i.replace("\n", "") for i in templates_array]



####    calaculate Primers

print('calculating Primers')

primer3_res_df = pd.DataFrame()

for n in range(0, len(templates_array)):
    # select sequence
    if (n % 2) != 0:
        continue
    else:
        # define deletion size and product size
        if n == 0:
            deletion_size = finished_selection['DEL_LEN'].iloc[0]
            product_size_min = deletion_size + primer_size_min
        else:
            deletion_size = finished_selection['DEL_LEN'].iloc[n//2]
            product_size_min = deletion_size + primer_size_min
        #print(n//2, deletion_size)
        n2 = n+1
        primer3_res = primer3.bindings.designPrimers(
        {
            'SEQUENCE_ID': templates_array[n],
            'SEQUENCE_TEMPLATE': templates_array[n2],
            'SEQUENCE_EXCLUDED_REGION': [int(primer_regions), int(deletion_size)]
        },
        {
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_PICK_LEFT_PRIMER': 1,
            'PRIMER_PICK_INTERNAL_OLIGO': 0,
            'PRIMER_PICK_RIGHT_PRIMER': 1,
            'PRIMER_MIN_SIZE': 19,
            'PRIMER_MAX_SIZE': 21,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_MAX_POLY_X': 100,
            'PRIMER_INTERNAL_MAX_POLY_X': 100,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 12,
            'PRIMER_MAX_SELF_END': 8,
            'PRIMER_PAIR_MAX_COMPL_ANY': 12,
            'PRIMER_PAIR_MAX_COMPL_END': 8,
            'PRIMER_PRODUCT_SIZE_RANGE': [[product_size_min, product_size_max]],
        })
        primer3_res_df = primer3_res_df.append(pd.DataFrame(primer3_res, index=[n, n2]))



##  get only relevant Primer info

primer3_res_df = primer3_res_df[['PRIMER_LEFT_EXPLAIN', 'PRIMER_RIGHT_EXPLAIN', 'PRIMER_PAIR_EXPLAIN', 'PRIMER_LEFT_0_SEQUENCE', 'PRIMER_RIGHT_0_SEQUENCE', 'PRIMER_LEFT_0', 'PRIMER_RIGHT_0', 'PRIMER_PAIR_0_PRODUCT_SIZE', 'PRIMER_LEFT_1_SEQUENCE', 'PRIMER_RIGHT_1_SEQUENCE', 'PRIMER_LEFT_1', 'PRIMER_RIGHT_1', 'PRIMER_PAIR_1_PRODUCT_SIZE', 'PRIMER_LEFT_2_SEQUENCE', 'PRIMER_RIGHT_2_SEQUENCE', 'PRIMER_LEFT_2', 'PRIMER_RIGHT_2', 'PRIMER_PAIR_2_PRODUCT_SIZE', 'PRIMER_LEFT_3_SEQUENCE', 'PRIMER_RIGHT_3_SEQUENCE', 'PRIMER_LEFT_3', 'PRIMER_RIGHT_3', 'PRIMER_PAIR_3_PRODUCT_SIZE', 'PRIMER_LEFT_4_SEQUENCE', 'PRIMER_RIGHT_4_SEQUENCE', 'PRIMER_LEFT_4', 'PRIMER_RIGHT_4', 'PRIMER_PAIR_4_PRODUCT_SIZE']]

primer3_res_df = primer3_res_df[::2]


print('Could not find any Primers (NAN) for', primer3_res_df['PRIMER_LEFT_0_SEQUENCE'].isnull().sum(), 'Positions')
print('Could not find two Primers (NAN) for', primer3_res_df['PRIMER_LEFT_1_SEQUENCE'].isnull().sum(), 'Positions')
print('Could not find three Primers (NAN) for', primer3_res_df['PRIMER_LEFT_2_SEQUENCE'].isnull().sum(), 'Positions')
print('Could not find four Primers (NAN) for', primer3_res_df['PRIMER_LEFT_3_SEQUENCE'].isnull().sum(), 'Positions')
print('Could not find five Primers (NAN) for', primer3_res_df['PRIMER_LEFT_4_SEQUENCE'].isnull().sum(), 'Positions')

primer3_res_df = primer3_res_df.reset_index(drop=True)

# primers to fasta
no_primers = primer3_res_df.PRIMER_LEFT_0_SEQUENCE.notnull()
primers_fasta = primer3_res_df[no_primers]

print(len(primer3_res_df), 'positions with Primers')

primers_fasta = primers_fasta.drop(['PRIMER_PAIR_EXPLAIN', 'PRIMER_PAIR_0_PRODUCT_SIZE', 'PRIMER_PAIR_1_PRODUCT_SIZE', 'PRIMER_PAIR_2_PRODUCT_SIZE', 'PRIMER_PAIR_3_PRODUCT_SIZE', 'PRIMER_PAIR_4_PRODUCT_SIZE', 'PRIMER_LEFT_EXPLAIN', 'PRIMER_RIGHT_EXPLAIN', 'PRIMER_LEFT_0', 'PRIMER_RIGHT_0', 'PRIMER_LEFT_1', 'PRIMER_RIGHT_1', 'PRIMER_LEFT_2', 'PRIMER_RIGHT_2', 'PRIMER_LEFT_3', 'PRIMER_RIGHT_3', 'PRIMER_LEFT_4', 'PRIMER_RIGHT_4'], axis=1)



# clear file if already filled
print('writing primer sequences to ', primer_fas)
f = open(primer_fas, "w")
f.close()

# write in file
f = open(primer_fas, "a")

for t in range(0, len(primers_fasta.index)):
    f.write('>' + str(finished_selection['ID'].iloc[t]) + '_PRIMER_LEFT_0_SEQUENCE' + '\n')
    f.write(primers_fasta['PRIMER_LEFT_0_SEQUENCE'].iloc[t] + '\n')
    f.write('>' + str(finished_selection['ID'].iloc[t]) + '_PRIMER_RIGHT_0_SEQUENCE' + '\n')
    f.write(primers_fasta['PRIMER_RIGHT_0_SEQUENCE'].iloc[t] + '\n')
    
    f.write('>' + str(finished_selection['ID'].iloc[t]) + '_PRIMER_LEFT_1_SEQUENCE' + '\n')
    f.write(primers_fasta['PRIMER_LEFT_1_SEQUENCE'].iloc[t] + '\n')
    f.write('>' + str(finished_selection['ID'].iloc[t]) + '_PRIMER_RIGHT_1_SEQUENCE' + '\n')
    f.write(primers_fasta['PRIMER_RIGHT_1_SEQUENCE'].iloc[t] + '\n')
    
    f.write('>' + str(finished_selection['ID'].iloc[t]) + '_PRIMER_LEFT_2_SEQUENCE' + '\n')
    f.write(primers_fasta['PRIMER_LEFT_2_SEQUENCE'].iloc[t] + '\n')
    f.write('>' + str(finished_selection['ID'].iloc[t]) + '_PRIMER_RIGHT_2_SEQUENCE' + '\n')
    f.write(primers_fasta['PRIMER_RIGHT_2_SEQUENCE'].iloc[t] + '\n')
    
    f.write('>' + str(finished_selection['ID'].iloc[t]) + '_PRIMER_LEFT_3_SEQUENCE' + '\n')
    f.write(primers_fasta['PRIMER_LEFT_3_SEQUENCE'].iloc[t] + '\n')
    f.write('>' + str(finished_selection['ID'].iloc[t]) + '_PRIMER_RIGHT_3_SEQUENCE' + '\n')
    f.write(primers_fasta['PRIMER_RIGHT_3_SEQUENCE'].iloc[t] + '\n')
    
    f.write('>' + str(finished_selection['ID'].iloc[t]) + '_PRIMER_LEFT_4_SEQUENCE' + '\n')
    f.write(primers_fasta['PRIMER_LEFT_4_SEQUENCE'].iloc[t] + '\n')
    f.write('>' + str(finished_selection['ID'].iloc[t]) + '_PRIMER_RIGHT_4_SEQUENCE' + '\n')
    f.write(primers_fasta['PRIMER_RIGHT_4_SEQUENCE'].iloc[t] + '\n')

f.close()


print('running blast on primers')

#bbx(cmd='blastn', query=blastinputfile, db=blastDBfile, evalue=0.001, outfmt=blastoutformat, out=blastoutfile)

blastx_cline = bbx(cmd='blastn', query=blastinputfile, db=blastDBfile, evalue=0.01, outfmt=blastoutformat, out=blastoutfile)
blastx_cline()



#get_ipython().run_cell_magic('bash', '', 'blastn -task "blastn" -query "/Users/corbinian/Desktop/test_primers.fasta" -db "/Users/corbinian/Documents/projects/primer_desgin/scikit_allel/data/Ahypochondriacus_459_v2.0_DB" -evalue 0.001 -outfmt 10 -out "/Users/corbinian/Desktop/test_primers_blast.csv"')

#get_ipython().run_cell_magic('bash', '', 'blastn -task "blastn" -query blastinputfile -db blastDBfile -evalue 0.001 -outfmt blastoutformat -out blastoutfile')



primer3_res_df.index = finished_selection.index
#primer3_res_df
finished_selection = pd.concat([finished_selection, primer3_res_df], axis=1)



# filter out Deletions with no Primers
no_primers = finished_selection.PRIMER_LEFT_0_SEQUENCE.notnull()
finished_selection_primer3_res_df = finished_selection[no_primers]
finished_selection_primer3_res_df = finished_selection_primer3_res_df.drop(['PRIMER_LEFT_EXPLAIN', 'PRIMER_RIGHT_EXPLAIN', 'PRIMER_LEFT_0', 'PRIMER_RIGHT_0', 'PRIMER_LEFT_1', 'PRIMER_RIGHT_1', 'PRIMER_LEFT_2', 'PRIMER_RIGHT_2', 'PRIMER_LEFT_3', 'PRIMER_RIGHT_3', 'PRIMER_LEFT_4', 'PRIMER_RIGHT_4'], axis=1)


## filter deletion dataframe for succesful Primers


## save Primer list to csv

finished_selection_primer3_res_df.to_csv(path_or_buf=primer3out, na_rep='NAN')

