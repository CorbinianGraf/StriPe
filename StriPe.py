#!/usr/bin/env python
# coding: utf-8



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
pd.options.mode.chained_assignment = None  # default='warn' , turns off dataframe overwrite warning
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

import os

######   default parameters

###   General

inputfile = '/Users/corbinian/Documents/projects/primer_desgin/scikit_allel/data/step5_germline_highcom.vcf'
outputfile = 'Primer3_out.csv'
referencefile = '/Users/corbinian/Documents/projects/primer_desgin/scikit_allel/data/Ahypochondriacus_459_v2.0.fa'
StriPe_error_file = '/Users/corbinian/Documents/projects/primer_desgin/scikit_allel/data/'
delsequencefile = 'struc_var_seqs.fa'
primer_fas = 'primer_desgin_primer_out.fasta'
blastDBfile = '/Users/corbinian/Documents/projects/primer_desgin/scikit_allel/data/Ahypochondriacus_459_v2.0'
blastoutformat = 10
blastoutfile = 'primer_desgin_blast_out.csv'
blastevalue = 0.1


subpopulation1 = ['PI642741', 'PI490518', 'PI608019'] #, 'PI490612']  # caudatus
subpopulation2 = ['PI511745', 'PI490466', 'PI652426'] #, 'PI669836']  # quitensis
#subpopulation2 = ['PI643058', 'PI576481', 'PI511717', 'PI511714'] # cruentus


###   Primer_info

# products need to be run on the same gel = cannot difffer in size to much
# for now: product size range = 100-1000bp;  > 1.4-1.6% Agarose
product_size_amax = 2000
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
referenceflag =  False

###    Filters

# additional columns/infos from the vcf file to keep (e.g. "Quality") next to the necessary ones (e.g. 'CHROM', 'POS', 'END', 'SVTYPE')
info_keep_opt = []

Fst_min = 1
del_size_min = 100
del_size_max = product_size_amax - 2*primer_regions



#####    parse arguments



parser = argparse.ArgumentParser()

parser.add_argument('--requirments', help='list of packages needed to run')
parser.add_argument('--input', required=True, help='path to input vcf file')
parser.add_argument('--out', required=True, help='desired output file and path')
parser.add_argument('--Serr', help='desired path for error file output, in case analysis can no longer be continued')
parser.add_argument('--ref', required=True, help='path to reference genome')
parser.add_argument('--dseq', help='desired output fasta file and path for the sequences of the structual variances')
parser.add_argument('--pfas', help='desired output fasta file and path for the sequences of the potential primers')
parser.add_argument('--sub1ref', help='raises reference as subpopulation1 flag as True. Choose one sample from --sub2 to be instead entered as --sub1 for this to work')
parser.add_argument('--bout')
parser.add_argument('--bDB', required=True)
parser.add_argument('--boutf', type=int)
parser.add_argument('--beval', type=float)
parser.add_argument('--sub1', nargs='+', required=True)   #### Textfile read
parser.add_argument('--sub2', nargs='+', required=True)
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
parser.add_argument('--opt_info', nargs='+')
parser.add_argument('--f_fst', type=int)
parser.add_argument('--f_dels_min', type=int)
parser.add_argument('--vcfe')

args = parser.parse_args()

# files and strings
if args.input:
    inputfile = args.input
if args.out:
    outputfile = args.out
if args.ref:
    referencefile = args.ref
if args.Serr:
    StriPe_error_file = args.Serr  
if args.dseq:
    delsequencefile = args.dseq
if args.pfas:
    primer_fas = args.pfas
if args.sub1ref:
    referenceflag = args.sub1ref
if args.bout:
    blastoutfile = args.bout
if args.bDB:
    blastDBfile = args.bDB
if args.boutf:
    blastoutformat = args.boutf
if args.beval:
    blastevalue = args.beval 
if args.vcfe:
    vcf_extract = args.vcfe
    
if args.sub1:
    subpopulation1 = args.sub1
if args.sub2:
    subpopulation2 = args.sub2

# Primer parameters   help = default values
if args.prod_max:
    product_size_amax = args.prod_max
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
del_size_max = product_size_amax - 2*primer_regions

blastinputfile = primer_fas

StriPe_error_file = str(StriPe_error_file + subpopulation1[0] + '.csv')

###   reading in populations

my_samples = subpopulation1 + subpopulation2

# later needed to reffere to each sample by position in genome array (calc Fst)
subpop1 = np.arange(len(subpopulation1))
subpop2 = np.arange(len(subpopulation1), len(my_samples))
#print(subpop1)
#print(subpop2)

###   data extraction

# needed for Genotype data
callset = allel.read_vcf(inputfile, fields='*', samples=my_samples, log=sys.stdout)
callset.keys()


##  what categories do exist in vcf file / what is their name
df = allel.vcf_to_dataframe(inputfile, fields='*')

# if X == True:  print(df)
if 'vcf_extract' in locals():
    df.to_csv(path_or_buf=vcf_extract, na_rep='NAN')


### reindex subpopulations to scikit-allels order

samples_df = pd.DataFrame(callset['samples'])
subpop1_reindex = []
subpop2_reindex = []
for i in range(len(my_samples)):
    if (i < len(subpopulation1)):
        subpop1_reindex.append(int(np.squeeze(samples_df.index[samples_df[0] == my_samples[i]].tolist())))
        
    if (i >= len(subpopulation1)):
        subpop2_reindex.append(int(np.squeeze(samples_df.index[samples_df[0] == my_samples[i]].tolist())))



#####   filtering

print('filtering vcf file \n')

###   filter out unimportant vcf info (dependent on what you want to do later)

# set info_keep to user specifications, when given as ARGV
info_keep = ['ID', 'CHROM', 'POS', 'END', 'SVTYPE', 'PRECISE']

info_keep = info_keep + info_keep_opt
    
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

# add reference genome to genotype data
if referenceflag == True:
    y = len(gt)
    x = 1
    z = 2
    R = [[[0 for k in range(z)] for j in range(x)] for i in range(y)]
    rt = allel.GenotypeArray(N, dtype='i1')
    rgt = np.append(nt,gt, axis=1)
    rgt = allel.GenotypeArray(ngt, dtype='i1')
    
    gt_del = rgt.subset(variant_selection)
else:
    gt_del = gt.subset(variant_selection)
# extract genotypes for all DEL variants

# find all homozygoth Genotypes (True / False)
homo = gt_del.is_hom()

# extract all variants that are homozygoth for all samples
homo = np.all(homo, axis=1)
gt_homo = gt_del.subset(homo)

#print("lengthomo ", len(gt_homo))

#print(np.count_nonzero(gt_homo), 'structual variances homozygous for all individuals')

# get all variant rows for gt_homo
variants_pass = variants_pass[homo]



###  find all variants at which the subpop/ sample Genotypes are segregating

if referenceflag == True:
    subpop1_reindex = [0]
    subpop2_reindex = list(range(1, len(my_samples)))

ac1 = gt_homo.count_alleles(subpop=subpop1_reindex)
ac2 = gt_homo.count_alleles(subpop=subpop2_reindex)
acu = allel.AlleleCountsArray(ac1 + ac2)

flt = acu.is_segregating()
gt_seg = gt_homo[flt]

print(np.count_nonzero(flt), 'positions are homoyzgous and segregating')

# get all variant data for segregating positions (flt=True)
variants_pass = variants_pass[flt]
variants_pass



###   calc hudson Fst
print('calculating Hudson Fst for each position')
ac1 = gt_seg.count_alleles(subpop=subpop1_reindex)
ac2 = gt_seg.count_alleles(subpop=subpop2_reindex)

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

lenge = list(range(0,len(gt_seg)))
gt_subpop1 = gt_seg[lenge,subpop1_reindex[0],1]
gt_subpop2 = gt_seg[lenge,subpop2_reindex[0],1]
gt_subpop2

finished_selection['GT_supop1'] = gt_subpop1
finished_selection['GT_supop2'] = gt_subpop2
finished_selection


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

if np.count_nonzero(finished_selection) == 0:
    print('No structual variation passed the PCR product restriction. List of considered SVs can be found in: ', StriPe_error_file)
    finished_selection.to_csv(path_or_buf=StriPe_error_file, na_rep='NAN')
    f = open(outputfile, "w")
    f.write("No structual variation passed the PCR product restriction. List of considered SVs can be found in: " + outputfile)
    f.close()
    
    exit()

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
del_pos = del_pos.sort_values('DEL_LEN')
del_pos = del_pos.drop(['SVTYPE', 'DEL_LEN'], axis=1)

if subpopulation1[0] == 'PI490518':
    print(finished_selection)

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

product_size_max = product_size_amax
for n in range(0, len(templates_array)):
    # select sequence
    if (n % 2) != 0:
        continue
    else:
        # define deletion size and product size
        if n == 0:
            deletion_size = finished_selection['DEL_LEN'].iloc[0]
            #product_size_max = deletion_size*10
            #if product_size_max >= product_size_amax:
                #product_size_max = product_size_amax
            product_size_min = deletion_size + primer_size_min
        else:
            deletion_size = finished_selection['DEL_LEN'].iloc[n//2]
            #product_size_max = deletion_size*10
            #if product_size_max >= product_size_amax:
               # product_size_max = product_size_amax
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

# control feature for primer3 output
primer3_out = str('/Users/corbinian/Documents/projects/primer_desgin/scikit_allel/data/test/' + subpopulation1[0] + '_primer3_output.csv')
#if subpopulation1[0] == 'PI576481':
    #primer3_res_df = primer3_res_df[primer3_res_df.PRIMER_LEFT_NUM_RETURNED != 0]

if len(primer3_res_df[primer3_res_df.PRIMER_LEFT_NUM_RETURNED == 0]) == len(primer3_res_df):
    print("no primers matching the Primer3 parameters could be found. Complet Primer3 output can be found in: ", primer3_out)
    primer3_res_df.to_csv(path_or_buf=primer3_out, na_rep='NAN')
    f = open(outputfile, "w")
    f.write("no primers matching the Primer3 parameters could be found. Complet Primer3 output can be found in: " + outputfile)
    f.close()
    
    exit()

##  get only relevant Primer info
#print(primer3_res_df.keys)
primer3_res_df = primer3_res_df[['PRIMER_LEFT_EXPLAIN', 'PRIMER_RIGHT_EXPLAIN', 'PRIMER_PAIR_EXPLAIN', 'PRIMER_LEFT_0_SEQUENCE', 'PRIMER_RIGHT_0_SEQUENCE', 'PRIMER_LEFT_0', 'PRIMER_RIGHT_0', 'PRIMER_PAIR_0_PRODUCT_SIZE', 'PRIMER_LEFT_1_SEQUENCE', 'PRIMER_RIGHT_1_SEQUENCE', 'PRIMER_LEFT_1', 'PRIMER_RIGHT_1', 'PRIMER_PAIR_1_PRODUCT_SIZE', 'PRIMER_LEFT_2_SEQUENCE', 'PRIMER_RIGHT_2_SEQUENCE', 'PRIMER_LEFT_2', 'PRIMER_RIGHT_2', 'PRIMER_PAIR_2_PRODUCT_SIZE', 'PRIMER_LEFT_3_SEQUENCE', 'PRIMER_RIGHT_3_SEQUENCE', 'PRIMER_LEFT_3', 'PRIMER_RIGHT_3', 'PRIMER_PAIR_3_PRODUCT_SIZE', 'PRIMER_LEFT_4_SEQUENCE', 'PRIMER_RIGHT_4_SEQUENCE', 'PRIMER_LEFT_4', 'PRIMER_RIGHT_4', 'PRIMER_PAIR_4_PRODUCT_SIZE']]

primer3_res_df = primer3_res_df[::2]


#print('Could not find any Primers (NAN) for', primer3_res_df['PRIMER_LEFT_0_SEQUENCE'].isnull().sum(), 'Positions')
#print('Could not find two Primers (NAN) for', primer3_res_df['PRIMER_LEFT_1_SEQUENCE'].isnull().sum(), 'Positions')
#print('Could not find three Primers (NAN) for', primer3_res_df['PRIMER_LEFT_2_SEQUENCE'].isnull().sum(), 'Positions')
#print('Could not find four Primers (NAN) for', primer3_res_df['PRIMER_LEFT_3_SEQUENCE'].isnull().sum(), 'Positions')
#print('Could not find five Primers (NAN) for', primer3_res_df['PRIMER_LEFT_4_SEQUENCE'].isnull().sum(), 'Positions')

primer3_res_df = primer3_res_df.reset_index(drop=True)

# primers to fasta
no_primers = primer3_res_df.PRIMER_LEFT_0_SEQUENCE.notnull()
primers_fasta = primer3_res_df[no_primers]

print(len(primer3_res_df), 'positions with Primers')

primers_fasta = primers_fasta.drop(['PRIMER_PAIR_EXPLAIN', 'PRIMER_PAIR_0_PRODUCT_SIZE', 'PRIMER_PAIR_1_PRODUCT_SIZE', 'PRIMER_PAIR_2_PRODUCT_SIZE', 'PRIMER_PAIR_3_PRODUCT_SIZE', 'PRIMER_PAIR_4_PRODUCT_SIZE', 'PRIMER_LEFT_EXPLAIN', 'PRIMER_RIGHT_EXPLAIN', 'PRIMER_LEFT_0', 'PRIMER_RIGHT_0', 'PRIMER_LEFT_1', 'PRIMER_RIGHT_1', 'PRIMER_LEFT_2', 'PRIMER_RIGHT_2', 'PRIMER_LEFT_3', 'PRIMER_RIGHT_3', 'PRIMER_LEFT_4', 'PRIMER_RIGHT_4'], axis=1)



# clear file if already filled
#print('writing primer sequences to ', primer_fas)
f = open(primer_fas, "w")
f.close()

# write in file
f = open(primer_fas, "a")

for t in range(0, len(primers_fasta.index)):
    f.write('>' + str(finished_selection['ID'].iloc[t]) + '_PRIMER_LEFT_0_SEQUENCE' + '\n')
    f.write(str(primers_fasta['PRIMER_LEFT_0_SEQUENCE'].iloc[t]) + '\n')
    f.write('>' + str(finished_selection['ID'].iloc[t]) + '_PRIMER_RIGHT_0_SEQUENCE' + '\n')
    f.write(str(primers_fasta['PRIMER_RIGHT_0_SEQUENCE'].iloc[t]) + '\n')
    
    f.write('>' + str(finished_selection['ID'].iloc[t]) + '_PRIMER_LEFT_1_SEQUENCE' + '\n')
    f.write(str(primers_fasta['PRIMER_LEFT_1_SEQUENCE'].iloc[t]) + '\n')
    f.write('>' + str(finished_selection['ID'].iloc[t]) + '_PRIMER_RIGHT_1_SEQUENCE' + '\n')
    f.write(str(primers_fasta['PRIMER_RIGHT_1_SEQUENCE'].iloc[t]) + '\n')
    
    f.write('>' + str(finished_selection['ID'].iloc[t]) + '_PRIMER_LEFT_2_SEQUENCE' + '\n')
    f.write(str(primers_fasta['PRIMER_LEFT_2_SEQUENCE'].iloc[t]) + '\n')
    f.write('>' + str(finished_selection['ID'].iloc[t]) + '_PRIMER_RIGHT_2_SEQUENCE' + '\n')
    f.write(str(primers_fasta['PRIMER_RIGHT_2_SEQUENCE'].iloc[t]) + '\n')
    
    f.write('>' + str(finished_selection['ID'].iloc[t]) + '_PRIMER_LEFT_3_SEQUENCE' + '\n')
    f.write(str(primers_fasta['PRIMER_LEFT_3_SEQUENCE'].iloc[t]) + '\n')
    f.write('>' + str(finished_selection['ID'].iloc[t]) + '_PRIMER_RIGHT_3_SEQUENCE' + '\n')
    f.write(str(primers_fasta['PRIMER_RIGHT_3_SEQUENCE'].iloc[t]) + '\n')
    
    f.write('>' + str(finished_selection['ID'].iloc[t]) + '_PRIMER_LEFT_4_SEQUENCE' + '\n')
    f.write(str(primers_fasta['PRIMER_LEFT_4_SEQUENCE'].iloc[t]) + '\n')
    f.write('>' + str(finished_selection['ID'].iloc[t]) + '_PRIMER_RIGHT_4_SEQUENCE' + '\n')
    f.write(str(primers_fasta['PRIMER_RIGHT_4_SEQUENCE'].iloc[t]) + '\n')

f.close()


print('running blast on primers')

#bbx(cmd='blastn', query=blastinputfile, db=blastDBfile, evalue=0.01, outfmt=blastoutformat, out=blastoutfile)

#blastx_cline = bbx(cmd='blastn', query=blastinputfile, db=blastDBfile, evalue=0.1, outfmt=blastoutformat, out=blastoutfile)
#blastx_cline()

#blastx_cline = bbx(cmd="blastn", query="primer_desgin_primer_out.fasta", db="Ahypochondriacus_459_v2.0", evalue=0.1, outfmt=10, out="primer_desgin_blast_out.csv")

blastevalue = str(blastevalue)
blastoutformat = str(blastoutformat)

#cmd = "blastn -task blastn -query primer_desgin_primer_out.fasta -db Ahypochondriacus_459_v2.0 -evalue 0.1 -outfmt 10 -out primer_desgin_blast_out.csv"
#os.system(cmd)

cmd = "blastn -task blastn -query " + blastinputfile + " -db " + blastDBfile + " -evalue " + blastevalue + " -outfmt " + blastoutformat + " -out " + blastoutfile
os.system(cmd)

#print("blastcommand: " + cmd)

print('writing summary')

primer3_res_df.index = finished_selection.index
#primer3_res_df
finished_selection = pd.concat([finished_selection, primer3_res_df], axis=1)



# filter out Deletions with no Primers
no_primers = finished_selection.PRIMER_LEFT_0_SEQUENCE.notnull()
finished_selection_primer3_res_df = finished_selection[no_primers]
finished_selection_primer3_res_df = finished_selection_primer3_res_df.drop(['PRIMER_LEFT_EXPLAIN', 'PRIMER_RIGHT_EXPLAIN', 'PRIMER_LEFT_0', 'PRIMER_RIGHT_0', 'PRIMER_LEFT_1', 'PRIMER_RIGHT_1', 'PRIMER_LEFT_2', 'PRIMER_RIGHT_2', 'PRIMER_LEFT_3', 'PRIMER_RIGHT_3', 'PRIMER_LEFT_4', 'PRIMER_RIGHT_4'], axis=1)


### use primer_fas headers to count number of times each primer alignes -> add to output

with open(primer_fas, 'r') as file:
    fastalines = file.read()
    fastahead = fastalines.split('\n')[0::2][:-1]

with open(blastoutfile, 'r') as file2:
    blastres = file2.read()

pcount = []
for i in range(len(fastahead)):
    primername = fastahead[i][1:]
    pcount.append(blastres.count(primername))

primer0l = pcount[0::10]
primer0r = pcount[1::10]
primer1l = pcount[2::10]
primer1r = pcount[3::10]
primer2l = pcount[4::10]
primer2r = pcount[5::10]
primer3l = pcount[6::10]
primer3r = pcount[7::10]
primer4l = pcount[8::10]
primer4r = pcount[9::10]

#col_num = len(finished_selection_primer3_res_df.columns)

finished_selection_primer3_res_df.insert(loc=23, column='PRIMER_RIGHT_4_blastcount', value=primer4r)
finished_selection_primer3_res_df.insert(loc=22, column='PRIMER_LEFT_4_blastcount', value=primer4l)
finished_selection_primer3_res_df.insert(loc=20, column='PRIMER_RIGHT_3_blastcount', value=primer3r)
finished_selection_primer3_res_df.insert(loc=19, column='PRIMER_LEFT_3_blastcount', value=primer3l)
finished_selection_primer3_res_df.insert(loc=17, column='PRIMER_RIGHT_2_blastcount', value=primer2r)
finished_selection_primer3_res_df.insert(loc=16, column='PRIMER_LEFT_2_blastcount', value=primer2l)
finished_selection_primer3_res_df.insert(loc=14, column='PRIMER_RIGHT_1_blastcount', value=primer1r)
finished_selection_primer3_res_df.insert(loc=13, column='PRIMER_LEFT_1_blastcount', value=primer1l)
finished_selection_primer3_res_df.insert(loc=11, column='PRIMER_RIGHT_0_blastcount', value=primer0r)
finished_selection_primer3_res_df.insert(loc=10, column='PRIMER_LEFT_0_blastcount', value=primer0l)


## filter deletion dataframe for succesful Primers


## save Primer list to csv

finished_selection_primer3_res_df.to_csv(path_or_buf=outputfile, na_rep='NAN')

