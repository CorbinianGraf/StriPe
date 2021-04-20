# StriPe
STructural vaRIation Primer dEsign

Usage:

StriPe.py --input path/to/structual_variation.vcf --sub1 Accession1 Accession2 --sub2 Accession3 Accession4 --out path/to/out.csv --ref path/to/reference.fasta --bDB path/to/blast_Database

Requirments:

This script uses Primer3, pyBedtools and Blast


StriPe.py [-h] --input INPUT --out OUT --ref REF --bDB BDB [--requirments REQUIRMENTS] [--dseq DSEQ] [--pfas PFAS]\
&nbsp;&nbsp;&nbsp;&nbsp;[--binput BINPUT] [--bout BOUT] --bDB BDB [--boutf BOUTF] [--p3out P3OUT] --sub1 SUB1 --sub2 SUB2\
&nbsp;&nbsp;&nbsp;&nbsp;[--pGC_min PGC_MIN] [--pGC_max PGC_MAX] [--pTM_min PTM_MIN] [--pTM_max PTM_MAX] [--pTM_opt PTM_OPT]\
&nbsp;&nbsp;&nbsp;&nbsp;[--opt_info OPT_INFO] [--f_fst F_FST] [--f_dels_min F_DELS_MIN]\

optional arguments:
  -h, --help                show this help message and exit\
  --requirments REQUIRMENTS list of packages needed to run\
  --input INPUT             path to input.vcf file\
  --out OUT                 desired output.csv file name and path\
  --ref REF                 path to reference.fasta genome\
  --dseq DSEQ               desired output.fasta file name and path for the sequences of the structual variances\
  --pfas PFAS               desired output.fasta file name and path for the sequences of the potential primers\
  --bout BOUT               desired output file name and path for the blast output\
  --bDB BDB                 path to the blast Database\
  --boutf BOUTF             desired blast outputformat (by default "10" i.e. csv)\
  --p3out P3OUT             optional: desired output.csv file name and path for the complet Primer3 output\
  --sub1 SUB1               Names / Accessions of subpopulation1 Individuals in VCF\
  --sub2 SUB2               Names / Accessions of subpopulation2 Individuals in VCF\
  --prod_max PROD_MAX       maximum size of PCR product of primers (default: 2000)\
  --p_regions P_REGIONS     size in bp of each region (left and right) next to a structual varaiance (default: 200)\
  --ps_min PS_MIN           minimum size for primers (defalt: 19)\
  --ps_max PS_MAX           maximum size for primers (default: 21)\
  --ps_opt PS_OPT           optimal size for primers (default: 20)\
  --pGC_min PGC_MIN         minimum GC content of primers (default: 20)\
  --pGC_max PGC_MAX         maximum GC content of primers (default: 80)\
  --pTM_min PTM_MIN         minimum annealing temperature of primers (default: 57)\
  --pTM_max PTM_MAX         maximum annealing temperature of primers (default: 60)\
  --pTM_opt PTM_OPT         optimal annealing temperature of primers (default: 63)\
  --opt_info OPT_INFO       optional information / columns to be carried over into the StriPe output file\
  --f_fst F_FST             minimum FST value for filtering (default: 1)\
  --f_dels_min F_DELS_MIN   minimum Structual Variation size to be allowed (default: 20)\
