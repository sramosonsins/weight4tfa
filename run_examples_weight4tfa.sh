##weight4tfa v0.1beta (20170906) Sebastian E. Ramos-Onsins.
#
#Flags:
#      -i [path and name of the tfa file (text or gz)]
#      -o [path and name of the output weighted file (text or gz)]
#      -g [path of the GFF_file]
#         [add also: coding,noncoding,synonymous,nonsynonymous,silent, others (whatever annotated)]
#         [if 'synonymous', 'nonsynonymous', 'silent' add: Genetic_Code: Nuclear_Universal,mtDNA_Drosophila,mtDNA_Mammals,Other]
#         [if 'Other', introduce the single letter code for the 64 triplets in the order UUU UUC UUA UUG ... etc.]
#      -c [in case use coding regions, criteria to consider transcripts (max/min/first/long)]. DEFAULT: long
#      -n [name of each scaffold to analyze. tfa can be a list separated by commas(ex. -n chr1,chr2,chr3]
#      -l [total length of each scaffold(s). if more than one, separated by commas]
#   OPTIONAL PARAMETERS:
#      -h [help and exit]
#      -G [number of samples in the outgroup (if exist, the last samples)]. DEFAULT: 0
#      -m [masking regions: file indicating the start and the end of regions to be masked by 0 weights]. DEFAULT: NONE
#      -C [coordinates of regions: file indicating the start and the end of regions to be weighted (rest would be weighted as 0 if the file is included)]. DEFAULT: NONE

#To compile:
gcc ./sources/*.c -lm -o ./bin/weight4tfa -Wall -O3 -g -O0 -lz
cd ./Examples

../bin/weight4tfa -h > ../weight4tfa_help.txt

echo --------------------------------------------------------------------------------------------------
echo weight4tfa: Made for generating weighting file from GFF file
echo --------------------------------------------------------------------------------------------------
echo
echo weight4tfa.ex01
echo ../bin/weight4tfa -i 3Kallchr.tfa.gz -o 3Kallchr.tfa_weights.gz -g ./3Kallchr.gtf nonsynonymous Nuclear_Universal -c max -n chr10,chr12,chr14 -l 3000,2000,3000 -G 2
../bin/weight4tfa -i 3Kallchr.tfa.gz -o 3Kallchr.tfa_weights.gz -g ./3Kallchr.gtf nonsynonymous Nuclear_Universal -c max -n chr10,chr12,chr14 -l 3000,2000,3000 -G 2
echo
echo weight4ta.ex02
echo ../bin/weight4tfa -i 100Kallchr.tfa.gz -o 100Kallchr.tfa_weights.gz -g ./100Kallchr.gtf nonsynonymous Nuclear_Universal -c max -n chr10,chr12,chr14 -l 50000,100000,75000 -G  2
../bin/weight4tfa -i 100Kallchr.tfa.gz -o 100Kallchr.tfa_weights.gz -g ./100Kallchr.gtf nonsynonymous Nuclear_Universal -c max -n chr10,chr12,chr14 -l 50000,100000,75000 -G  2
echo
