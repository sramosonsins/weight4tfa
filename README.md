# weight4tfa
Calculates weights for tfasta files (must be .gz and indexed)

#weight4tfa v0.1beta (20180205) Sebastian E. Ramos-Onsins.

Flags:
      -i [path and name of the tfa file (gz file indexed)]
      -g [path of the GFF_file]
         [add also: coding,noncoding,synonymous,nonsynonymous,silent, others (whatever annotated)]
         [if 'synonymous', 'nonsynonymous', 'silent' add: Genetic_Code: Nuclear_Universal,mtDNA_Drosophila,mtDNA_Mammals,Other]
         [if 'Other', introduce the single letter code for the 64 triplets in the order UUU UUC UUA UUG ... etc.]
      -c [in case use coding regions, criteria to consider transcripts (max/min/first/long)]. DEFAULT: long
      -n [name of each scaffold to analyze. tfa can be a list separated by commas(ex. -n chr1,chr2,chr3]
      -l [total length of each scaffold(s). if more than one, separated by commas]
   OPTIONAL PARAMETERS:
      -h [help and exit]
      -o [path and name of the output weighted file (must be ending with .gz)]
      -G [number of samples in the outgroup (if exist. Only allowed the last samples in the list)]. DEFAULT: 0
      -m [masking regions: file indicating the start and the end of regions to be masked by 0 weights]. DEFAULT: NONE
      -C [coordinates of regions: file indicating the start and the end of regions to be weighted (rest would be weighted as 0 if the file is included)]. DEFAULT: NONE

