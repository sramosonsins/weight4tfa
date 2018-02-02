//
//  main_weight4tfa.c
//  weighttfa
//
//  Created by Sebastian Ramos-Onsins on 28/08/17.
//  Copyright © 2017 Sebastian Ramos-Onsins. All rights reserved.
//

#include "main_weight4tfa.h"

int main(int argc, const char * argv[]) {

    char *chr_name;
    char chr_name_all[ MSP_MAX_NAME];
    char **chr_name_array;
    int nscaffolds;

    long int lentotal;
    char chr_length_all[MSP_MAX_NAME];
    char **chr_length_array;
    int nscaffolds2;
 
    char file_in[ MSP_MAX_FILENAME];
    char file_out[MSP_MAX_FILENAME];
    char file_GFF[MSP_MAX_FILENAME];
    char file_Wcoord[MSP_MAX_FILENAME];
    char file_masked[MSP_MAX_FILENAME];
    char file_log[MSP_MAX_FILENAME];

    FILE *file_input	= 	0;
    FILE *file_output	=	stdout;
    FILE *file_wcoor    =   0;
    FILE *file_msk   	= 	0;
    FILE *file_logerr   =   stdout;
    
    SGZip file_input_gz;
    SGZip file_output_gz;
    SGZip file_wcoor_gz;
    SGZip file_msk_gz;
    SGZip file_logerr_gz;

    struct SGZIndex file_input_gz_index;          /* This is the index for the input gz file. */
    struct SGZIndex file_output_gz_index;          /* This is the index for the weight gz file. */

    /* GFF variables */
    int 	gfffiles			= 0;
    char 	subset_positions[ MSP_MAX_GFF_WORDLEN ];
    char 	code_name[ MSP_GENETIC_CODETYPE_LEN ];
    char 	genetic_code[ MSP_GENCODE_COMBINATIONS ];
    char 	criteria_transcript[ MSP_GFF_CRITERIA_MSG_LEN ];

    int include_unknown = 1;
    int outgroup = 0;
    
    long int x;
    int i,j,k;
    
    long int nwindows,masked_nwindows;
    long int *wgenes;
    long int *masked_wgenes;

    /*finally the two important vectors to print (weights)*/
    double *vector_sizepos; /*weights*/
    double *vector_segrpos; /*reject syn/nsyn/sil variants when no assigned*/

    memset( file_out, 0, MSP_MAX_FILENAME);

    memset( chr_name_all, 0, MSP_MAX_NAME);
    memset( chr_length_all, 0, MSP_MAX_NAME);
    
    memset( file_Wcoord, 0, MSP_MAX_FILENAME);
    memset( file_masked, 0, MSP_MAX_FILENAME);

    /*collect arguments*/
    if(collect_arguments(argc,argv,file_in,file_out,file_log,file_GFF,file_Wcoord,file_masked,
                          &gfffiles,subset_positions,code_name,genetic_code,criteria_transcript,
                          &outgroup,chr_name_all,chr_length_all) == 0) {
        printf("\nError collecting arguments.\n");
        exit(1);
    }
    
    /* Opening files */
    if( file_in[0] == '\0' ) {
        fprintf(stdout,"\n It is not possible to read the input file %s\n", file_in);
        exit(1);
    }
    else {
        if( (file_input = fzopen( file_in, "r", &file_input_gz)) == 0) {
            fzprintf(file_logerr,&file_logerr_gz,"\n It is not possible to open the input file %s\n", file_in);
            exit(1);
        }
    }
    
    /*define the file for weigth for positions*/
    if(file_out[0] == '\0') {
        strcpy(file_out,file_in);
        strcat(file_out,"_WEIGHTS.gz");
    }
    if( (file_output = fzopen( file_out, "w", &file_output_gz)) == 0) {
        fprintf(stdout,"\n It is not possible to write in the output weight file %s\n", file_out);
        exit(1);
    }
    
    /*print all the argv. Header*/
    fzprintf(file_output,&file_output_gz,"#weight4tfa ");
    for(i=1;i<argc;i++) {
        fzprintf(file_output,&file_output_gz,"%s ",argv[i]);
    }
    fzprintf(file_output,&file_output_gz,"\n");
    
    /*initialize index file pointer*/
    /*read input index*/
    load_index_from_file(file_input_gz.index_file_name, &file_input_gz_index);
    /*write output index*/
    init_gzindex_structure(&file_output_gz_index);
    file_output_gz.index = &file_output_gz_index;
    
    strcpy(file_log, file_out);
    strcat(file_log,".log");
    if( (file_logerr = fzopen( file_log, "w", &file_logerr_gz)) == 0) {
        fprintf(stdout,"\n It is not possible to write the log file %s.", file_log);
        exit(1);
    }
    
    /*printf("\nOpen log file...");*/
    fzprintf(file_logerr,&file_logerr_gz,"\nOpen log file...");

    /*separate all values of the list chr_name_all in chr_name_array: */
    /* Only do the list if input and output is tfa*/
    nscaffolds = 1;
    chr_name_array = (char **)calloc(nscaffolds,sizeof(char *));
    chr_name_array[0] = (char *)calloc(MSP_MAX_NAME,sizeof(char));
    j=0;
    while(chr_name_all[j] != '\0') {
        k=0;
        while(chr_name_all[j] != ',' && chr_name_all[j] != '\0' && j < MSP_MAX_NAME) {
            chr_name_array[nscaffolds-1][k] = chr_name_all[j];
            j++; k++;
        }
        if(chr_name_all[j] == ',') {
            chr_name_array[nscaffolds-1][k] = '\0';
            nscaffolds += 1;
            chr_name_array = (char **)realloc(chr_name_array,nscaffolds*sizeof(char *));
            chr_name_array[nscaffolds-1] = (char *)calloc(MSP_MAX_NAME,sizeof(char));
            j++;
        }
    }
    nscaffolds2 = 1;
    chr_length_array = (char **)calloc(nscaffolds2,sizeof(char *));
    chr_length_array[0] = (char *)calloc(MSP_MAX_NAME,sizeof(char));
    j=0;
    while(chr_length_all[j] != '\0') {
        k=0;
        while(chr_length_all[j] != ',' && chr_length_all[j] != '\0' && j < MSP_MAX_NAME) {
            chr_length_array[nscaffolds2-1][k] = chr_length_all[j];
            j++; k++;
        }
        if(chr_length_all[j] == ',') {
            chr_length_array[nscaffolds2-1][k] = '\0';
            nscaffolds2 += 1;
            chr_length_array = (char **)realloc(chr_length_array,nscaffolds2*sizeof(char *));
            chr_length_array[nscaffolds2-1] = (char *)calloc(MSP_MAX_NAME,sizeof(char));
            j++;
        }
    }
    if(nscaffolds != nscaffolds2) {
        fzprintf(file_logerr,&file_logerr_gz,"\nERROR: the number of scaffolds is not coincident in options -l and -n.\n");
        exit(1);
    }
    
    /**************************************/
    /**************************************/
    /**************************************/
    /**************************************/
    /*               LOOP                 */
    /**************************************/
    /**************************************/
    /**************************************/
    /**************************************/
    
    nwindows = 0;
    masked_nwindows = 0;
    wgenes = 0;
    masked_wgenes = 0;

    for(i=0;i<nscaffolds;i++) {
        /* Open coordinates file */
        if( file_Wcoord[0] == '\0' ) {
            file_wcoor = 0;
            nwindows = 0;
        }
        else {
            if( (file_wcoor = fzopen( file_Wcoord, "r", &file_wcoor_gz)) == 0) {
                fzprintf(file_logerr,&file_logerr_gz,"\n It is not possible to open the coordinates file %s\n", file_Wcoord);
                exit(1);
            }
        }
        /* Opening mask coordinates file */
        if( file_masked[0] == '\0' ) {
            file_msk = 0;
            masked_nwindows = 0;
        }
        else {
            if( (file_msk = fzopen( file_masked, "r", &file_msk_gz)) == 0) {
                fzprintf(file_logerr,&file_logerr_gz,"\n It is not possible to open the masked coordinates file %s\n", file_masked);
                exit(1);
            }
        }
        
        /*scaffold to work*/
        lentotal = atol(chr_length_array[i]);
        chr_name = chr_name_array[i];
        
        /*define vectors for weigthts*/
        if((vector_sizepos = (double *)calloc((unsigned long)lentotal,sizeof(double))) == NULL ) {
            fzprintf(file_logerr,&file_logerr_gz,"\nError: memory not reallocated. vecor_sizepos \n");
            exit(1);
        }
        if((vector_segrpos = (double *)calloc((unsigned long)lentotal,sizeof(double))) == NULL ) {
            fzprintf(file_logerr,&file_logerr_gz,"\nError: memory not reallocated. vector_segrpos \n");
            exit(1);
        }
        for(x=0;x<lentotal;x++) {
            vector_sizepos[x] = (double)1;
            vector_segrpos[x] = (double)1;
        }


        /* read coordinates file */
        if( file_Wcoord[0] != '\0' ) {
            if(read_coordinates(file_wcoor,&file_wcoor_gz,file_output,&file_output_gz, file_logerr, &file_logerr_gz,&wgenes,&nwindows,chr_name) == 0) {
                fzprintf(file_logerr,&file_logerr_gz,"Error processing coordinates file %s\n", file_Wcoord);
                exit(1);
            }
            fzclose(file_wcoor, &file_wcoor_gz);
        }
        /* read mask coordinates file */
        if( file_masked[0] != '\0' ) {
            if(read_coordinates(file_msk,&file_msk_gz,file_output,&file_output_gz, file_logerr, &file_logerr_gz,&masked_wgenes,&masked_nwindows,chr_name) == 0) {
                fzprintf(file_logerr,&file_logerr_gz,"Error processing masked coordinates file %s\n", file_masked);
                exit(1);
            }
            fzclose(file_msk, &file_msk_gz);
        }
        
        /**************************************/
        /**************************************/
        /**************************************/
        /**************************************/
        
        /*read sorted GFF*/
        /*assign final weights according to tfa data at each specific gene*/
        if(assign_weights_from_gff(file_GFF,
                   file_input,&file_input_gz,&file_input_gz_index,
                   file_output,&file_output_gz,
                   file_logerr,&file_logerr_gz,
                   vector_sizepos,vector_segrpos,
                   include_unknown,subset_positions,genetic_code,criteria_transcript,
                   outgroup,chr_name,i,lentotal) == 0) {
            fzprintf(file_logerr,&file_logerr_gz,"\nError reading GFF file. No complete output file released.\n");
            unload_all_index_positions(&file_input_gz_index);
            fzclose(file_input, &file_input_gz);
            fzprintf(file_logerr,&file_logerr_gz,"\nProgram Ended\n");
            fzclose(file_logerr, &file_logerr_gz);
           exit(1);
        }
        
        /**************************************/
        /**************************************/
        /**************************************/
        /**************************************/

        /*print weights file*/
        if(i==0) fzprintf(file_output,&file_output_gz,"#CHR:POSITION\tWEIGHTPOS\tWEIGHTVAR_%s\n",subset_positions);
        for(x=0;x<lentotal;x++) {
            fzprintf(file_output,&file_output_gz,"%s:%ld\t",chr_name,x+1);
            fzprintf(file_output,&file_output_gz, "%.3f\t",vector_sizepos[x]);
            fzprintf(file_output,&file_output_gz,"%.3f\t",vector_segrpos[x]);
            fzprintf(file_output,&file_output_gz,"\n");
        }
        
        /*free all arrays*/
        if(file_wcoor) free(wgenes);
        if(file_msk) free(masked_wgenes);
        free(vector_sizepos);
        free(vector_segrpos);
    }

    fzclose(file_output, &file_output_gz);
    unload_all_index_positions(&file_input_gz_index);
    fzclose(file_input, &file_input_gz);
    
    fzprintf(file_logerr,&file_logerr_gz,"\nProgram Ended\n");
    fzclose(file_logerr, &file_logerr_gz);
    
    /**************************************/
    /**************************************/
    /**************************************/
    /**************************************/
    /*           END LOOP                 */
    /**************************************/
    /**************************************/
    /**************************************/
    /**************************************/

    for(i=0;i<nscaffolds;i++) {
        free(chr_name_array[i]);
        free(chr_length_array[i]);
    }
    free(chr_name_array);
    free(chr_length_array);

    return(0);
}

void usage(void)
{
    printf(WEIGHT4TFA);
    printf("\nFlags:\n");
    printf("      -i [path and name of the tfa file (gz file indexed)]\n");
    printf("      -g [path of the GFF_file]\n");
    printf("         [add also: coding,noncoding,synonymous,nonsynonymous,silent, others (whatever annotated)]\n");
    printf("         [if 'synonymous', 'nonsynonymous', 'silent' add: Genetic_Code: Nuclear_Universal,mtDNA_Drosophila,mtDNA_Mammals,Other]\n");
    printf("         [if 'Other', introduce the single letter code for the 64 triplets in the order UUU UUC UUA UUG ... etc.]\n");
    printf("      -c [in case use coding regions, criteria to consider transcripts (max/min/first/long)]. DEFAULT: long\n");
    printf("      -n [name of each scaffold to analyze. tfa can be a list separated by commas(ex. -n chr1,chr2,chr3]\n");
    printf("      -l [total length of each scaffold(s). if more than one, separated by commas]\n");
    printf("   OPTIONAL PARAMETERS:\n");
    printf("      -h [help and exit]\n");
    printf("      -o [path and name of the output weighted file (must be ending with .gz)]\n");
    printf("      -G [number of samples in the outgroup (if exist. Only allowed the last samples in the list)]. DEFAULT: 0\n");
    printf("      -m [masking regions: file indicating the start and the end of regions to be masked by 0 weights]. DEFAULT: NONE\n");
    printf("      -C [coordinates of regions: file indicating the start and the end of regions to be weighted (rest would be weighted as 0 if the file is included)]. DEFAULT: NONE\n");
    printf("\n");
}

int collect_arguments(int argc, const char **argv,
                      char *file_in,
                      char *file_out,
                      char *file_log,
                      char *file_GFF,
                      char *file_Wcoord,
                      char *file_masked,
                      int   *gfffiles,
                      char 	*subset_positions,
                      char 	*code_name,
                      char 	*genetic_code,
                      char 	*criteria_transcript,
                      int   *outgroup,
                      char  *chr_name_all,
                      char *chr_length_all) {
    
    int arg;
    int x;
    
    if( argc > 1 )
    {
        arg = 1;
        while(arg < argc)
        {
            if( argv[arg][0] != '-' )
            {
                if(argv[arg][0] == '>')
                    break;
                printf(" argument should be -%s ?\n", argv[arg]);
                usage();
                exit(1);
            }
            
            switch (argv[arg][1])
            {
                case 'i' : /* input file*/
                    arg++;
                    strcpy( file_in, argv[arg] );
                    break;
                case 'o' : /* output file */
                    arg++;
                    strcpy( file_out, argv[arg] );
                    break;
                case 'g': /* g GFF file name, AND more words
                           2nd : synonymous, nonsynonymous, silent or whatever
                           3rd : selected genetic code name or "Other"
                           next 64th's : in case of 'Other', provide 64 IUPAC letters of each codon.
                           * Check order.
                           */
                    arg++;
                    strcpy( file_GFF, argv[arg]);
                    arg++;
                    strcpy( subset_positions, argv[arg] );
                    
                    *gfffiles = 1;
                    
                    /* Go if known coding option - */
                    if(( strcmp(subset_positions,"synonymous")==0 ||
                         strcmp(subset_positions,"nonsynonymous")==0 ||
                         strcmp(subset_positions,"silent")==0))
                    {
                        arg++;
                        strcpy( code_name,argv[arg] );
                        
                        if( strcmp(code_name, "Nuclear_Universal") == 0 )
                        {
                            genetic_code[0] = 'F';
                            genetic_code[1] = 'F';
                            genetic_code[2] = 'L';
                            genetic_code[3] = 'L';
                            genetic_code[4] = 'S';
                            genetic_code[5] = 'S';
                            genetic_code[6] = 'S';
                            genetic_code[7] = 'S';
                            genetic_code[8] = 'Y';
                            genetic_code[9] = 'Y';
                            genetic_code[10] = '*';
                            genetic_code[11] = '*';
                            genetic_code[12] = 'C';
                            genetic_code[13] = 'C';
                            genetic_code[14] = '*';
                            genetic_code[15] = 'W';
                            genetic_code[16] = 'L';
                            genetic_code[17] = 'L';
                            genetic_code[18] = 'L';
                            genetic_code[19] = 'L';
                            genetic_code[20] = 'P';
                            genetic_code[21] = 'P';
                            genetic_code[22] = 'P';
                            genetic_code[23] = 'P';
                            genetic_code[24] = 'H';
                            genetic_code[25] = 'H';
                            genetic_code[26] = 'Q';
                            genetic_code[27] = 'Q';
                            genetic_code[28] = 'R';
                            genetic_code[29] = 'R';
                            genetic_code[30] = 'R';
                            genetic_code[31] = 'R';
                            genetic_code[32] = 'I';
                            genetic_code[33] = 'I';
                            genetic_code[34] = 'I';
                            genetic_code[35] = 'M';
                            genetic_code[36] = 'T';
                            genetic_code[37] = 'T';
                            genetic_code[38] = 'T';
                            genetic_code[39] = 'T';
                            genetic_code[40] = 'N';
                            genetic_code[41] = 'N';
                            genetic_code[42] = 'K';
                            genetic_code[43] = 'K';
                            genetic_code[44] = 'S';
                            genetic_code[45] = 'S';
                            genetic_code[46] = 'R';
                            genetic_code[47] = 'R';
                            genetic_code[48] = 'V';
                            genetic_code[49] = 'V';
                            genetic_code[50] = 'V';
                            genetic_code[51] = 'V';
                            genetic_code[52] = 'A';
                            genetic_code[53] = 'A';
                            genetic_code[54] = 'A';
                            genetic_code[55] = 'A';
                            genetic_code[56] = 'D';
                            genetic_code[57] = 'D';
                            genetic_code[58] = 'E';
                            genetic_code[59] = 'E';
                            genetic_code[60] = 'G';
                            genetic_code[61] = 'G';
                            genetic_code[62] = 'G';
                            genetic_code[63] = 'G';
                        }
                        else if(strcmp(code_name,"mtDNA_Drosophila")==0)
                        {
                            genetic_code[0] = 'F';
                            genetic_code[1] = 'F';
                            genetic_code[2] = 'L';
                            genetic_code[3] = 'L';
                            genetic_code[4] = 'S';
                            genetic_code[5] = 'S';
                            genetic_code[6] = 'S';
                            genetic_code[7] = 'S';
                            genetic_code[8] = 'Y';
                            genetic_code[9] = 'Y';
                            genetic_code[10] = '*';
                            genetic_code[11] = '*';
                            genetic_code[12] = 'C';
                            genetic_code[13] = 'C';
                            genetic_code[14] = 'W';
                            genetic_code[15] = 'W';
                            genetic_code[16] = 'L';
                            genetic_code[17] = 'L';
                            genetic_code[18] = 'L';
                            genetic_code[19] = 'L';
                            genetic_code[20] = 'P';
                            genetic_code[21] = 'P';
                            genetic_code[22] = 'P';
                            genetic_code[23] = 'P';
                            genetic_code[24] = 'H';
                            genetic_code[25] = 'H';
                            genetic_code[26] = 'Q';
                            genetic_code[27] = 'Q';
                            genetic_code[28] = 'R';
                            genetic_code[29] = 'R';
                            genetic_code[30] = 'R';
                            genetic_code[31] = 'R';
                            genetic_code[32] = 'I';
                            genetic_code[33] = 'I';
                            genetic_code[34] = 'M';
                            genetic_code[35] = 'M';
                            genetic_code[36] = 'T';
                            genetic_code[37] = 'T';
                            genetic_code[38] = 'T';
                            genetic_code[39] = 'T';
                            genetic_code[40] = 'N';
                            genetic_code[41] = 'N';
                            genetic_code[42] = 'K';
                            genetic_code[43] = 'K';
                            genetic_code[44] = 'S';
                            genetic_code[45] = 'S';
                            genetic_code[46] = 'S';
                            genetic_code[47] = 'S';
                            genetic_code[48] = 'V';
                            genetic_code[49] = 'V';
                            genetic_code[50] = 'V';
                            genetic_code[51] = 'V';
                            genetic_code[52] = 'A';
                            genetic_code[53] = 'A';
                            genetic_code[54] = 'A';
                            genetic_code[55] = 'A';
                            genetic_code[56] = 'D';
                            genetic_code[57] = 'D';
                            genetic_code[58] = 'E';
                            genetic_code[59] = 'E';
                            genetic_code[60] = 'G';
                            genetic_code[61] = 'G';
                            genetic_code[62] = 'G';
                            genetic_code[63] = 'G';
                        }
                        else if( strcmp(code_name,"mtDNA_Mammals") == 0 )
                        {
                            genetic_code[0] = 'F';
                            genetic_code[1] = 'F';
                            genetic_code[2] = 'L';
                            genetic_code[3] = 'L';
                            genetic_code[4] = 'S';
                            genetic_code[5] = 'S';
                            genetic_code[6] = 'S';
                            genetic_code[7] = 'S';
                            genetic_code[8] = 'Y';
                            genetic_code[9] = 'Y';
                            genetic_code[10] = '*';
                            genetic_code[11] = '*';
                            genetic_code[12] = 'C';
                            genetic_code[13] = 'C';
                            genetic_code[14] = 'W';
                            genetic_code[15] = 'W';
                            genetic_code[16] = 'L';
                            genetic_code[17] = 'L';
                            genetic_code[18] = 'L';
                            genetic_code[19] = 'L';
                            genetic_code[20] = 'P';
                            genetic_code[21] = 'P';
                            genetic_code[22] = 'P';
                            genetic_code[23] = 'P';
                            genetic_code[24] = 'H';
                            genetic_code[25] = 'H';
                            genetic_code[26] = 'Q';
                            genetic_code[27] = 'Q';
                            genetic_code[28] = 'R';
                            genetic_code[29] = 'R';
                            genetic_code[30] = 'R';
                            genetic_code[31] = 'R';
                            genetic_code[32] = 'I';
                            genetic_code[33] = 'I';
                            genetic_code[34] = 'M';
                            genetic_code[35] = 'M';
                            genetic_code[36] = 'T';
                            genetic_code[37] = 'T';
                            genetic_code[38] = 'T';
                            genetic_code[39] = 'T';
                            genetic_code[40] = 'N';
                            genetic_code[41] = 'N';
                            genetic_code[42] = 'K';
                            genetic_code[43] = 'K';
                            genetic_code[44] = 'S';
                            genetic_code[45] = 'S';
                            genetic_code[46] = '*';
                            genetic_code[47] = '*';
                            genetic_code[48] = 'V';
                            genetic_code[49] = 'V';
                            genetic_code[50] = 'V';
                            genetic_code[51] = 'V';
                            genetic_code[52] = 'A';
                            genetic_code[53] = 'A';
                            genetic_code[54] = 'A';
                            genetic_code[55] = 'A';
                            genetic_code[56] = 'D';
                            genetic_code[57] = 'D';
                            genetic_code[58] = 'E';
                            genetic_code[59] = 'E';
                            genetic_code[60] = 'G';
                            genetic_code[61] = 'G';
                            genetic_code[62] = 'G';
                            genetic_code[63] = 'G';
                        }
                        else if( strcmp(code_name,"Other") == 0 ) {
                            for(x=0;x<64;x++) {
                                arg++;
                                if(argv[arg][0] == '-') {
                                    printf("\n Error in -g argument: In case use \"Other\", include the genetic code of the 64 aa values.");
                                    usage();
                                    exit(1);
                                }
                                genetic_code[x] = atoi(argv[arg]);
                            }
                        }
                        else {
                            printf(" %s: Unknown code, sorry", code_name);
                            exit(1);
                        }	
                    }
                    break;
                case 'c' : /* c Criteria used for analyzing the transcripts */
                    /* Basically, max or min */
                    arg++;
                    strcpy(criteria_transcript,argv[arg]);
                    if(strcmp( criteria_transcript, "max")!=0  && 
                       strcmp( criteria_transcript, "min")!=0  && 
                       strcmp( criteria_transcript, "first")!=0   && 
                       strcmp( criteria_transcript, "long")!=0  ) 
                    {
                        printf("\n Error: the argument -c has only the choices 'max', 'min', 'first' or 'long'.");
                        usage();
                        exit(1);
                    }
                    break;
                case 'G' : /* outgroup present or not */
                    arg++;
                    *outgroup = (int)atoi(argv[arg]);
                    if(*outgroup < 0) {
                        printf("\n Error in -o argument: only  or positive values are allowed.");
                        usage();
                        exit(1);
                    }
                    break;
                case 'W' : /* file with the coordinates of each window [init end](overwrite options -w and -s)*/
                    arg++;
                    strcpy(file_Wcoord, argv[arg] );					
                    break;
                case 'm' : /* file with the coordinates of each window [init end] to be masked by Ns*/
                    arg++;
                    strcpy(file_masked, argv[arg] );
                    break;
                case 'n' : /* name of the scaffold to analyze*/
                    arg++;
                    strcpy(chr_name_all, argv[arg] );
                    break;
                case 'l' : // total lengths of the scaffold(s) to analyze
                    arg++;
                    strcpy( chr_length_all, argv[arg] );
                    break;
                case 'h' : /* h HELP */
                    usage();
                    exit(0);
                    break;
            }
            arg++;
        }
    }
    else {
        usage();
        exit(0);
    }

    return(1);
}

int read_coordinates(FILE *file_wcoor, SGZip *file_wcoor_gz, FILE *file_output, SGZip *file_output_gz, FILE *file_logerr, SGZip *file_logerr_gz, long int **wgenes, long int *nwindows,char *chr_name) {
    
    char *valn=0;
    int c;
    int x;
    long int xx;
    long int dd;
    double ee;
    int inside_chr;
    
    /*printf("\nReading coordinates file...");*/
    fflush(stdout);
    fzprintf(file_logerr,file_logerr_gz,"\nReading coordinates file...");
    
    if((valn = (char *)calloc(100,sizeof(char))) == 0) {
        fzprintf(file_logerr,file_logerr_gz,"\nError: memory not reallocated. read_coordinates.00 \n");
        return(0);
    }
    if((*wgenes = (long int *)calloc(10000,sizeof(long int))) == 0) {
        fzprintf(file_logerr,file_logerr_gz,"\nError: memory not reallocated. read_coordinates.0 \n");
        free(valn);
        return(0);
    }
    
    *nwindows = 0;
    c = fzgetc(file_wcoor, file_wcoor_gz);
    
    if(check_comment(&c,file_wcoor,file_wcoor_gz) == 0) {
        fzprintf(file_logerr,file_logerr_gz,"\nWarning: no coordinates assigned. \n");
        free(valn);
        return(0);
    }
    inside_chr = 0;
    /*
     if(c==0 || c==-1 || c < 1 || c=='\xff' || c=='\xfe') {
     fzprintf(file_logerr,file_logerr_gz,"\nWarning: no coordinates assigned. \n");
     free(*wgenes);
     file_wcoor=0;
     return(0);
     }
     else {*/
    xx=0;
    while (c != 0 && c != -1 && c!='\xff' && c!='\xfe') {
        /*now keep all values: two columns, only numbers*/
        /*first column is the name of the scaffold*/
        while(c == 32 || c == 9 || c == 13 || c == 10)
            c = fzgetc(file_wcoor, file_wcoor_gz);
        if(check_comment(&c,file_wcoor, file_wcoor_gz) == 0) {
            c=0;break;
        }
        while(c == 32 || c == 9 || c == 13 || c == 10)
            c = fzgetc(file_wcoor, file_wcoor_gz);
        x=0;
        while(c != 32 && c != 9 && c != 13 && c != 10 && c!=0 && c!=-1 && x < 100  && c!='\xff' && c!='\xfe' && c>0) {
            valn[x] = c;
            c = fgetc(file_wcoor);
            x++;
        }
        valn[x] = '\0';/*scaffold name*/
        if(c==0 || c==-1 || c < 1  || c=='\xff' || c=='\xfe') {
            c=0;break;
        }
        while(strcmp(chr_name,valn) != 0) { /*discard those rows having different scaffold than chr_name*/
            if(inside_chr == 1) {
                inside_chr = -1; /*we passed target region*/
                break;
            }
            while(c != 13 && c != 10 && c!=0 && c!=-1  && c!='\xff' && c!='\xfe')
                c = fzgetc(file_wcoor, file_wcoor_gz);
            if(check_comment(&c,file_wcoor, file_wcoor_gz) == 0){
                free(valn);*nwindows = (xx)/2;return(1);
            }
            while(c == 32 || c == 9 || c == 13 || c == 10)
                c = fzgetc(file_wcoor, file_wcoor_gz);
            x=0;
            while(c != 32 && c != 9 && c != 13 && c != 10 && c!=0 && c!=-1 && x < 100  && c!='\xff' && c!='\xfe' && c>0) {
                valn[x] = c;
                c = fgetc(file_wcoor);
                x++;
            }
            valn[x] = '\0';/*scaffold name*/
            if(c==0 || c==-1 || c < 1  || c=='\xff' || c=='\xfe') {
                c=0;break;
            }
        }
        if(c==-1 || c == 0)
            break;
        if(inside_chr == -1)
            break; /*we go out*/
        inside_chr = 1;
        /*KEEP POSITIONS (first initial position, then end, next region and so on)*/
        while(c == 32 || c == 9 || c == 13 || c == 10)
            c = fzgetc(file_wcoor, file_wcoor_gz);
        if(check_comment(&c,file_wcoor, file_wcoor_gz) == 0) {
            c=0;free(valn);return(0);
        }
        while(c == 32 || c == 9 || c == 13 || c == 10)
            c = fzgetc(file_wcoor, file_wcoor_gz);
        x=0;
        while(c != 32 && c != 9 && c != 13 && c != 10 && c!=0 && c!=-1 && x < 100 && c!='\xff' && c!='\xfe' && c>0) {
            valn[x] = c;
            c = fzgetc(file_wcoor, file_wcoor_gz);
            x++;
        }
        valn[x] = '\0';
        wgenes[0][xx] = (long int)atof(valn);
        
        xx++;
        while(c == 32 || c == 9 || c == 13 || c == 10)
            c = fgetc(file_wcoor);
        if(check_comment(&c,file_wcoor, file_wcoor_gz) == 0) {
            c=0;free(valn);return(0);
        }
        while(c == 32 || c == 9 || c == 13 || c == 10)
            c = fzgetc(file_wcoor, file_wcoor_gz);
        x=0;
        while(c != 32 && c != 9 && c != 13 && c != 10 && c!=0 && c!=-1 && x < 100 && c!='\xff' && c!='\xfe' && c>0) {
            valn[x] = c;
            c = fgetc(file_wcoor);
            x++;
        }
        valn[x] = '\0';
        wgenes[0][xx] = (long int)round((double)atof(valn));
        
        xx++;
        dd = (long int)floor((double)xx/(double)10000);
        ee = (double)xx/(double)10000;
        if(dd==ee) {
            if((*wgenes = realloc(*wgenes,((long int)(dd+1)*(long int)10000*(long int)sizeof(long int)))) == 0) {
                puts("Error: realloc error read_coordinates.1\n");
                free(valn);/*free(*wgenes);*/
                return(0);
            }
        }
    }
    if(xx == 0) {
        fzprintf(file_logerr,file_logerr_gz,"\nError: no coordinates assigned, read_coordinates.2 \n");
        /*free(*wgenes);
         file_wcoor=0;*/
        return(0);
    }
    *nwindows = (xx)/2;
    /*}*/
    free(valn);
    return 1;
}

int  check_comment(int *c, FILE *file_input, SGZip *file_input_gz) {
    if(*c=='#') {
        /*parse comments*/
        while(*c!=10 && *c!=13 && *c != 0 && *c!=-1 && *c!='\xff' && *c!='\xfe') {
            *c = fzgetc(file_input, file_input_gz);
        }
        if(*c==0 || *c==-1 || *c =='\xff' || *c=='\xfe')
            return(0);
    }
    if(*c==0 || *c==-1 || *c =='\xff' || *c=='\xfe')
        return(0);
    
    return(1);
}

