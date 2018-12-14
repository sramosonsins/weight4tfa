//
//  main_weight4tfa.h
//  weighttfa
//
//  Created by Sebastian Ramos-Onsins on 28/08/17.
//  Copyright © 2017 Sebastian Ramos-Onsins. All rights reserved.
//

#ifndef main_weight4tfa_h
#define main_weight4tfa_h

#include <stdio.h>

#ifdef	__cplusplus
extern "C" {
    #endif
    
    #include "common.h"
    #include "zutil.h"
    #include "zindex.h"
    
    void usage(void);
    int collect_arguments(int argc, const char **argv,
                          char *file_in,
                          char *file_out,
                          char *file_log,
                          char *file_GFF,
                          char *file_Wcoord,
                          char *file_masked,
                          int  *gfffiles,
                          char *subset_positions,
                          char *code_name,
                          char *genetic_code,
                          char *criteria_transcript,
                          int  *outgroup,
                          char *chr_name_all,
                          char *chr_length_all);
    int read_coordinates(FILE *file_wcoor,
                         SGZip *file_wcoor_gz,
                         FILE *file_output,
                         SGZip *file_output_gz,
                         FILE *file_logerr,
                         SGZip *file_logerr_gz,
                         long int **wgenes,
                         long int *nwindows,
                         char *chr_name);
    int assign_weights_from_gff(char *name_fileinputgff,
                               FILE *file_input,
                               SGZip *file_input_gz,
                               struct SGZIndex *file_input_gz_index,
                               FILE *file_output,
                               SGZip *file_output_gz,
                               FILE *file_logerr,
                               SGZip *file_logerr_gz,
                               double *vector_sizepos,
                               double *vector_segrpos,
                               int include_unknown,
                               char *subset_positions,
                               char *genetic_code,
                               char *criteria_transcript,
                               int outgroup,
                               char *chr_name,
                               unsigned long i,
                               long int n_site);
    int check_comment(int *c, FILE *file_input, SGZip *file_input_gz);
    int read_index_file(char *chr_name_all, unsigned long *nscaffolds,char ***chr_name_array,char ***chr_length_array);
    #ifdef	__cplusplus
}
#endif

#endif /* main_weight4tfa_h */
