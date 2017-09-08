/*
 * assign_weights_from_gff.h
 *
*/

#ifndef WEIGHTS_FROM_GFF_H_
#define WEIGHTS_FROM_GFF_H_

#include "common.h"

#ifdef	__cplusplus
extern "C" {
#endif

#include <stdio.h> /* para BUFSIZ??*/
#include "zutil.h"
#include "zindex.h"

struct valuesgff
{
	char filename[256];
	char source[256];
	char feature[256];
	char strand[1];
	long int start;
	long int end;
	char score[256];
	char frame[1];
	char seqname[256];
	char gene_id[256];
	char transcript_id[256];
};

int assign_weights_from_gff(char *name_fileinputgff,
                            FILE *file_input,
                            SGZip *file_input_gz,
                            struct SGZIndex *file_input_gz_index,
                            FILE *file_output,
                            SGZip *file_output_gz,
                            FILE *file_logerr,
                            SGZip *file_logerr_gz,
                            double *matrix_sizepos,
                            double *matrix_segrpos,
                            int include_unknown,
                            char *subset_positions,
                            char *genetic_code,
                            char *criteria_transcript,
                            int outgroup_presence,
                            char *chr_name,
                            int first,
                            long int n_site);

int tripletnsamp(char *cod3n,char *DNA_matr,char strand,double *cmat,
				 int n_samp,long int n_site,long int end,long int ii/*,
				 FILE *file_output*//*,int mainargc*/,int include_unknown,int type_output,
				 /*long int *nmhits, long int *mhitbp,*/ int outgroup_presence, int nsamoutg,
                 FILE *file_logerr, SGZip *file_logerr_gz);
	
int comp_trcpt_id(const void *a,const void *b);
int comp_start_id(const void *a,const void *b);
int comp_end_id(const void *a,const void *b);
int comp_gene_id(const void *a,const void *b);
/*do a function that read the genetic code and add the degeneration to each position at each codon*/
int function_do_nfold_triplets(int n_fold[64][3], char *genetic_code, char tripletsN[64][3]);

int function_read_tfasta(FILE *file_input,
                         SGZip *file_input_gz,
                         struct SGZIndex *index_input,
                         FILE *file_logerr,
                         SGZip *file_logerr_gz,
                         long int init_site,
                         long int end_site,
                         int n_sam,
                         long int *n_site,
                         char **DNA_matr,
                         char *chr_name,
                         int first);
int transform_beg_chr(char *ID, char *chr_name, long int beg,int nchstr, int count0s);
int check_comment(char *c, FILE *file_input, SGZip *file_input_gz);
    
#ifdef	__cplusplus
}
#endif



#endif /* WEIGHTS_FROM_GFF_H_ */
