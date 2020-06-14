/*
 *  common.h
 *
 *  Created by Sebastian E. Ramos Onsins on 27/11/2012.
 *
 */

#ifndef COMMON_
#define COMMON_

#ifdef	__cplusplus
extern "C" {
	#endif
	
	#define _CRT_SECURE_NO_DEPRECATE
	
	#include <stdio.h>
	#include <string.h>
	#include <stdlib.h>
	#include <math.h>
	
	#define WEIGHT4TFA "#weight4tfa v0.1beta (20200614) Sebastian E. Ramos-Onsins.\n"

	#define MSP_MAX_FILENAME			(unsigned long) 4096 /**< @brief Maximum Filename Length allowed */
	#define MSP_MAX_GFF_WORDLEN         (unsigned long) 20
	#define MSP_GENETIC_CODETYPE_LEN	(unsigned long) 50	/* e.g. "Nuclear universal" */
	#define MSP_GENCODE_COMBINATIONS    (unsigned long) 64 /* 4^3 */
	#define MSP_GFF_CRITERIA_MSG_LEN    (unsigned long) 20 /* e.g. "MIN" */
    #define MSP_MAX_FILELINE_LEN		(unsigned long) 4096
    #define MSP_MAX_NAME                (unsigned long) 4096
	#define CHI_INTERVAL                (unsigned long) 10
    #define MSP_MAX_COL_LINE            (unsigned long) 32767*2
    #define MAXLEN                      (unsigned long) 1000
    
	#ifndef NULL
	#define NULL	0
	#endif
	
	#ifdef	__cplusplus
	}
#endif

#endif /* COMMON */
