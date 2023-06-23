/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zGeometricTail.h - part of the ZOE library for genomic analysis
 
 chaochun wei 

\******************************************************************************/

#ifndef ZOE_GEOMETRICTAIL_H
#define ZOE_GEOMETRICTAIL_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "zTools.h"

struct zGeoTailCell {   /* for explicit intron length distribtion 
						   with geometric distribution tail */ 
	score_t score;
	int from_state;  /* j */
	int state;    
	int length;   /* length = k > r, r is the length limit for explicit part */
};

typedef struct zGeoTailCell zGeoTailCell; 





struct zGeoTailCellVec {
	zGeoTailCell **elem;
	int size; 
	int limit; 
};
typedef struct zGeoTailCellVec zGeoTailCellVec;

void zInitGeoTailCellVec (zGeoTailCellVec*, int);
void zFreeGeoTailCellVec (zGeoTailCellVec*);
void zPushGeoTailCellVec (zGeoTailCellVec*, zGeoTailCell*);
void zSetGeoTailCellVec (zGeoTailCellVec*, int, zGeoTailCell*);
void zSetGeoTailCellVec_v (zGeoTailCellVec*, int, score_t, coor_t, int, int);


struct zGeoTailVec {

    zGeoTailCellVec **elem; 
    int size;
    int limit;
};

typedef struct zGeoTailVec zGeoTailVec;


void zInitGeoTailVec (zGeoTailVec*, int);
void zFreeGeoTailVec (zGeoTailVec*);
void zPushGeoTailVec (zGeoTailVec*, zGeoTailCellVec*);

void zCutGeoTailVec (zGeoTailVec*, int);
#endif
