/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
 zGeometricTail.c - part of the ZOE library for genomic analysis
 
 chaochun wei 

\******************************************************************************/


#ifndef ZOE_GEOMETRICTAIL_C
#define ZOE_GEOMETRICTAIL_C

#include "zGeometricTail.h" 

void zInitGeoTailCellVec (zGeoTailCellVec* gtcv, int limit){
    int i = 0;

    gtcv->size  = 0;
    gtcv->limit = 0;


    if (limit > 0) {
		gtcv->elem  = zCalloc(limit, sizeof(zGeoTailCell*), "zInitGeoTailCellVec\n");
		gtcv->size = limit;
		gtcv->limit = limit;
		for(i = 0; i < limit; i ++) {
			gtcv->elem[i]= zMalloc(sizeof(zGeoTailCell), "zMalloc GeoTailCell\n");
			zSetGeoTailCellVec_v(gtcv, i, MIN_SCORE, -1, 0, -1);
		}
    }
    else 
		gtcv->elem = NULL;
	  

}

void zFreeGeoTailCellVec (zGeoTailCellVec* gtcv){
	int i;
	for(i = 0; i < gtcv->size; i ++){
		zFree(gtcv->elem[i]);
	}
    zFree(gtcv->elem);
	gtcv->elem = NULL;
	zFree(gtcv);
	gtcv = NULL;
}

void zSetGeoTailCellVec (zGeoTailCellVec* gtcv, int i,  zGeoTailCell* s ) {
    if(gtcv->limit < i || gtcv->size < i ) { zDie(" array offset out of limit\n");}
    else {
		gtcv->elem[i]->score = s->score;
		gtcv->elem[i]->from_state = s->from_state;
		gtcv->elem[i]->length = s->length;
    }

}

void zSetGeoTailCellVec_v (zGeoTailCellVec *gtcv, int i, score_t score, coor_t length, int from_state, int state ) {
    if(gtcv->limit < i || gtcv->size < i) { zDie(" array offset out of limit\n");}
    else {
		gtcv->elem[i]->score      = score;
		gtcv->elem[i]->from_state = from_state;
		gtcv->elem[i]->length     = length;
		gtcv->elem[i]->state      = state;
    }

}

void zPushGeoTailCellVec (zGeoTailCellVec* gtcv, zGeoTailCell *thing ) {
	if (gtcv->limit == gtcv->size) {
		if (gtcv->limit == 0) gtcv->limit  = 1;
		else                 gtcv->limit *= 2;
		gtcv->elem = zRealloc(gtcv->elem, gtcv->limit * sizeof(zGeoTailCell *), 
				      "zPushGeoTailCellVec");
	}
	gtcv->elem[gtcv->size] = thing;
	gtcv->size++;
}




void zInitGeoTailVec (zGeoTailVec* gtv, int limit){
	gtv->size  = limit;
	gtv->limit = limit;
	if (limit > 0) gtv->elem = zMalloc(limit * sizeof(zGeoTailCellVec*), "zInitGeoTailVec");
	else                gtv->elem = NULL;
}

void zFreeGeoTailVec (zGeoTailVec* gtv){ /* only free the space in gtv.elem */
	int i;
	if(gtv != NULL) {
		if(gtv->elem != NULL) {
			for(i = 0; i < gtv->size; i ++){
				zFreeGeoTailCellVec (gtv->elem[i]);
				gtv->elem[i] = NULL;
			}
			zFree(gtv->elem);
			fprintf(stderr, "free: %d\n",(int) ((gtv->size - 1) * sizeof(zGeoTailCellVec*)));
			gtv->size = 0;
		}
	}
	
}


void zPushGeoTailVec (zGeoTailVec* gtv, zGeoTailCellVec *thing) {

	if (gtv->limit == gtv->size) {
		if (gtv->limit == 0) gtv->limit  = 1;
		else                 gtv->limit *= 2;
		fprintf(stderr,  "realloc %d\n", (int)(gtv->limit*sizeof(zGeoTailCellVec*)));
		gtv->elem = zRealloc(gtv->elem, gtv->limit * sizeof(zGeoTailCellVec*), "zPushGeoTailVec");
	}
	gtv->elem[gtv->size] = thing;
	gtv->size++;
}




void zCutGeoTailVec(zGeoTailVec* gtv,  int len) {
    int p;
    if(gtv->size > len ) {
		if(len > 0) {
			for(p = len ; p < gtv->size; p ++ ){
				zFreeGeoTailCellVec(gtv->elem[p]);
				gtv->elem[p]=NULL;
			}

			gtv->size = len;

		}
		else{
			if (gtv->elem != NULL) {
				gtv->elem[0] = NULL;  /* keep elem[0] , free all the other space */
				for(p = 1 ; p < gtv->size; p ++ ){
					zFree(gtv->elem[p]);
					gtv->elem[p]=NULL;
				}
				gtv->size = 0;
			}
		}
    }
	

}
#endif
