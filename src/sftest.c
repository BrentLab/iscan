/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/*****************************************************************************\
sftest.c

sequence feature test program for zoe
\*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ZOE.h"

int main (int argc, char *argv[]) {
	
	FILE      *stream;
	zSfeature f;
	zSFVec    sfv;
	int       i;
	zVec      vec;
	
	/* set the program name */
	zSetProgramName(argv[0]);
	
	/* process commandline */
	if (argc < 2) {
		(void)fprintf(stderr, "usage: %s <file of zSfeature>\n", argv[0]);
		exit(1);
	}
	
	/* set up a vector with initial size 10 */
	zInitSFVec(&sfv, 10);
	
	/* read sequence features into a vector, keep pointers too */
	if ((stream = fopen(argv[1], "r")) == NULL)
		zDie("file error (%s)", argv[1]);
	while (zReadSfeature(stream, &f)) {
		zPushSFVec(&sfv, &f);
		zFreeSfeature(&f);
	}
	fclose(stream);
	
	/* write in same order as they came in */
	printf("Original order\n");
	for (i = 0; i < sfv.size; i++) {
		zWriteSfeature(stdout, &sfv.elem[i]);
	}
	
	/* sort for fun */
	printf("Sorted order\n");
	zInitVec(&vec, sfv.size);
	for (i = 0; i < sfv.size; i++) {
		zPushVec(&vec, &sfv.elem[i]);
	}
	qsort(vec.elem, vec.size, sizeof(void*), zSfPtrCmp);
	for (i = 0; i < vec.size; i++) {
		zWriteSfeature(stdout, vec.elem[i]);
	}
	
	zFreeSFVec(&sfv);

	return 0;
}
