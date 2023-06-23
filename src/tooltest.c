/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/*****************************************************************************\
tooltest.c

testing some of the generic tools

\*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "ZOE.h"

void usage (void);
void testInt (FILE*);
void testFloat (FILE*);
void testText (FILE*);
void testTable (FILE*);
void testIntern (FILE*);
void testOpen (void);

/*****************************************************************************\
 Main Program
\*****************************************************************************/

int main (int argc, char *argv[]) {
	FILE* stream;
	
	zSetProgramName(argv[0]);
	zParseOptions(&argc, argv);

	if (argc != 2) usage();
	if ((stream= fopen(argv[2], "r")) == NULL) zDie("error opening %s", argv[2]);

	if (zOption("int")) testInt(stream);
	if (zOption("float")) testFloat(stream);
	if (zOption("text")) testText(stream);
	if (zOption("table")) testTable(stream);
	if (zOption("intern")) testIntern(stream);

	fclose(stream);

	return 0;
}

void usage (void) {
	(void)fprintf(stderr, "usage: tooltest <type> <file>\n");
	(void)fprintf(stderr, "type: -int, -float, -text, -table\n");
	exit(1);
}

void testInt (FILE *stream) {
	int i;
	int val;
	zIVec vec;
	
	zInitIVec(&vec, 0);
	while (fscanf(stream, "%d", &val) == 1) {
		zPushIVec(&vec, val);
	}
	qsort(vec.elem, vec.size, sizeof(int), zIcmp);
	for (i = 0; i < vec.size; i++) {
		(void)printf("%d\n", vec.elem[i]);
	}
	zFreeIVec(&vec);
}

void testFloat (FILE *stream) {
	int i;
	float val;
	zFVec vec;
	
	zInitFVec(&vec, 0);
	while (fscanf(stream, "%f", &val) == 1) {
		zPushFVec(&vec, val);
	}
	qsort(vec.elem, vec.size, sizeof(float), zFcmp);
	for (i = 0; i < vec.size; i++) {
		(void)printf("%f\n", vec.elem[i]);
	}
	zFreeFVec(&vec);
}

void testText (FILE *stream) {
	int i;
	char  val[256];
	zTVec vec;
	
	zInitTVec(&vec, 0);
	while (fscanf(stream, "%s", val) == 1) {
		zPushTVec(&vec, val);
	}
	qsort(vec.elem, vec.size, sizeof(int), zTcmp);
	for (i = 0; i < vec.size; i++) {
		(void)printf("%s\n", vec.elem[i]);
	}
	zFreeTVec(&vec);
}

void testTable (FILE *stream) {
	int    i;
	char   line[256], key[64], val[64];
	zTVec  kvec;
	zTVec  vvec;
	zHash  hash;
	zVec  *keys;
	char  *k, *v;
	
	zInitTVec(&kvec, 10);
	zInitTVec(&vvec, 10);
	zInitHash(&hash);
	
	while (fgets(line, sizeof(line), stream)) {
		if (sscanf(line, "%s %s", key, val) != 2) {
			zDie("table format violated");
		}
		zPushTVec(&kvec, key);
		zPushTVec(&vvec, val);
		zSetHash(&hash, kvec.last, vvec.last);
	}
	
	keys = zKeysOfHash(&hash);
	qsort(keys->elem, keys->size, sizeof(void*), zTcmp);
	for (i = 0; i < keys->size; i++) {
		k = keys->elem[i];
		v = zGetHash(&hash, k);
		printf("%s => %s\n", k, v);
	}
	
	zFreeVec(keys);
	zFree(keys);
	zFreeHash(&hash);
	zFreeTVec(&kvec);
	zFreeTVec(&vvec);
}

void testIntern (FILE *stream) {
	char *s1 = (char*)"foo";
	char *s2 = (char*)"bar";
	char *r1 = (char*)"foo";
	char *r2 = (char*)"bar";
	char *a1, *a2;

	char *val;
	int   i;
	
	/*** Intern ***/
	assert(strcmp(s1,r1) == 0); assert(strcmp(s2,r2) == 0);
	
	a1 = malloc(strlen(s1)+1); 
	assert((a1 != s1) && (a1 != r1));
	a2 = malloc(strlen(s2)+1); 
	assert((a2 != s2) && (a2 != r2));

	strcpy(a1, s1); assert(strcmp(a1, r1) == 0);
	strcpy(a2, s2); assert(strcmp(a2, r2) == 0);	
	
	s1 = zStrIntern(s1, 0); assert(strcmp(s1, r1) == 0);
	s2 = zStrIntern(s2, 0); assert(strcmp(s2, r2) == 0);
	
	assert(s1 != s2); assert(strcmp(s1, s2));
	
	r1 = zStrIntern(r1, 0); assert(r1 == s1); assert(strcmp(r1, a1) == 0);
	r2 = zStrIntern(r2, 0); assert(r2 == s2); assert(strcmp(r2, a2) == 0);

	assert(r1 != r2); assert(strcmp(r1, r2));
	
	a1 = zStrIntern(a1, 1); assert(a1 == r1);
	a2 = zStrIntern(a2, 1); assert(a2 == r2);
	
	assert(strcmp(r1, a1) == 0);
	assert(strcmp(r2, a2) == 0);
	
	printf("%s %s %s\n", s1, r1, a1);
	printf("%s %s %s\n", s2, r2, a2);
	
	val = malloc(65*sizeof(char));

	/*** String Pool ***/
	while (fscanf(stream, "%64s", val) == 1) {
		i = zChar2StrIdx(val);
		assert(0 == strcmp(zStrIdx2Char(i), val));
		assert(zChar2StrIdx(val) == i);
	}

	free(a1); free(a2);
	free(val);
}
