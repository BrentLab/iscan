/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/*****************************************************************************\
conseq.c

tests conservation sequence stuff
\*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ZOE.h"

int main (int argc, char *argv[]) {

	FILE      *stream;
	zFastaFile fasta;
	zDNA       dna;
	zConseq    conseq;
	zHash hash;	
	zTVec text;
  	zVec* keys;
	coor_t i;
	int    j;
	int counter;
	void* val;
	int* iptr;
	char charptr[2];
	char c;


	zInitHash(&hash);
	zInitTVec(&text, 100);

	if (argc != 3) zDie("Not enough arguments");

	stream = fopen(argv[1], "r");
	zReadFastaFile (stream, &fasta);

	zFastaToDNA(&fasta, &dna);

	zWriteFastaFile (stdout, (zFastaFile*)&dna);


	for(i=0;i<dna.length;i++) {

	  c = dna.seq[i];

	  charptr[0] = c;
	  charptr[1] = '\0';
	  
	  val = zGetHash(&hash, charptr);

	  if(val == NULL) {
	    zPushTVec(&text, charptr);
	    iptr = malloc(sizeof(int));
	    *iptr = 1;
	    zSetHash(&hash, text.last, iptr);
	  }
	  else {
	    iptr = val;
	    (*iptr)++;
	  }

	}

	counter = 0;
	keys = zKeysOfHash(&hash);
	qsort(keys->elem, keys->size, sizeof(void*), zTcmp);
	for(j=0;j< keys->size; j++) {
	  printf("%s\t%d\n", (char*)keys->elem[j], *(int*)zGetHash(&hash, (char*)keys->elem[j]));
	  counter += *(int*)zGetHash(&hash, (char*)keys->elem[j]);
	}

	printf("total\t%d\n", counter);

	fclose(stream);

	zInitHash(&hash);
	zInitTVec(&text, 100);

	stream = fopen(argv[2], "r");
	zReadFastaFile (stream, &fasta);


	zFastaToConseq(&fasta, &conseq, 3);

	zWriteFastaFile (stdout, (zFastaFile*)&conseq);

	for(i=0;i<conseq.length;i++) {

	  c = conseq.seq[i];

	  charptr[0] = c;
	  charptr[1] = '\0';
	  
	  val = zGetHash(&hash, charptr);

	  if(val == NULL) {
	    zPushTVec(&text, charptr);
	    iptr = malloc(sizeof(int));
	    *iptr = 1;
	    zSetHash(&hash, text.last, iptr);
	  }
	  else {
	    iptr = val;
	    (*iptr)++;
	  }

	}

	for(i=0;i<conseq.length;i++) {
	  printf("%i\n", conseq.s10[i]);
	}

	counter = 0;

	keys = zKeysOfHash(&hash);
	qsort(keys->elem, keys->size, sizeof(void*), zTcmp);
	for(j=0;j< keys->size; j++) {
	  printf("%s\t%d\n", (char*)keys->elem[j], *(int*)zGetHash(&hash, (char*)keys->elem[j]));
	  counter += *(int*)zGetHash(&hash, (char*)keys->elem[j]);
	}
	printf("total\t%d\n", counter);

	fclose(stream);

	zFreeVec(keys);
	zFree(keys);
	zFreeHash(&hash);
	zFreeTVec(&text);

	return(0);
}

