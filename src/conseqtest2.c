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
	zConseq    conseq;
	int i;

	if (2 != argc) zDie("wrong number of parameters");
	stream = fopen(argv[1], "r");
	zReadFastaFile (stream, &fasta);

	zFastaToConseq(&fasta, &conseq, 8);

	zWriteFastaFile (stdout, (zFastaFile*)&conseq);

	for(i=0;i<conseq.conseq_length;i++) {
	  printf("%i\n", conseq.s10[i]);
	}

	printf("%i\t%i\n", conseq.bits, conseq.digits);

	return(0);
}

