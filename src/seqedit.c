/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/*****************************************************************************\
seqtest.c

tests sequence stuff
\*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ZOE.h"

void usage (void);

int main (int argc, char *argv[]) {

	/* commandline options */
	int ANTI          = 0;
	int REVERSE       = 0;
	int COMPLEMENT    = 0;
	int SUBSEQ        = 0;
	int SELECT        = 0;
	int SUBSEQ_FROM   = 0;
	int SUBSEQ_LENGTH = 0;
	int TRANSLATE     = 0;
	char *filename;
	
	/* sequence variables */
	FILE      *stream;
	zFastaFile fasta;
	zDNA       dna;
	zDNA       sub;
	zProtein   pro;
	
	/* set the program name */
	zSetProgramName(argv[0]);
			
	/* process commandline */
	if (argc == 1) usage();
	zParseOptions(&argc, argv);
	filename = argv[1];
	if (zOption("r")) REVERSE     = 1;
	if (zOption("c")) COMPLEMENT = 1;
	if (zOption("a")) ANTI       = 1;
	if (zOption("x")) TRANSLATE  = 1;
	if (zOption("S")) SELECT     = 1;
	if (zOption("s")) SUBSEQ_FROM   = atoi(zOption("s"));
	if (zOption("l")) SUBSEQ_LENGTH = atoi(zOption("l"));
	
	if (SUBSEQ_FROM && SUBSEQ_LENGTH) {
		SUBSEQ = 1;
		SUBSEQ_FROM--; /* sequence is zero-based internally */
	} else if (SUBSEQ_FROM || SUBSEQ_LENGTH) {
		(void)fprintf(stderr, "You must select both s and l for subsequence\n");
		exit(2);
	} else {
		SUBSEQ = 0;
	}
		
	/* file or stdin? */
	if ((strcmp(filename, "-") == 0)) stream = stdin;
	else if ((stream = fopen(filename, "r") ) == NULL)
			zDie("file error (%s)", filename);
			
	zReadFastaFile(stream, &fasta);
	zFastaToDNA(&fasta, &dna);
		
	/* perform sequence operations */
	if (REVERSE)    zReverseDNA(&dna);	
	if (COMPLEMENT) zComplementDNA(&dna);
	if (ANTI)       zAntiDNA(&dna);
	if (SUBSEQ) {
		if (SELECT) {
			zSelectDNA(&dna, &sub, SUBSEQ_FROM, SUBSEQ_LENGTH);
			zFreeDNA(&dna);
			dna = sub;
		} else {
			zSubseqDNA(&dna, &sub, SUBSEQ_FROM, SUBSEQ_LENGTH);
			dna = sub; /* can't free parent with subseq - OS cleanup */
		}
	}
	if (TRANSLATE) zTranslateDNA(&dna, &pro, 0);
	
	/* output */
	if (TRANSLATE) zWriteFastaFile(stdout, (zFastaFile *) &pro);
	else           zWriteFastaFile(stdout, (zFastaFile *) &dna);
	
	exit(0);
}

void usage (void) {
	(void)fprintf(stderr, "seqmunge - a tool for manipulating DNA sequences\n\n");
	(void)fprintf(stderr, "usage: seqmunge [options] <fasta file>\n");
	(void)fprintf(stderr, "options:\n");
	(void)fprintf(stderr, "  -r         reverse\n");
	(void)fprintf(stderr, "  -c         complement\n");
	(void)fprintf(stderr, "  -a         anti (reverse-complement)\n");
	(void)fprintf(stderr, "  -s <int>   subseq from i\n");
	(void)fprintf(stderr, "  -l <int>   subseq/select length\n");
	(void)fprintf(stderr, "  -S         use select instead of subseq\n");
	(void)fprintf(stderr, "  -x         translate DNA to protein\n");
	(void)fprintf(stderr, "  -          in place of file name for stdin\n");
	exit(1);
}
