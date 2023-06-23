/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/*****************************************************************************\
scannertest.c

test program for scanners
\*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ZOE.h"

int main (int argc, char *argv[]) {

	/* commandline */
	char *fasta_file;
	char *model_file;
	size_t begin;
	size_t end;
	
	/* sequence variables */
	FILE      *stream;
	zFastaFile fasta;
	zDNA       dna, exon;
	zProtein   pro;
	zModel     model;
	zScanner   scanner;
	zSfeature   f;
	score_t    score;
	coor_t     i, length;
	char       ss[16];
	
	/* set the program name */
	zSetProgramName(argv[0]);
		
	/* process commandline */
	if (argc < 5) {
		(void)fprintf(stderr, "usage: %s <model> <fasta> <begin> <end>\n", argv[0]);
		exit(1);
	}
	model_file = argv[1];
	fasta_file = argv[2];
	begin = atoi(argv[3]);
	end   = atoi(argv[4]);
	
	/* use zero-based coordinates internally */
	begin--;
	end--;
		
	/* dna */
	if ((stream = fopen(fasta_file, "r") ) == NULL)
			zDie("file error (%s)", fasta_file);
	if (!zReadFastaFile(stream, &fasta)) zDie("%s fasta error", argv[0]);
	(void)fclose(stream);
	zFastaToDNA(&fasta, &dna);
	zFreeFastaFile(&fasta);
	
	/* model & scanner */
	if ((stream = fopen(model_file, "r")) == NULL)
		zDie("file open (%s)", model_file);
	if (!zReadModel(stream, &model)) zDie("%s zReadModel failed", argv[0]);
	zAmbiguateModel(&model, -1);
	zInitScanner(&scanner, &dna, NULL, &model);
	
	/* score sequence */
	if (model.type == CDS) {
		length = end - begin + 1;
		zSubseqDNA(&dna, &exon, begin, length - (length % 3));
		zTranslateDNA(&exon, &pro, 0);
		zWriteFastaFile(stdout, (zFastaFile*)&exon);
		zWriteFastaFile(stdout, (zFastaFile*)&pro);
		f.name  = zChar2StrIdx("Exon");
		f.start = begin;
		f.end   = end;
		for (i = 0; i < 3; i++) {
			f.lfrag = i;
			score = scanner.scoref(&scanner, &f);
			zScore2Text(score, ss);
			printf("left_frag = %d, score = %s\n", i, ss);
		}
		
	} else {
		for (i = begin; i <= end; i++) {
			score = scanner.score(&scanner, i);
			zScore2Text(score, ss);
			(void)printf("%c pos = %d, score = %s\n", dna.seq[i], i, ss);
		}
	}
	
	/* clean up */
	zFreeScanner(&scanner);
	zFreeDNA(&dna);
	zFreeModel(&model);
	return(0);
}

