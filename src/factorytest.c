/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/*****************************************************************************\
testexonfactory.c

yet another test program for zoe
\*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ZOE.h"

/* declarations for functions internal to this program */
void getDNA(zDNA*, const char*);
void getModel(zModel*, const char*);

int main (int argc, char *argv[]) {

	/* commandline */
	char *fasta_file;
	
	/* objects */
	zDNA            dna;
	zFeatureFactory efactory[4];
	zSFList        *exons;
	zSfeature      *exon;
	zModel          cds_model,   sig_model,
	                accpt_model, donor_model,
		start_model, stop_model;
/* 		            polya_model, pbcap_model, */
/* 		            tata_model; */
	zScanner      **scanner;
	int 			feature_count;
	
	/* iterators */
	coor_t i;
	int    j;
	
	/* set the program name */
	zSetProgramName(argv[0]);
	
	/* process commandline */
	if (argc < 2) {
		(void)fprintf(stderr, "%s - an exon factory test program\n", argv[0]);
		(void)fprintf(stderr, "usage: %s <fasta file>\n", argv[0]);
		exit(1);
	}
	fasta_file = argv[1];	
	
	/* parse files */
	getDNA(&dna, fasta_file);
	getModel(&accpt_model, "T/model.acceptor");
	getModel(&donor_model, "T/model.donor");
	getModel(&start_model, "T/model.start");
	getModel(&stop_model,  "T/model.stop");
	getModel(&cds_model,   "T/model.coding");
	getModel(&sig_model,   "T/model.sigpep");
/* 	getModel(&polya_model, "T/model.polya"); */
/* 	getModel(&pbcap_model, "T/model.pbcap"); */
/* 	getModel(&tata_model, "T/model.tata"); */
		
	/* create scanners */
	feature_count = zStringPoolCount();
	scanner = zCalloc(feature_count, sizeof(zScanner*), "bar");
	i = zChar2StrIdx("Coding");
	scanner[i] = zMalloc(sizeof(zScanner), "foo");
	zInitScanner(scanner[i],   &dna, NULL, &cds_model);
	i = zChar2StrIdx("Acceptor");
	scanner[i] = zMalloc(sizeof(zScanner), "foo");
	zInitScanner(scanner[i], &dna, NULL, &accpt_model);
	i = zChar2StrIdx("Donor");
	scanner[i] = zMalloc(sizeof(zScanner), "foo");
	zInitScanner(scanner[i], &dna, NULL, &donor_model);
	i = zChar2StrIdx("Start");
	scanner[i] = zMalloc(sizeof(zScanner), "foo");
	zInitScanner(scanner[i], &dna, NULL, &start_model);
	i = zChar2StrIdx("Stop");
	scanner[i] = zMalloc(sizeof(zScanner), "foo");
	zInitScanner(scanner[i],  &dna, NULL, &stop_model);
	i = zChar2StrIdx("SIG_PEP");
	scanner[i] = zMalloc(sizeof(zScanner), "foo");
	zInitScanner(scanner[i],  &dna, NULL, &sig_model);
/* 	i = zChar2StrIdx("PolyA"); */
/* 	scanner[i] = zMalloc(sizeof(zScanner), "foo"); */
/* 	zInitScanner(scanner[i],  &dna, NULL, &polya_model); */
/* 	i = zChar2StrIdx("PB_CAP"); */
/* 	scanner[i] = zMalloc(sizeof(zScanner), "foo"); */
/* 	zInitScanner(scanner[i],  &dna, NULL, &pbcap_model); */
/* 	i = zChar2StrIdx("TATA"); */
/* 	scanner[i] = zMalloc(sizeof(zScanner), "foo"); */
/* 	zInitScanner(scanner[i],  &dna, NULL, &tata_model); */
		
	/* create exon factory */
	zGetFactory("Esngl")(&efactory[0], &dna, scanner, 0, feature_count, '+', false, NULL);
	zGetFactory("Eterm")(&efactory[1], &dna, scanner, 0, feature_count, '+', false, NULL);
	zGetFactory("Einit")(&efactory[2], &dna, scanner, 0, feature_count, '+', false, NULL);
	zGetFactory("Exon" )(&efactory[3], &dna, scanner, 0, feature_count, '+', false, NULL);
			
	/* find all exons, report those longer than 100 */
	for (i = 0; i < dna.length; i++) {
		for (j = 0; j < 4; j++) {
			efactory[j].create(&efactory[j], i, exons);
			zSFListMoveFirst(exons);
			exon = zSFListGetNext(exons);
			while(exon != 0){
				if ((exon->end - exon->start + 1) > 100) {
					zWriteSfeature(stdout, exon);
				}
				exon = zSFListGetNext(exons);
			}
			zFreeSFList(exons);
			zFree(exons);
		}
	}
	
	/* clean up */
	for (j = 0; j < feature_count; j++) {
		if (NULL == scanner[j]) continue;
		zFreeScanner(scanner[j]);
	}
	zFree(scanner);
	zFreeModel(&cds_model);
	zFreeModel(&accpt_model);
	zFreeModel(&donor_model);
	zFreeModel(&start_model);
	zFreeModel(&stop_model);
	for (i = 0; i < 4; i++) {
		zFreeFeatureFactory(&efactory[i]);
	}
	zFreeDNA(&dna);
	
	exit(0);
}

/******************************************************************************\
	Functions for this program
\******************************************************************************/

void getDNA (zDNA *dna, const char *file) {
	FILE*      stream;
	zFastaFile fasta;
	if ((stream = fopen(file, "r")) == NULL) zDie("file error (%s)", file);
	if (!zReadFastaFile(stream, &fasta)) zDie("fasta parse error");
	fclose(stream);
	zFastaToDNA(&fasta, dna);
	zFreeFastaFile(&fasta);
}

void getModel (zModel *model, const char *file) {
	FILE* stream;
	if ((stream = fopen(file, "r")) == NULL) zDie("file error (%s)", file);
	if (!zReadModel(stream, model)) zDie("model parse error");
	fclose(stream);
	zChar2StrIdx(model->name);
	zAmbiguateModel(model, -1);
}


