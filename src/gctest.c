/************************************************************************\
 gctest.c  - tests the GC functions of zDNA.h

\************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "ZOE.h"

int main(int argc, char** argv) {
    zFastaFile ff;
    zDNA       dna;
    char*      ff_filename;
    FILE*      stream;
    size_t     win;
    coor_t     i;
    float*     gcs;

    zSetProgramName(argv[0]);
    zParseOptions(&argc, argv);

    if (argc != 2) {
	printf("Usage: gctest [-w window_size] fasta_file\n");
	exit(1);
    }
    ff_filename = argv[1];
    if (zOption("w")) {
	win = atoi(zOption("w"));
    } else {
	win = 10000;
    }

    if ((stream = fopen(ff_filename, "r")) == NULL) {
	fprintf(stderr, "Couldn't open Fasta File (%s)\n", ff_filename);
	exit(2);
    }

    if (!zReadFastaFile(stream, &ff)) {
	fprintf(stderr, "Error reading fasta file (%s)\n", ff_filename);
	exit(3);
    }
    
    zFastaToDNA(&ff, &dna);

    gcs = zCalcSlidingGCLevel(&dna, win);

    for (i = 0; i < dna.length; i++) {
	printf("%d:\t%g\n", i, gcs[i]);
    }

    return 0;
}
