/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */  
/*****************************************************************************\
gtftest.c

GTF feature test program for zoe
\*****************************************************************************/

#include <stdio.h>
#include "ZOE.h"

void openHMM(zHMM *hmm, const char* hmm_file) {
	FILE *hmm_stream;

	if ((hmm_stream = fopen(hmm_file, "r")) == NULL)
		zDie("can't open hmm file (%s)", zOption("c"));
	if (!zReadHMM(hmm_stream, hmm)) zDie("error reading hmm");
	fclose(hmm_stream);
}

int main(int argc, char** argv) {
	FILE *stream;
	zGTFVec *gtfvec;

    /* set the program name */
    zSetProgramName(argv[0]);
    
	/* parse options */
	zParseOptions(&argc, argv);

    /* process commandline */
    if (argc < 2) {
		(void)fprintf(stderr, "usage: %s <file of zSfeature>\n", argv[0]);
		exit(1);
    }
	
	if ((stream = fopen(argv[1], "r")) == NULL)
		zDie("file error (%s)", argv[1]);

    gtfvec = zReadGTFVec(stream);

    if (zOption("s")) {
        zSortGTFVec(gtfvec);
    }

    zWriteGTFVec(stdout, gtfvec);

    zFreeGTFVec(gtfvec);
    zFree(gtfvec);
    gtfvec = NULL;

    fclose(stream);
    return 0;
}
