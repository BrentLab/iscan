/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/*****************************************************************************\
modeltest.c

model test program for zoe
\*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ZOE.h"

int main (int argc, char *argv[]) {

	/* commandline */
	char *model_file;
	
	/* sequence variables */
	FILE      *stream;
	zModel     model;
	
	/* set the program name */
	zSetProgramName(argv[0]);
	
	/* process commandline */
	if (argc < 2) {
		(void)fprintf(stderr, "usage: %s <model>\n", argv[0]);
		exit(1);
	}
	model_file = argv[1];
	
	/* get model */
	if ((stream = fopen(model_file, "r")) == NULL)
		zDie("file error (%s)", model_file);
	
	/* do stuff */
	if (!zReadModel(stream, &model))
		zDie("%s could not open file %s", argv[0], model_file);
	
	zWriteModel(stdout, &model);
	zAmbiguateModel(&model, -1);
	zWriteModel(stdout, &model);
	zDeambiguateModel(&model);
	zWriteModel(stdout, &model);
	
	/* clean up */
	zFreeModel(&model);
	return 0;
}
