/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/*****************************************************************************\
conseqscannertest.c

test program for conseqscanners
\*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ZOE.h"

int main (int argc, char *argv[]) {
  
  /*  	 commandline */ 
  char *fasta_file;
  char *model_file;
  size_t begin;
  size_t end;
  
  /*	 sequence variables */
  FILE      *stream;
  zFastaFile fasta;
  zConseq    conseq, exon;
  zModel     model;
  zConseqScanner   conseqscanner;
  zSfeature  f;
  score_t    score;
  coor_t     i, length;
  char       ss[16];
  
  /*  	 set the program name */ 
  zSetProgramName(argv[0]);
  
  /*  	 process commandline */ 
  if (argc < 5) {
    (void)fprintf(stderr, "usage: %s <model> <fasta> <begin> <end>\n", argv[0]);
    exit(1);
  }
  model_file = argv[1];
  fasta_file = argv[2];
  begin = atoi(argv[3]);
  end   = atoi(argv[4]);
  
  /*   use zero-based coordinates internally */ 
  begin--;
  end--;
  
  /*  	 conseq */ 
  if ((stream = fopen(fasta_file, "r") ) == NULL)
    zDie("file error (%s)", fasta_file);
  if (!zReadFastaFile(stream, &fasta)) zDie("%s fasta error", argv[0]);
  (void)fclose(stream);
  zFastaToConseq(&fasta, &conseq, 3);
  zFreeFastaFile(&fasta);
  
  /*	 model & conseqscanner */ 
  if ((stream = fopen(model_file, "r")) == NULL)
    zDie("file open (%s)", model_file);
  if (!zReadModel(stream, &model)) zDie("%s zReadModel failed", argv[0]);
  /*	zAmbiguateModel(&model, -1); */
  
  zWriteModel (stdout, &model);
  
  zInitConseqScanner(&conseqscanner, &conseq, &model);
  
  /*	 score sequence */
  if (model.type == CDS) {
    length = end - begin + 1;
    zSubseqConseq(&conseq, &exon, begin, length - (length % 3));
    zWriteFastaFile(stdout, (zFastaFile*)&exon);
    f.name  = zChar2StrIdx("Exon");
    f.start = begin;
    f.end   = end;
    for (i = 0; i < 3; i++) {
      f.lfrag = i;
      score = conseqscanner.scoref(&conseqscanner, &f);
      zScore2Text(score, ss);
      printf("left_frag = %d, score = %s\n", i, ss);
    }
    
  } else {
    for (i = begin; i <= end; i++) {
      score = conseqscanner.score(&conseqscanner, i);
      zScore2Text(score, ss);
      (void)printf("%c pos = %d, score = %s\n", conseq.seq[i], i, ss);
    }
  }
  
  /*  	 clean up */ 
  zFreeConseqScanner(&conseqscanner);
  zFreeConseq(&conseq);
  zFreeModel(&model);
  return(0);
  
}

