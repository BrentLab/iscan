/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */  
/*****************************************************************************\
zoe2gtf2.c

Convert GTF file to a list of zoe features;
\*****************************************************************************/

#include "ZOE.h"

static const char* usage =
    "Usage: \n\
     zoe2gtf hmm_file dna_file gtf_file\n\
The dna is required to determine the phases  1T, 2TA, & 2TG\n";

void openHMM(zHMM *hmm, const char* hmm_file, bool est_para_mode_bool) {
	FILE *hmm_stream;

	if ((hmm_stream = fopen(hmm_file, "r")) == NULL)
		zDie("can't open hmm file (%s)", zOption("c"));
	if (!zReadHMM(hmm_stream, hmm, est_para_mode_bool)) zDie("error reading hmm");
	fclose(hmm_stream);
}

void openDNAFromFasta(zDNA* dna, char* fasta_file) {
    FILE      *fasta_stream;
    zFastaFile fasta;
    
    if ((fasta_stream = fopen(fasta_file, "r")) == NULL)
        zDie("can't open FASTA file (%s)", fasta_file);
  
    zReadFastaFile(fasta_stream, &fasta);
    fclose(fasta_stream);
    zFastaToDNA(&fasta, dna);
    zFreeFastaFile(&fasta);
}

int main(int argc, char** argv) {
    char*      hmm_file;
    char*      dna_file;
    char*      gtf_file;
    FILE      *gtf_stream;
    zHMM       hmm;
    zDNA       dna;
    zDNA       rdna;
    zGTFVec   *gtfvec;
    zSFVec     sfv;
    bool       est_para_mode_bool = false;

    zSetProgramName(argv[0]);
    zParseOptions(&argc, argv);
  
    if (argc != 4) {
        zDie("wrong number of parameters\n%s", usage);
    }

	if( zOption("-pe")) {
		est_para_mode_bool = true;
	}

    hmm_file = argv[1];
    dna_file = argv[2];
    gtf_file = argv[3];

    openHMM(&hmm, hmm_file, est_para_mode_bool);
    openDNAFromFasta(&dna, dna_file);
    zCopyDNA(&dna, &rdna);
    zAntiDNA(&rdna);
    
    if (NULL == (gtf_stream = fopen(gtf_file, "r"))) {
        zDie("Couldn't open gtf file (%s)", gtf_file);
    }
    gtfvec = zReadGTFVec(gtf_stream);
    fclose(gtf_stream);
    
    zInitSFVec(&sfv, gtfvec->size);
    zGTFVec2SFVec(&hmm, &dna, &rdna, gtfvec, &sfv);
    zWriteSFVec(stdout, &sfv);

    zFreeSFVec(&sfv);
    zFreeGTFVec(gtfvec);
    zFree(gtfvec);

    zFreeDNA(&dna);
    zFreeHMM(&hmm);

    return 0;
}
