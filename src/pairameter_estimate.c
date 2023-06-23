/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */  
/*****************************************************************************\
pairagon2estgen.c

Convert a zoe list of features to a GTF file;
\*****************************************************************************/

#include <stdio.h>
#include <assert.h>
#include "ZOE.h"
#include "zHardCoding.h"

#define NUCLEOTIDES 5

#define ABS(x) ((x<0)?x*(-1):x)

#define SAMPLE_SIZE 22000
#define STATES 17
#define DURATION_MAX ((coor_t) pow(2,24))
#define MAX_CHAR_LENGTH (64)
#define EXTRA_TRANSITIONS (7)
#define EXTRA_INITS (5)
#define EXPLICIT_MAX (30)

static const char* usage =
    "Usage: \n\n parameter_estimate pairagon_list [-m] -cdna=cdna_list -genomic=genomic_list\n";

struct zParameter {
	zHMM_StateType stype;
	zStrIdx name;
	zStrIdx type;
	zStrIdx seq_type;
	zModelType mtype;
	int   length;
	long  count; /* # of times the state was seen */
	long  init_count;
	score_t init;
	long *transition_count;
	score_t *transition;
	long  duration;
	int  *frequency;
	long **nuc_freq;
	int   scoring;
	int   increments[2];
};
typedef struct zParameter zParameter;

extern int zDumpAFVec(FILE*, zAFVec*);
extern int zIsIntronEntry(zStrIdx, strand_t);
extern int zIsIntronExit(zStrIdx, strand_t);
extern int zMapHMM(zHMM*);
extern int FORWARD, REVERSE, BOTH;

int zInitParameterEstimate(zParameter** parameter);
int zFreeParameter(zParameter *p);
int zParameterEstimate(zParameter** parameter, zAFVec* afv, zDNA* genomic, zDNA* cdna, int transform);
int zParameterEstimateState(zParameter *parameter, int state_length, coor_t genomic_start, coor_t genomic_end, coor_t cdna_start, coor_t cdna_end, zDNA *genomic, zDNA *cdna); 
int zConsolidateParameters(zParameter **parameter); 
int zMakeParameters(zParameter **parameter, zHMM *hmm); 
void zReportParameters(zParameter **parameter); 
void zKernelSmooth(int *count, float *smooth, coor_t max, int m); 

int states = STATES;
int *smap;
long *multi_nuc_freqs = NULL;

long  random_genomic_count[NUCLEOTIDES];
long  random_cdna_count[NUCLEOTIDES];
score_t random_genomic[NUCLEOTIDES];
score_t random_cdna[NUCLEOTIDES];

int alignment_count = 0;
int file_count = 0;

/* from MGC clones 
score_t random_genomic[NUCLEOTIDES] = {0.269, 0.215, 0.223, 0.293};
score_t random_cdna[NUCLEOTIDES] = {0.2525, 0.2522, 0.2563, 0.2389};
*/

int main(int argc, char** argv) {
	char      *pairagon_list;
	char      *cdna_list;
	char      *genomic_list;

	FILE      *pairagon_list_stream;
	FILE      *cdna_list_stream;
	FILE      *genomic_list_stream;

	char      pairagon_filename[256];
	char      cdna_filename[256];
	char      genomic_filename[256];

	FILE      *stream;
	zAFVec    *afv;
	zDNA      *genomic, **multi_cdna;
	zVec      *multi_cdna_vec = NULL;

	zParameter **parameter;
	zHMM       hmm;

	int cdna_entries;
	int i;
	int cdna_orientation;
	int splice_orientation;

	/* flags for hack in estimation */
	int init_hack;
	int unify_models_hack;
	int polya_hack;
	int overhang_hack;
	int splice_entry_hack;
	int splice_exit_hack;
	int intron_length_hack;

	zChar2StrIdx("parameter_estimate");
	zSetProgramName(argv[0]);
	zParseOptions(&argc, argv);

	if (argc != 2 ||
	    zOption("cdna") == NULL ||
	    zOption("genomic") == NULL)
		zDie(usage);

	pairagon_list = argv[1];
	genomic_list = zOption("genomic");
	cdna_list = zOption("cdna");

	/* Test cDNA list file */
	if ((cdna_list_stream = fopen(cdna_list, "r")) == NULL) {
		zDie("fasta file error (%s)", cdna_list);
	}
	/* Test genomic list file */
	if ((genomic_list_stream = fopen(genomic_list, "r")) == NULL) {
		zDie("fasta file error (%s)", genomic_list);
	}
	/* Test pairagon list file */
	if ((pairagon_list_stream = fopen(pairagon_list, "r")) == NULL) {
		zDie("Couldn't open pairagon list file %s", pairagon_list);
	}

	init_hack = unify_models_hack = polya_hack = overhang_hack = splice_entry_hack = splice_exit_hack = intron_length_hack = 0;
	if (zOption("-hinit") != NULL) {
		init_hack = 1;
	}
	if (zOption("-hunify_models") != NULL) {
		unify_models_hack = 1;
	}
	if (zOption("-hpolya") != NULL) {
		polya_hack = 1;
	}
	if (zOption("-hoverhang") != NULL) {
		overhang_hack = 1;
	}
	if (zOption("-hsplice_exit") != NULL) {
		splice_exit_hack = 1;
	}
	if (zOption("-hintron_length") != NULL) {
		intron_length_hack = 1;
	}

	if (zOption("-hack_usual") != NULL) {
		init_hack = unify_models_hack = polya_hack = overhang_hack = splice_exit_hack = intron_length_hack = 1;
	}

	/* Not part of the usual hack */
	if (zOption("-hsplice_entry") != NULL) {
		splice_entry_hack = 1;
	}

	parameter = (zParameter**) zMalloc(states*sizeof(zParameter*), "main: parameter");
	zInitParameterEstimate(parameter);

	fprintf(stderr, "Reading files\n");
	while (fgets(cdna_filename, 256, cdna_list_stream) != NULL &&
	       fgets(genomic_filename, 256, genomic_list_stream) != NULL &&
	       fgets(pairagon_filename, 256, pairagon_list_stream) != NULL) {

		char buffer[256];
		sscanf(cdna_filename, "%s", buffer);
		strcpy(cdna_filename, buffer);
		sscanf(genomic_filename, "%s", buffer);
		strcpy(genomic_filename, buffer);
		sscanf(pairagon_filename, "%s", buffer);
		strcpy(pairagon_filename, buffer);

		/* Read in a cDNA fasta file with multiple sequences */
		if ((stream = fopen(cdna_filename, "r")) == NULL) {
			zDie("fasta file error (%s)", cdna_filename);
		}
		fclose(stream);
		multi_cdna_vec = (zVec*) zMalloc(sizeof(zVec), "main: multi_cdna_vec");
		zInitVec(multi_cdna_vec, 2);
		cdna_entries = zLoadMultiDNAFromMultiFasta(multi_cdna_vec, cdna_filename, NULL);

		multi_cdna = (zDNA**) zMalloc(cdna_entries*sizeof(zDNA*), "main: multi_cdna");
		for (i = 0; i < cdna_entries; i++) {
			multi_cdna[i] = (zDNA*) multi_cdna_vec->elem[i];
			zSetDNAPadding(multi_cdna[i], PADDING);
		}

		/* Read in genomic fasta file */
		if ((stream = fopen(genomic_filename, "r")) == NULL) {
			zDie("fasta file error (%s)", genomic_filename);
		}
		fclose(stream);
		genomic = (zDNA*) zMalloc(sizeof(zDNA), "main: genomic");
		zInitDNA(genomic);
		zLoadDNAFromFasta(genomic, genomic_filename, NULL);
		zSetDNAPadding(genomic, PADDING);

		if ((stream = fopen(pairagon_filename, "r")) == NULL) {
			zDie("Couldn't open pairagon file %s", pairagon_filename);
		}
		while ((afv = zReadAFVec(stream, genomic, NULL, &cdna_orientation, &splice_orientation)) != NULL) {
	/* Pass NULL so that afv->elem[n].cdna->def will be set to what's on the output file */
			zDNA *cdna = zMalloc(sizeof(zDNA), "main:cdna");
			char cdna_def[256];
			zInitDNA(cdna);
			for (i = 0; i < cdna_entries; i++) {
				strcpy(cdna_def, multi_cdna[i]->def);
				zChopString(cdna_def, 21); /* Do the exact same thing done to the def when printing the state sequence out and compare */
							   /* Add 1 since there is a '>' in the header, but comparison starts after the '>' - see next line */
				zCopyDNA(multi_cdna[i], cdna);
				if (zChar2StrIdx(cdna_def + 1) == afv->elem[0].cdna_def) {
					if (cdna_orientation == REVERSE) zAntiDNA(cdna);
					break;
				}
			}

			if (cdna == NULL) {
				zDie("Matching cDNA sequence cannot be found for '%s' in %s", zStrIdx2Char(afv->elem[0].cdna_def), cdna_filename);
			}
			if (afv->size == 0) zWarn("Did not read any features");
			for (i = 0; i < afv->size; i++) {
				afv->elem[i].cdna = cdna;
				afv->elem[i].genomic = genomic;
			}
			
			/* zWriteAlignment(stdout, afv, 0); */
			if (afv->elem[0].strand == '+') {
				if (zParameterEstimate(parameter, afv, genomic, cdna, zOption("t") != NULL) == -1) {
					/*zWarn("N found in %s", pairagon_filename);*/
				}
			} else {
				zAntiDNA(cdna);
				zAntiDNA(genomic);
				zAntiAFVec(afv);
				for (i = 0; i < afv->size; i++) {
					afv->elem[i].cdna = cdna;
					afv->elem[i].genomic = genomic;
				}
				if (zParameterEstimate(parameter, afv, genomic, cdna, zOption("t") != NULL) == -1) {
					/*zWarn("N found in %s", pairagon_filename);*/
				}
				zAntiDNA(genomic);
			}

			for (i = PADDING; i < (int)cdna->length - PADDING; i++) {
				int c = (int) zGetDNAS5(cdna, i);
				if (c < NUCLEOTIDES) {
					random_cdna_count[c]++;
				}
			}

			zFreeDNA(cdna);
			zFree(cdna);
			alignment_count++;
			if (afv != NULL) zFreeAFVec(afv);
			zFree(afv);
		}
		fclose(stream);

		file_count++;
		for (i = PADDING; i < (int)genomic->length - PADDING; i++) {
			int g = (int) zGetDNAS5(genomic, i);
			if (g < NUCLEOTIDES) {
				random_genomic_count[g]++;
			}
		}
		zFreeDNA(genomic);
		zFree(genomic);
		for (i = 0; i < cdna_entries; i++) {
			zFreeDNA(multi_cdna[i]);
			zFree(multi_cdna[i]);
		}
		zFree(multi_cdna);
		zFreeVec(multi_cdna_vec);
		zFree(multi_cdna_vec);
		for (i = 0; i < cdna_entries; i++) {
		}
	}
	fclose(cdna_list_stream);
	fclose(genomic_list_stream);
	fclose(pairagon_list_stream);

	/* Hacks before making parameters */

	fprintf(stderr, "Handling pre estimation hacks\n");

	/* Hack for making the IntronU2 and IntronU12 equi-length by merging the two*/
	if (intron_length_hack == 1) {
		int intronu2  = smap[zChar2StrIdx("IntronU2")];
		int intronu12 = smap[zChar2StrIdx("IntronU12")];

		parameter[intronu2]->duration += parameter[intronu12]->duration;
		parameter[intronu2]->count    += parameter[intronu12]->count;
		parameter[intronu12]->duration = parameter[intronu2]->duration;
		parameter[intronu12]->count    = parameter[intronu2]->count;
	}
	/* End Hack for making the IntronU2 and IntronU12 equi-length */


	/* Hack for U12 --> other transition probabilities */
	if (splice_exit_hack == 1) {
		/*
		Exit2		StateA		should reflect Exit --> StateA
		*/

		int exit, exit2;
		exit   = smap[zChar2StrIdx("ExitU2")];
		exit2  = smap[zChar2StrIdx("ExitU12")];

		for (i = 0; i < states; i++) {
			parameter[exit2]->transition_count[i] = parameter[exit]->transition_count[i];
		}
	}
	/* End Hack for U12 --> other transition probabilities */

	/* Hack for other --> Intron transition probabilities */
	/* GT/AG: 98.9, GC/AG: 01.0, AT/AC: 00.1 */
	if (splice_entry_hack == 1) {
		/*
		StateA		Entry2		should be 0.1/99.9*StateA-->Entry
		*/

		int entry, entry2;
		entry  = smap[zChar2StrIdx("EntryU2")];
		entry2 = smap[zChar2StrIdx("EntryU12")];

		for (i = 0; i < states; i++) {
			if (parameter[i]->transition_count[entry] > 0) {
				parameter[i]->transition_count[entry2] = (long)ceil(0.1/99.9*parameter[i]->transition_count[entry]);
			}
		}


		/* Hack to let GC/AG introns in by 1:98.9 ratio */
		parameter[entry]->nuc_freq[0][NUCLEOTIDES*2+1] =  10; /* GC */
		parameter[entry]->nuc_freq[0][NUCLEOTIDES*2+3] = 989; /* GT */
	}
	/* End Hack for other --> Intron transition probabilities */

	/* Hack for unaligned states */
	if (overhang_hack == 1) {
		/*
		manually added transitions are:
		RGenomic1       Match           0.9
		RGenomic1       RCDna1          0.1   0.9:0.1 will try to extend the alignment more
		RCDna1          Match           1.000000
		Match           RGenomic2       90% of the sequences should do this
		Match           RCDna2          10% of the sequences should do this
		RGenomic2       RCDna2          1.000000           
		RCDna1		RGenomic1	0
		*/

		char from[EXTRA_TRANSITIONS][MAX_CHAR_LENGTH] = {"RGenomic1", "RGenomic1", "RCDna1", "Match",     "Match",  "RGenomic2", "RCDna1"};
		char to[EXTRA_TRANSITIONS][MAX_CHAR_LENGTH]   = {"Match",     "RCDna1",    "Match",  "RGenomic2", "RCDna2", "RCDna2",    "RGenomic1"};
		long extra[EXTRA_TRANSITIONS]                 = { 9,           1,           1,        0,           0,        1,           0};  
		extra[4] = alignment_count/10;         /* Match --> RCDna2 */
		extra[3] = alignment_count - extra[4]; /* Match --> RGenomic2 */

		for (i = 0; i < EXTRA_TRANSITIONS; i++) {
			int from_state = smap[zChar2StrIdx(from[i])];
			int to_state   = smap[zChar2StrIdx(to[i])];
			parameter[from_state]->transition_count[to_state] = extra[i];
		}
	}
	/* End Hack for unaligned states */

	/* Hack for avoiding polyA in match state */
	if (polya_hack == 1) {
		/* Boost the number of occurrences of A by 50% */
		int rcdna2;
		rcdna2   = smap[zChar2StrIdx("RCDna2")];
		parameter[rcdna2]->nuc_freq[0][0] += (parameter[rcdna2]->nuc_freq[0][0]/2);
		parameter[rcdna2]->scoring = 1;
	}
	/* End Hack for avoiding polyA in match state */

	fprintf(stderr, "Consolidating parameters\n");
	zConsolidateParameters(parameter); 

	fprintf(stderr, "Making parameters\n");
	zMakeParameters(parameter, &hmm); 
	if (zOption("m") != NULL) {
		for (i = 0; i < hmm.models; i++) {
			zAmbiguateModel(&hmm.model[i], -1);
			if (zIsIntronEntry(zChar2StrIdx(hmm.model[i].name), '+') || zIsIntronExit(zChar2StrIdx(hmm.model[i].name), '+')) { /* Donor/Acceptor */
			} else {
				zMarginalizeModel(&hmm.model[i], 5);
			}
		}
	}

	fprintf(stderr, "Handling post estimation hacks\n");
	/* Hack for initial probabilities */
	if (init_hack == 1) {
		/* Override estimated initial probabilities with RG1=0.985, Match=0.010, RC1 = 0.005 */

		hmm.state[zGetFivePrimeGenomic(&hmm)].init[0] = hmm.state[zGetThreePrimeGenomic(&hmm)].init[0] = zFloat2Score(0.985);
		hmm.state[zGetFivePrimeCDna(&hmm)].init[0]    = hmm.state[zGetThreePrimeCDna(&hmm)].init[0]    = zFloat2Score(0.005);
		hmm.state[zGetMatch(&hmm)].init[0]                                                             = zFloat2Score(0.010);
	}

	/* End Hack for initial probabilities */

	/* Hack for fixing the durations of overhangs the same as NULL model */
	if (unify_models_hack == 1) {
		for (i = 0; i < hmm.orig_states; i++) {
			zHMM_State *state = &hmm.orig_state[i];
			if (zIsGenomicOverhang(state->name)) {
				state->duration = zChar2StrIdx(GENOMIC_NULL_DURATION_NAME);
			} else if (zIsCDnaOverhang(state->name)) {
				state->duration = zChar2StrIdx(CDNA_NULL_DURATION_NAME);
			}
		}
	}
	/* End Hack for fixing the durations of overhangs the same as NULL model */

	fprintf(stderr, "Done");
	zMapHMM(&hmm);
	zWriteHMM(stdout, &hmm);
	zFreeHMM(&hmm);
	for (i = 0; i < states; i++) {
		zFreeParameter(parameter[i]);
		zFree(parameter[i]);
	}
	zFree(parameter);
	zFree(smap);
	zStringPoolFree();
	return 0;
}

int zFreeParameter(zParameter *p) {
	int j;
	zFree(p->transition_count);
	zFree(p->transition);
	for (j = 0; j < p->length; j++) zFree(p->nuc_freq[j]);
	zFree(p->nuc_freq);
	return 1;
}

int zInitParameterEstimate(zParameter** parameter) {
	int i, j;
	char state_names[STATES][32]       = {"RGenomic1", "RCDna1", "Match",  "EntryU2", "ExitU2", "EntryU12", "ExitU12",  "IntronU2", "IntronU12", "BranchU2", "BrAccU2", "BranchU12", "BrAccU12",  "Genomic", "CDna",   "RGenomic2", "RCDna2"};
	zHMM_StateType state_types[STATES] = {INTERNAL,    INTERNAL, INTERNAL, INTERNAL,  INTERNAL, INTERNAL,   INTERNAL,   INTERNAL,   INTERNAL,    INTERNAL,   EXPLICIT,  INTERNAL,    EXPLICIT,    INTERNAL,  INTERNAL, INTERNAL,    INTERNAL};
	int increments[STATES][2]          = {{1,0},       {0,1},    {1,1},    {8,0},     {6,0},    {8,0},      {6,0},      {1,0},      {1,0},       {8,0},      {1,0},     {8,0},       {1,0},       {1,0},     {0,1},    {1,0},       {0,1}};
	zModelType model_types[STATES]     = {LUT,         LUT,      LUT,      WMM,       WMM,      WMM,        WMM,        LUT,        LUT,         WMM,        LUT,       WMM,         LUT,         LUT,       LUT,      LUT,         LUT};
	/* Set the scoring status: 0 if the model is equal to the null model */
	int scoring[STATES]                = {0,           0,        1,        1,         1,        1,          1,          0,          0,           0,          0,         1,           0,           0,         0,        0,           0};
	for (i = 0; i < states; i++) {
		parameter[i]                = (zParameter*) zMalloc(sizeof(zParameter), "zInitParameterEstimate(): p");
		parameter[i]->name          = zChar2StrIdx(state_names[i]);
		parameter[i]->stype         = state_types[i];
		parameter[i]->mtype         = model_types[i];
		parameter[i]->increments[0] = increments[i][0];
		parameter[i]->increments[1] = increments[i][1];
		parameter[i]->scoring       = scoring[i];
	}
	/* Split it into two since I need to call zChar2StrIdx before I could use the StringPoolCount, but I need them both inside this function */
	smap = (int*) zCalloc(zStringPoolCount() + 1, sizeof(int), "main:smap");
	for (i = 0; i < states; i++) {
		zParameter *p = parameter[i];
		smap[p->name] = i;
		p->count = 0;

		/* Set model lengths */

		p->length = p->increments[0] + p->increments[1];

		if (p->stype == EXPLICIT) {
			p->frequency = (int*) zMalloc(sizeof(int)*(EXPLICIT_MAX+1), "zInitParameterEstimate: p->frequency");
			for (j = 0; j <= EXPLICIT_MAX; j++) {
				p->frequency[j] = 0;
			}
		}

		p->init_count = 0;
		p->duration = 0;
		p->transition_count = zMalloc(sizeof(long)*states, "zInitParameterEstimate(): parameter->transition_count");
		p->transition       = zMalloc(sizeof(score_t)*states, "zInitParameterEstimate(): parameter->transition");
		for (j = 0; j < states; j++) {
			p->transition_count[j] = 0;
			p->transition[j] = 0.0;
		}
		p->nuc_freq = (long**) zMalloc(p->length*sizeof(long*), "zInitParameterEstimate(): nuc_freq");
		for (j = 0; j < p->length; j++) {
			p->nuc_freq[j] = (long*) zMalloc(zPOWER[NUCLEOTIDES][ABS(p->length)]*sizeof(long), "zInitParameterEstimate(): nuc_freq[j]");
		}

		/* Initialize pseudocounts */

		for (j = 0; j < zPOWER[NUCLEOTIDES][ABS(p->length)]; j++) {
			int k;
			for (k = 0; k < p->length; k++) p->nuc_freq[k][j] = 0;
		}
	}
	return 1;
}

int zCountMultiplets(int order, zDNA* dna) {
	int j, index;
	coor_t i;
	if (multi_nuc_freqs == NULL) {
		multi_nuc_freqs = (long*) zMalloc(sizeof(long)*(zPOWER[NUCLEOTIDES][order]), "zCountMultiplets:multi_nuc_freqs");
		for (i = 0; i < (coor_t)zPOWER[NUCLEOTIDES][order]; i++) {
			multi_nuc_freqs[i] = 0;
		}
	}
	for (i = PADDING + order - 1; i < dna->length - PADDING; i++) {
		index = 0;
		for (j = 0; j < order; j++) {
			index += zPOWER[NUCLEOTIDES][j] * (int)(zGetDNAS5(dna, i-j));
		}
		multi_nuc_freqs[index]++;
	}
	return 1;
}

int zTransformAlnFeature(zAlnFeature *in, zAlnFeature *out) {
	char name[256];
	zCopyAlnFeature(in, out);
	strcpy(name, zStrIdx2Char(out->name));
	/* 2 becomes 8, and hence the +6 and -6 */
	if (strcmp(name, "Entry") == 0) {
		out->name = zChar2StrIdx("EntryU2");
		out->genomic_end += 6;
		out->length += 6;
	/* 2 becomes 6, and hence the +4 and -4 */
	} else if (strcmp(name, "Exit") == 0) {
		out->name = zChar2StrIdx("ExitU2");
		out->genomic_start -= 4;
		out->length += 4;
	} else if (strcmp(name, "IntronC") == 0) {
		out->name = zChar2StrIdx("IntronU2");
		out->genomic_start += 6;
		out->genomic_end -= 4;
		out->length -= 10;
	}
	return 1;
}

int zTransformAFVec(zAFVec *in, zAFVec *out) {
	int i;
	zAlnFeature* af = zMalloc(sizeof(zAlnFeature), "zTransformAFVec: af");
	for (i = 0; i < in->size; i++) {
		zTransformAlnFeature(&in->elem[i], af);
		zPushAFVec(out, af);
	}
	zFree(af);
	return 1;
}

int zParameterEstimate(zParameter** parameter, zAFVec* afv1, zDNA* genomic, zDNA* cdna, int transform) {
        int i;
        int previous_state = -1;
	zAFVec* afv;
	zAFVec* afv2 = zMalloc(sizeof(zAFVec), "zParameterEstimate: afv2");
	zInitAFVec(afv2, 2);
        /*fprintf(stderr, "Scanning %s vs %s\n", cdna->def, genomic->def);*/
	afv = afv1;
	if (transform == 1) {
		zTransformAFVec(afv1, afv2);
		afv = afv2;
	}
        for (i = 0; i < afv->size; i++) {
                zAlnFeature *af = &afv->elem[i];
                int state = smap[af->name];

		/* Process emission parameters */

		if (af->genomic_end > genomic->length || af->cdna_end > cdna->length) {
			zDie("Wrong sequence for alignment %s vs %s", genomic->def, cdna->def);
		}
		if (zParameterEstimateState(parameter[state], af->length, af->genomic_start, af->genomic_end, af->cdna_start, af->cdna_end, genomic, cdna) == -1) return -1;

		/* Process transition parameters */

		if (i > 0) {
			(parameter[previous_state]->transition_count[state])++;
		} else { /* Init Prob */
			parameter[state]->init_count++;
		}
		previous_state = state;
        }
	zFreeAFVec(afv2);
	zFree(afv2);
        return 1;
}

int zParameterEstimateState(zParameter *parameter, int state_length, coor_t genomic_start, coor_t genomic_end, coor_t cdna_start, coor_t cdna_end, zDNA *genomic, zDNA *cdna) {
	int length = parameter->length;
	coor_t i, j;
	zStrIdx name = parameter->name;
	zModelType mtype = parameter->mtype;

	if (parameter->stype == EXPLICIT && state_length <= EXPLICIT_MAX) {
		parameter->frequency[state_length]++;
	}
	parameter->duration += state_length;
	parameter->count++;
	if (zIsCDnaOnly(name) && mtype == LUT) { /* All CDna states are length 1. If that changes, so will this! */
		for (i = cdna_start + 1; i <= cdna_end; i ++ ) {
			int c = (int) zGetDNAS5(cdna, i);
			if (c < NUCLEOTIDES) parameter->nuc_freq[0][c]++;
		}
	} else if (zIsGenomicOnly(name) && mtype == LUT) { /* Genomic stuff, Entry, Exit, Introns */
			for (i = genomic_start + 1; i <= genomic_end; i += length) {
				int index = 0;
/*
if (length == 2) {
printf("Found %c%c at %u in %s for %s vs %s\n", zGetDNASeq(genomic, i), zGetDNASeq(genomic, i+1), i, zStrIdx2Char(name), genomic->def, cdna->def);
}
*/
				for (j = 0; j < (coor_t) length; j++) {
					int g = (int) zGetDNAS5(genomic, i+j);
					if (g < 4) {
						index += g*zPOWER[NUCLEOTIDES][length - 1 - j];
					} else if (length == 2) {
						zDie("N found in Splice sites at %u: %s", genomic->def, i+j);
						/*return -1;*/
					} /* N in any state is skipped */
				}
				parameter->nuc_freq[0][index]++;
				/* 5*cdna_i + cdna_i+1 will give the single-D location of the 5 X 5 matrix */
			}
	} else if (zIsGenomicOnly(name) && mtype == WMM) { /* WMM */
			for (j = 0; j < (coor_t) length; j++) {
				int index = 0;
				int g = (int) zGetDNAS5(genomic, genomic_start+1+j);
				if (g < 4) {
					index = g;
				} /* N in any state is skipped */
				parameter->nuc_freq[j][index]++;
			}
	} else if (zIsMatch(name) && mtype == LUT) {
		switch(length) {
			case 2: /* Match */
				for (i = cdna_start + 1, j = genomic_start + 1; i <= cdna_end && j <= genomic_end; i++, j++) {
					int c = (int) zGetDNAS5(cdna, i);
					int g = (int) zGetDNAS5(genomic, j);
					if (g < NUCLEOTIDES && c < NUCLEOTIDES) {
						parameter->nuc_freq[0][NUCLEOTIDES*g + c]++;
					} else {
						/*zWarn("N found in Match: %s", genomic->def);*/
						/*return -1;*/
					}
					/* 5*genomic + cdna will give the single-D location of the 5 X 5 matrix */
				}
				break;
			default:
				zDie("Model length %u not found", length);
				break;
		}
	}
	return 1;
}

void zReportParameters(zParameter **parameter) {
	int i, j;
	for (i = 0; i < states; i++) {
		zParameter *p = parameter[i];
		long        length;
		fprintf(stdout, "%s\t%ld\t%.8f\t%ld\n", zStrIdx2Char(p->name), p->count, p->init, p->duration);
		for (j = 0; j < states; j++) {
			fprintf(stdout, "   \t%s\t%s\t%ld\n", zStrIdx2Char(p->name), zStrIdx2Char(parameter[j]->name), p->transition_count[j]);
		}
		if (p->mtype == WMM) {
			length = ABS(p->length);
			for (j = 0; j < length; j++) {
				int k;
				for (k = 0; k < NUCLEOTIDES; k++) {
					fprintf(stdout, "   \t%ld", p->nuc_freq[j][k]);
				}
				fprintf(stdout, "\n");
			}
		} else { /* LUT */
			length = zPOWER[NUCLEOTIDES][ABS(p->length)];
			for (j = 0; j < length; j++) {
				fprintf(stdout, "   \t%ld", p->nuc_freq[0][j]);
				if ((j+1)%NUCLEOTIDES == 0) {
					fprintf(stdout, "\n");
				}
			}
		}
	}
}

int zConsolidateParameters(zParameter **parameter) {
	int i, j;
	long init_count = 0, count = 0;
	for (i = 0; i < states; i++) {
		zParameter *p = parameter[i];
		init_count += p->init_count;
	}
		
	for (i = 0; i < states; i++) {
		zParameter *p = parameter[i];
		long transition_count = 0;
		p->init = 1.0 * p->init_count / init_count;
		for (j = 0; j < states; j++) {
			transition_count += p->transition_count[j];
		}
		for (j = 0; j < states; j++) {
			p->transition[j] = 1.0 * p->transition_count[j] / transition_count;
		}
	}

	count = 0;
	for (j = 0; j < NUCLEOTIDES; j++) 
		count += random_genomic_count[j];
	for (j = 0; j < NUCLEOTIDES; j++) 
		random_genomic[j] = 1.0*random_genomic_count[j]/count;

	count = 0;
	for (j = 0; j < NUCLEOTIDES; j++) 
		count += random_cdna_count[j];
	for (j = 0; j < NUCLEOTIDES; j++) 
		random_cdna[j] = 1.0*random_cdna_count[j]/count;

	return 1;
}

int zMakeParameters(zParameter **parameter, zHMM *hmm) {

	int i, j, count;
	long random_count;

	
        /* Ininilize HMM object's pointers */

	hmm->name       = NULL;
	hmm->orig_state = NULL;                                                                  
	hmm->state      = NULL;                                                                  
	hmm->transition = NULL; 
	hmm->duration   = NULL;                                                                  
	hmm->model      = NULL;
	hmm->est_model  = NULL;
	hmm->est_models = 0;
	hmm->cons_model = NULL;
	hmm->cons_models = 0;
	hmm->phylo_model = NULL;
	hmm->phylo_models = 0;
	hmm->gtf_conv = NULL;
	hmm->feature_count = 0;
	hmm->dmap  = NULL;
	hmm->mmap  = NULL;
	hmm->somap = NULL;                                                                       
	hmm->simap = NULL;                                                                       
	hmm->tmap = NULL;                                                                        
	hmm->jmap = NULL;
	hmm->fmap = NULL;                                                                        

	hmm->name = (char*) zMalloc(sizeof(char)*MAX_CHAR_LENGTH, "zMakeParameters: hmm->name");
	hmm->name = strcpy(hmm->name, "pairagon");

	/* States <STATES> */

	fprintf(stderr, "\tMaking states\n");
	hmm->states = states;
	hmm->orig_states = states;
	hmm->orig_state  = (zHMM_State*) zMalloc(sizeof(zHMM_State)*hmm->orig_states, "zMakeParameters: hmm->state");
	hmm->iso_states = 1;
	hmm->iso_state = (float*) zMalloc(sizeof(float)*hmm->iso_states, "zMakeParameters: hmm->iso_state");
	hmm->iso_state[0] = 100;
	for (i = 0; i < hmm->orig_states; i++) {
		zParameter *p = parameter[i];
		zHMM_State *state = &hmm->orig_state[i];
		state->type = p->stype;
		state->name = p->name;
		state->strand = '+';
		state->phase  = 0;
		state->ffactory = NULL; 
		state->init = (score_t*) zMalloc(sizeof(score_t)*hmm->iso_states, "zMakeParameters: state->init");
		state->init[0] = zFloat2Score(p->init);
		state->duration = p->name;
		state->model = p->name;

	}

	hmm->states = states;
	hmm->state  = (zHMM_State*) zMalloc(sizeof(zHMM_State)*hmm->orig_states, "zMakeParameters: hmm->state");

	hmm->simap      = zMalloc(states * sizeof(zIVec), "zReadHMM simap");
	for (i = 0; i < hmm->orig_states; i++) {
		memcpy(&hmm->state[i], &hmm->orig_state[i], sizeof(zHMM_State));
		zInitIVec(&hmm->simap[i], 2);
		zPushIVec(&hmm->simap[i], i);
	}

	/* Transitions <HMM_TRANSITIONS> */

	fprintf(stderr, "\tMaking transitions\n");
	hmm->transitions = hmm->states * hmm->states;
	hmm->transition  = (zTransition*) zMalloc(sizeof(zTransition)*hmm->transitions, "zMakeParameters: hmm->transition");
	hmm->iso_transitions = 1;
	hmm->iso_transition = (float*) zMalloc(sizeof(float)*hmm->iso_transitions, "zMakeParameters: hmm->iso_transition");
	hmm->iso_transition[0] = 100;
	count = 0;
	for (i = 0; i < hmm->orig_states; i++) {
		for (j = 0; j < hmm->orig_states; j++) {
			if (parameter[i]->transition_count[j] != 0) {
				zTransition *transition = &hmm->transition[count];
				transition->from = parameter[i]->name;
				transition->to   = parameter[j]->name;
				transition->prob = zMalloc(sizeof(score_t)*hmm->iso_transitions, "zMakeParameters: transition->prob");
				transition->prob[0] = parameter[i]->transition[j]; 
				transition->score = zMalloc(sizeof(score_t)*hmm->iso_transitions, "zMakeParameters: transition->score");
				transition->score[0] = zFloat2Score(parameter[i]->transition[j]); 
				count++;
			}
		}
	}
	hmm->transitions = count;

	/* Durations <STATE_DURATIONS> */

	fprintf(stderr, "\tMaking durations\n");
	hmm->durations = hmm->orig_states + 2;
	hmm->duration  = (zDurationGroup*) zMalloc(sizeof(zDurationGroup)*hmm->durations, "zMakeParameters: hmm->duration");
	for (i = 0; i < hmm->orig_states; i++) {
		zParameter *p = parameter[i];
		zDurationGroup *dg = &hmm->duration[i];
		zDuration *dur;
		zDistribution *d;

		dg->name = p->name;
		dg->durations = 1;
		dg->duration  = (zDuration*) zMalloc(sizeof(zDuration)*dg->durations, "zMakeParameters: dg->duration");
		dg->iso_bound = (float*) zMalloc(sizeof(float)*dg->durations, "zMakeParameters: dg->iso_bound");
		dg->iso_bound[0] = 100;

		dur = &dg->duration[0];
		dur->min = 1;
		dur->max = DURATION_MAX;

		if (hmm->state[i].type == EXPLICIT) { /* explicit comparison to p's member. fix it soon */
			float *smooth;
			dur->distributions = 2;
			dur->distribution = (zDistribution*) zMalloc(sizeof(zDistribution)*dur->distributions, "zMakeParameters: dur->distribution");
			d = &dur->distribution[0];
			d->type = DEFINED;
			d->start = 1;
			d->end   = EXPLICIT_MAX;
			d->params = EXPLICIT_MAX;
			smooth = (float*) zMalloc(sizeof(float)*d->params, "zMakeParameters: smooth");
			zKernelSmooth(p->frequency, smooth, d->end, 8);
			d->param = smooth;
			d = &dur->distribution[1];
			d->type = CONSTANT;
			d->start = EXPLICIT_MAX+1;
			d->end   = DURATION_MAX;
			d->params = 1;
			d->param = (float*) zMalloc(sizeof(float)*d->params, "zMakeParameters: d->param");
			d->param[0] = -300;
		} else if (p->mtype == WMM) {
			dur->distributions = 2;
			dur->distribution = (zDistribution*) zMalloc(sizeof(zDistribution)*dur->distributions, "zMakeParameters: dur->distribution");
			d = &dur->distribution[0];
			d->type = CONSTANT;
			d->start = 1;
			d->end   = 1;
			d->params = 1;
			d->param = (float*) zMalloc(sizeof(float)*d->params, "zMakeParameters: d->param");
			d->param[0] = 0;
			d = &dur->distribution[1];
			d->type = CONSTANT;
			d->start = 2;
			d->end   = DURATION_MAX;
			d->params = 1;
			d->param = (float*) zMalloc(sizeof(float)*d->params, "zMakeParameters: d->param");
			d->param[0] = -300;
		} else {
			dur->distributions = 1;
			dur->distribution = (zDistribution*) zMalloc(sizeof(zDistribution)*dur->distributions, "zMakeParameters: dur->distribution");
			d = &dur->distribution[0];
			d->type = GEOMETRIC;
			d->start = 1;
			d->end   = DURATION_MAX;
			d->params = 2;
			d->param = (float*) zMalloc(sizeof(float)*d->params, "zMakeParameters: d->param");
			d->param[1] = 1.;
			if (p->count == 0) {
				d->param[0] = 0;
			} else {
				d->param[0] = (int)ceil(1.0 * p->duration / p->count);
			}
		}
	}
	for (i = hmm->orig_states; i < hmm->orig_states + 2; i++) {
		zDurationGroup *dg = &hmm->duration[i];
		zDuration *dur;
		zDistribution *d;

		dg->durations = 1;
		dg->duration  = (zDuration*) zMalloc(sizeof(zDuration)*dg->durations, "zMakeParameters: dg->duration");
		dg->iso_bound = (float*) zMalloc(sizeof(float)*dg->durations, "zMakeParameters: dg->iso_bound");
		dg->iso_bound[0] = 100;

		dur = &dg->duration[0];
		dur->min = 1;
		dur->max = DURATION_MAX;
		dur->distributions = 1;
		dur->distribution = (zDistribution*) zMalloc(sizeof(zDistribution)*dur->distributions, "zMakeParameters: dur->distribution");
		d = &dur->distribution[0];
		d->type = GEOMETRIC;
		d->start = 1;
		d->end   = DURATION_MAX;
		d->params = 2;
		d->param = (float*) zMalloc(sizeof(float)*d->params, "zMakeParameters: d->param");
	}

	random_count = 0;
	for (i = 0; i < NUCLEOTIDES; i++) random_count += random_genomic_count[i];
	hmm->duration[hmm->orig_states].name = zChar2StrIdx(GENOMIC_NULL_DURATION_NAME);
	hmm->duration[hmm->orig_states].duration[0].distribution[0].param[0] = random_count/file_count;
	hmm->duration[hmm->orig_states].duration[0].distribution[0].param[1] = 1.;

	random_count = 0;
	for (i = 0; i < NUCLEOTIDES; i++) random_count += random_cdna_count[i];
	hmm->duration[hmm->orig_states + 1].name = zChar2StrIdx(CDNA_NULL_DURATION_NAME);
	hmm->duration[hmm->orig_states + 1].duration[0].distribution[0].param[0] = random_count/alignment_count;
	hmm->duration[hmm->orig_states + 1].duration[0].distribution[0].param[1] = 1.;

	/* Sequence Models <SEQUENCE_MODELS> */

	fprintf(stderr, "\tMaking models\n");
	hmm->models = 3*hmm->orig_states;
	hmm->model  = (zModel*) zMalloc(sizeof(zModel)*hmm->models, "zMakeParameters: hmm->model");
	
	count = 0;
	for (i = 0; i < hmm->orig_states; i++) {
		zParameter *p = parameter[i];
		zModel *model = &hmm->model[2*i];
		zModel *null  = &hmm->model[2*i+1];

		long nuc_freq_sum = 0;

		char name[MAX_CHAR_LENGTH];
		int mem;

		model->name = (char*) zMalloc(sizeof(char)*MAX_CHAR_LENGTH, "zMakeParameters: model->name");
		strcpy(model->name, zStrIdx2Char(p->name));

		null->name = (char*) zMalloc(sizeof(char)*MAX_CHAR_LENGTH, "zMakeParameters: null->name");
		strcpy(name, zStrIdx2Char(p->name));
		strcat(name, NULL_MODEL_SUFFIX);
		strcpy(null->name, name);
		zChar2StrIdx(null->name);

		model->length    = null->length    = ABS(p->length);
		model->symbols   = null->symbols   = NUCLEOTIDES;
		model->submodels = null->submodels = 0;
		model->submodel  = null->submodel  = NULL;
		model->type      = null->type      = p->mtype;

		mem = zPOWER[model->symbols][model->length];
		if (model->type == WMM) {
			mem = model->symbols*model->length;
		}
		model->data = zMalloc(mem * sizeof(score_t), "zMakeParameters: model->data");
		null->data  = zMalloc(mem * sizeof(score_t), "zMakeParameters: null->data");

		if (zIsGenomicOnly(hmm->orig_state[i].name)) {
			model->seq_type = null->seq_type = GENOMIC;
			model->focus    = null->focus    = model->length - 1;
		} else if (zIsCDnaOnly(hmm->orig_state[i].name)) {
			model->seq_type = null->seq_type = DNA;
			model->focus    = null->focus    = model->length - 1;
		} else if (zIsMatch(hmm->orig_state[i].name)) {
			model->seq_type = null->seq_type = PAIR;
			model->focus    = null->focus    = model->length/2 - 1;
		}

		/* Fill the data */

		fprintf(stderr, "\t\t%s\n", model->name);
		if (p->mtype == LUT) {
			for (j = 0; j < zPOWER[null->symbols][null->length]; j++) {
				char string[50];
				int stringlength;
				int k, lut_position;
				score_t probability;

				/* Get the number as base NUCLEOTIDES string */

				zDecToBase(j, null->symbols, string);
				stringlength = (int)strlen(string);

				/* Pad it with 0s in the left and make it model->length long */
				string[null->length]='\0';
				for (k = stringlength - 1; k >= 0; k--) {
					string[k+null->length-stringlength] = string[k];
				}
				for (k = 0; k < (int)null->length-stringlength; k++) {
					string[k] = '0';
				}
				stringlength = (int)strlen(string);

				assert(stringlength == (int)null->length);

				probability = 1;
				for (lut_position = 0; lut_position < (int)null->length; lut_position++) {
					int base = string[lut_position] - '0';
					if (null->seq_type == GENOMIC) {
						probability *= random_genomic[base];
					} else if (null->seq_type == DNA) {
						probability *= random_cdna[base];
					} else if (null->seq_type == PAIR) {
						if (lut_position < (int)null->length/2) {
							probability *= random_genomic[base];
						} else {
							probability *= random_cdna[base];
						}
					}
				}
				null->data[j] = zFloat2Score(probability);
			}
			switch (p->scoring) {
				case 0: /* non scoring states */
					memcpy(model->data, null->data, sizeof(score_t)*mem);
					break;
				case 1: /* scoring states */
					nuc_freq_sum = 0;
					for (j = 0; j < zPOWER[model->symbols][model->length]; j++) {
						nuc_freq_sum += p->nuc_freq[0][j];
					}
					for (j = 0; j < zPOWER[model->symbols][model->length]; j++) {
						if (p->nuc_freq[0][j] > 0) {
							model->data[j] = zFloat2Score(1.0 * p->nuc_freq[0][j] / nuc_freq_sum);
						} else {
							model->data[j] = MIN_SCORE;
						}
					}
					break;
				default:
					zDie("scoring should be 0 or 1");
			}
		} else if (p->mtype == WMM) {
			coor_t wmm;
			int mcount = 0;
			for (wmm = 0; wmm < null->length; wmm++) {
				for (j = 0; j < null->symbols; j++) {
					if (null->seq_type == GENOMIC) {
						null->data[mcount++] = zFloat2Score(random_genomic[j]);
					} else if (null->seq_type == DNA) {
						null->data[mcount++] = zFloat2Score(random_cdna[j]);
					} else {
						zDie("Can only estimate WMM parameters for GENOMIC and DNA types");
					}
				}
			}
			switch (p->scoring) {
				case 0: /* non scoring states */
					memcpy(model->data, null->data, sizeof(score_t)*mem);
					break;
				case 1: /* scoring states */
					mcount = 0;
					for (wmm = 0; wmm < model->length; wmm++) {
						nuc_freq_sum = 0;
						for (j = 0; j < model->symbols; j++) {
							nuc_freq_sum += p->nuc_freq[wmm][j];
						}
						for (j = 0; j < model->symbols; j++) {
							if (p->nuc_freq[wmm][j] > 0) {
								model->data[mcount++] = zFloat2Score(1.0 * p->nuc_freq[wmm][j] / nuc_freq_sum);
							} else {
								model->data[mcount++] = MIN_SCORE;
							}
						}
					}
					break;
				default:
					zDie("scoring should be 0 or 1");
			}
		} else {
			zDie("Can only estimate WMM and LUT models");
		}
	}
	count = hmm->orig_states*2;
	hmm->models = count;

	hmm->increments = zMalloc(sizeof(int*)*hmm->orig_states, "zMakeParameters:hmm->increments");
	for (i = 0; i < hmm->orig_states; i++) {
		hmm->increments[i] = zMalloc(2*sizeof(int), "zMakeParameters:hmm->increments[i]");
		hmm->increments[i][0] = parameter[i]->increments[0];
		hmm->increments[i][1] = parameter[i]->increments[1];
	}

	/* Freeze the hmm's features */
	hmm->feature_count = zStringPoolCount();

	return 1;
}

score_t zKernelFunction(score_t sigma, int distance) {
        score_t exponent = -0.5*(distance/sigma)*(distance/sigma);
        score_t denominator = sqrt(2*3.1415926)*sigma;
        return exp(exponent)/denominator;
}

void zKernelSmooth(int *count, float *smooth, coor_t max, int m) {
	int min = 0;
	int num_points = 0;
	int i, j;
	float a = 0.5;
	int len_limit = (int) max;
	int sample_size = -1;
	float normalization_sum = 0;

	float *sigma;
	float *pi_prime = (float*) zMalloc(sizeof(float)*(len_limit+1), "zKernelSmooth: pi_prime");
	zIVec *lengths = (zIVec*) zMalloc(sizeof(zIVec), "zKernelSmooth: lengths");
	zIVec *counts  = (zIVec*) zMalloc(sizeof(zIVec), "zKernelSmooth: counts");

	zInitIVec(lengths, 100);
	zInitIVec(counts,  100);

	for (i = min + 1; i <= len_limit; i++) {
		if (count[i] >= 1) {
			zPushIVec(lengths, i);
			zPushIVec(counts, count[i]);
			num_points += (count[i]);
		}
	}

	sample_size = counts->size;

	sigma    = (float*) zMalloc(sizeof(float)*sample_size, "zKernelSmooth: sigma");

	for (i = 0; i < sample_size; i++) {
		int left_neighbor = -1, right_neighbor = -1;
		if (i - m >= 0) {
			left_neighbor = 1 + lengths->elem[i] - lengths->elem[i - m];
		}
		if (i + m < sample_size) {
			right_neighbor= 1 + lengths->elem[i + m] - lengths->elem[i];
		}
		if (left_neighbor > -1 && right_neighbor > -1) {
			sigma[i] = MIN(left_neighbor, right_neighbor);
		} else if (left_neighbor == -1 && right_neighbor == -1) {
			zWarn("Both neighbors are undefined");
			zWarn("Switching to m=%d", m-1);
			zFree(sigma);
			zFree(pi_prime);
			zFreeIVec(counts);
			zFreeIVec(lengths);
			zKernelSmooth(count, smooth, max, m-1);
			return;
		} else if (left_neighbor == -1) {
			sigma[i] = right_neighbor;
		} else if (right_neighbor == -1) {
			sigma[i] = left_neighbor;
		}
		sigma[i] = MAX(sigma[i], a*lengths->elem[i]/pow(num_points, 0.2));
	}

	for (i = 1; i <= len_limit; i++) {
		float sum = 0;
		for (j = 0; j < lengths->size; j++) {
			sum += (counts->elem[j] * zKernelFunction(sigma[j], i - lengths->elem[j]));
		}
		sum = 1.0*sum/num_points;
		pi_prime[i] = sum;
		normalization_sum += sum;
	}

	for (i = 1; i <= len_limit; i++) {
		smooth[i-1] = zFloat2Score(1.0 * pi_prime[i] / normalization_sum);
	}
	zFree(sigma);
	zFree(pi_prime);
	zFreeIVec(counts);
	zFreeIVec(lengths);
}
