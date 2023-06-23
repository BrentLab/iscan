/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\

  zViterbi

  This holds the data needed for the memory optimized viterbi run.
  The same structure is used by both GHMM and GPAIRHMM. zViterbi.c
  implements all the functions that are used in both, and some that are
  used in GHMM exclusively, and zPairViterbi.c implements all the 
  functions that are used in GPAIRHMM exclusively.

  See also:
	zViterbi.c, zPairViterbi.c
  
  Copyright (C) 2006 Evan Keibler and Manimozhiyan Arumugam
\******************************************************************************/

#include "zSfeature.h"
#include "zAlnFeature.h"
#include "zHMM.h"
#include "zTrellis.h"
#include "zPairTrellis.h"
#include "zTBTree.h"

/* this struct holds all the traceback pointers/scores */
typedef struct zViterbiTraceback{
 	zTBTree        tbtree;           /* main traceback tree */
 	zPtrList       live_nodes;       /* future live nodes for states which have no fixed
										upper bound on their length */   
	/* the short_live_nodes stores future nodes with short/medium, dense tb pointers 
       in a memory efficient implementation than a full array but more computationally 
       efficient than just using the TB tree */
 	zTBTreeNode*** short_live_nodes; /* future live nodes for states which have fixed
										upper bounds */
	/* the cache stores past nodes and was added specifically for GPAIRHMM.  It is used to 
       store leaf nodes in a GHMM */
	zTBTreeNode**** cache;
	/* cache[cache_level][cdna_pos][state] is a zTBTreeNode*.
	   cdna_pos=0 for regular GHMM */
} zViterbiTraceback;

/* zViterbi structure holds the trellis, traceback tree and tbnodes */ 
struct zViterbi{

	/* The trellis pointers: only one of these can be NON-NULL at any time */

	zTrellis*      trellis;
	zPairTrellis*  pair_trellis;

	zHMM*          hmm;

 	zViterbiTraceback** tb;
 	int*           sln_map;
 	int*           sln_size;
 	int            sln_count;

	zTBTreeNode**  next_tbtn;
	int*           next_tbti;
	short*         tbt_use_count;
	
	int            cache_length;
	/* How many previous positions should be cached? *
	 * For geometric only, you just need 2, if you   *
	 * have fixed length states then you need that   *
	 * length + 1, if you have explicit states, then *
	 * you need the max_explicit_length + 1          */
	
	/* For SNP stuff */
	short*         snp_count;
	short*         active_snp_count;
	int*           first_snp;
	int**          snp_idx;
	int*           temp;
	short**        snp_vals;
	int*           snp_int_vals;

	zSFList*       sfl;
	zAFList*       afl;
	int            max_trellis_count;
	short          max_snp_count;
	short          max_active_snp_count;
	float          hap_threshold;
	int            trellis_count;
	zPtrList*      traceback;
	coor_t         pos;
	coor_t         cpos;

	/* These are added for compatibility with PairHMMs */

	coor_t         cdna_pos;
	coor_t         cdna_min;
	coor_t         cdna_max;
	
	int            active_count;
	int            dead_count;
	short*         active;
	coor_t*        start_pos;
	int*           next;

	float***       hapmap;
	short**        hapmax;

	int heuristic;
 	coor_t heuristic_spacing;
 	int trimmed_now;
 	int trimmed_later;
};
typedef struct zViterbi zViterbi;

void         zInitViterbiTraceback(zViterbiTraceback* tb, zViterbi* v);
void         zFreeViterbiTraceback(zViterbiTraceback* tb, zViterbi* v);
void         zResetViterbiTraceback(zViterbiTraceback* tb, zViterbi* v);
zTBTreeNode* zFindViterbiLiveNode(zViterbi *v,int allele, coor_t pos,int state);
void         zRemoveViterbiLiveNode(zViterbi *v,int t, zTBTreeNode* rem);
int          zViterbiIncorporateLiveNodes(zViterbi* v, int allele);
int          zBackwardViterbiIncorporateLiveNodes(zViterbi* v, int allele);
void         zSNPCopyAllele(zViterbi* v, int a1, int a2);
void         zSNPSwapAlleles(zViterbi* v, int a1,int a2);
void         zSNPReleaseAllele(zViterbi* v, int a);
void         zPinReleaseTrellis(zViterbi* v, int a);
void         zSnpViterbiMergeSFLists(zSFList* base, zSFList* ext);
void         zSNPCollapseAllele(zViterbi* v, int allele);
void         zSNPCleanUpDecodes(zViterbi *v);

struct zAlleleDecode{
	char*   header;
	zSFVec* sfv;
	zAFVec* afv;
	zSFList* sfl;
	zAFList* afl;
	coor_t  start;
	coor_t  stop;
	coor_t  diff_start;
	coor_t  diff_stop;
	short*  vals;
	int     snp_count;
	int     first_snp;
	score_t score_diff;
};
typedef struct zAlleleDecode zAlleleDecode;

void         zFreeAlleleDecode(zAlleleDecode*);
void         zSNPFillAlleleHeader(zViterbi* v,zAlleleDecode* ad);
