/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
  Traceback Tree

  This is a tree structure which can store all live traceback paths 
  from a trellis more efficiently than retaining the full trellis.
  
  Copyright (C) 2006 Evan Keibler and Manimozhiyan Arumugam
\******************************************************************************/

#ifndef ZTBTREE
#define ZTBTREE

#include "zTools.h"

typedef struct zTBTree {
	struct zTBTreeNode*  root;   /* this is the tree root */
	struct zTBTreeNode*  cpoint; /* the cpoint of this tree */
	struct zTBTreeNode*  dead_nodes;  /* this faux tree holds the dead nodes stored for reuse
								  the nodes are stored in a single linked list using the 
								  child pointer */
	int           size;   /* total number of nodes in tree */
	short         new_cpoint; /*flag set to 1 when cpoint changes */
	short         id;
	short         check_count;
} zTBTree;

typedef struct zTBTreeNode {
	coor_t           pos;
	coor_t           cdna_pos;
	int              state;
	score_t          score;
	strand_t         strand;
	frame_t          frame_data;
	zPhase_t         phase;
	struct zTBTreeNode*     parent;
	struct zTBTreeNode*     child;
	struct zTBTreeNode*     lsib;
	struct zTBTreeNode*     rsib;
	int              children;
	int              id;
	int              lock;
	struct zTBTree*         tree;
	
} zTBTreeNode;

/* each node sees its leftmost child only.  other accessed through the sibling
   pointers.  each child knows its parent. */

void zInitTBTree(zTBTree *t);
void zFreeTBTree(zTBTree *t);
void zResetTBTree(zTBTree* tree);
bool zTBTreeCheckNewCpoint(zTBTree* tree);
bool zTBTreeClearNewCpoint(zTBTree* tree);
void zReleaseDeadTBTreeNode(zTBTree* t, zTBTreeNode* n);
void zTBTreeLockNode(zTBTreeNode* n);
void zTBTreeClearSubTree(zTBTree* t, zTBTreeNode* n);
void zCopyTBTreeNode(zTBTree* t1, zTBTreeNode* n1, zTBTree* t2, zTBTreeNode* n2, int depth);
void zAppendTBTreeFromNode(zTBTree* t1, zTBTreeNode* n1, zTBTree* t2, zTBTreeNode* n2);
void zCopyTBTreeFromCPoint(zTBTree* t1,zTBTree* t2,coor_t min_pos);
zTBTreeNode* zGetTBTreeNode(zTBTree *t);
void zTBTreeSetChild(zTBTreeNode* p,zTBTreeNode* c);
void zClearTBTreeNodeChildren(zTBTree* t, zTBTreeNode* n);
void zReleaseRedundantTBTreeNode(zTBTree* t,zTBTreeNode* n);
void zReleaseTBTreeNode(zTBTree* t,zTBTreeNode* n);

#endif
