/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/*****************************************************************************\
mathtest.c

math test program for zoe
\*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ZOE.h"

void baseTest (void);
void funcTest (void);
void distTest (char*);
void duraTest (char*);
void addTest  (void);

int main (int argc, char *argv[]) {

	/* set the program name */
	zSetProgramName(argv[0]);

	if (argc == 1) {
		fprintf(stderr, "usage: %s <test> [test_arg]\n", argv[0]);
		fprintf(stderr, "tests: {func, base, add, dist [file], dura [file]}\n");
		exit(1);
	}
	
	     if (strcmp(argv[1], "func") == 0) funcTest();
	else if (strcmp(argv[1], "dist") == 0) distTest(argv[2]);
	else if (strcmp(argv[1], "dura") == 0) duraTest(argv[2]);
	else if (strcmp(argv[1], "base") == 0) baseTest();
	else if (strcmp(argv[1], "add")  == 0) addTest();
	else zDie("mathtest unknown test (%s)", argv[1]);
	
	return(0);
}

void baseTest (void) {
	char b[64];
	char s[10] = "101";
	printf("%d %d %d\n", zBaseToDec(2, s), zBaseToDec(4, s), zBaseToDec(5, s));
	
	zDecToBase(5, 2, b);
	printf("%s\n", b);
	zDecToBase(17, 4, b);
	printf("%s\n", b);
	zDecToBase(26, 5, b);
	printf("%s\n", b);
	
}

void funcTest (void) {
	int i;
	float f;
	
	/* score-float conversion functions */
	for (i = 10; i <= 25; i+=5) {
		f = (float)i/100;
		printf("zFloat2Score(%f) = %f\n", f, zFloat2Score(f));
	}
	for (i = 10; i <= 25; i+=5) {
		printf("zScore2Float(%d) = %f\n", i, zScore2Float(i));
	}
	
	/* zLog2 */
	for (i = 2; i <= 8; i+= 2) {
		printf("zLog2(%d) %f\n", i, zLog2(i));
	}
	
	/* zLnFactorial */
	for (i = 2; i <= 8; i+= 2) {
		printf("zLnFactorial(%d) %f\n", i, zLnFactorial(i));
	}

	/* zScorePoisson */
	for (i = 0; i <= 8; i+=2) {
		printf("zScorePoisson(4,%d) %f\n", i, zScorePoisson(4, i));
	}
	
	/* zScoreGeometric */
	for (i = 1; i <= 5; i++) {
		printf("zScoreGeometric(2,%d) %f\n", i, zScoreGeometric(2, i));
	}
}

void distTest (char *filename) {
	FILE *stream;
	zDistribution d;
	int i;
	score_t s;
	
	if ((stream = fopen(filename, "r")) == NULL)
		zDie("File (%s) not found", filename);
	if (!zReadDistribution(stream, &d))
		zDie("mathtest: zReadDistribution (%s) failed", filename);
	zWriteDistribution(stdout, &d);
	for (i = 1; i < 10; i++) {
		s = zScoreDistribution(&d, i);
		printf("%d %f\n", i, s);
	}
	zFreeDistribution(&d);
	fclose(stream);
}

void duraTest (char *filename) {
	FILE *stream;
	zDurationGroup d;
	int i, j;
	score_t s;
	
	if ((stream = fopen(filename, "r")) == NULL)
		zDie("File (%s) not found", filename);
	if (!zReadDurationGroup(stream, &d))
		zDie("mathtest: zReadDurationGroup (%s) failed", filename);
	zWriteDurationGroup(stdout, &d);
	for (j = 0; j < d.durations; j++) {
	  printf("Isochore upper bound: %f\n", d.iso_bound[j]);
	  for (i = 1; i < 40; i++) {
			s = zScoreDurationGroup(&d, i, 0.99*d.iso_bound[j]);
			printf("%d %f\n", i, s);
		}
	}
	zFreeDurationGroup(&d);
	fclose(stream);
}

score_t realAdd(score_t p, score_t q) {
	return zFloat2Score(zScore2Float(p) + zScore2Float(q));
}
void addTest(void) {
	score_t p,q;

	p = 0; q =5;
	printf("p:%g q:%g  ~p+~q:%g ~q+~p:%g q+p:%g\n", p, q, 
		   zFloatwiseScoreAdd(p,q), zFloatwiseScoreAdd(q,p), realAdd(p,q));
	p = 5; q =5;
	printf("p:%g q:%g  ~p+~q:%g ~q+~p:%g q+p:%g\n", p, q, 
		   zFloatwiseScoreAdd(p,q), zFloatwiseScoreAdd(q,p), realAdd(p,q));
	p = 0; q =0;
	printf("p:%g q:%g  ~p+~q:%g ~q+~p:%g q+p:%g\n", p, q, 
		   zFloatwiseScoreAdd(p,q), zFloatwiseScoreAdd(q,p), realAdd(p,q));
	p = -50; q =50;
	printf("p:%g q:%g  ~p+~q:%g ~q+~p:%g q+p:%g\n", p, q, 
		   zFloatwiseScoreAdd(p,q), zFloatwiseScoreAdd(q,p), realAdd(p,q));
	p = -2035; q = -2036;
	printf("p:%g q:%g  ~p+~q:%g ~q+~p:%g q+p:%g\n", p, q, 
		   zFloatwiseScoreAdd(p,q), zFloatwiseScoreAdd(q,p), realAdd(p,q));
	p = -1; q =-3;
	printf("p:%g q:%g  ~p+~q:%g ~q+~p:%g q+p:%g\n", p, q, 
		   zFloatwiseScoreAdd(p,q), zFloatwiseScoreAdd(q,p), realAdd(p,q));
	p = -10.203211; q =-25.8282;
	printf("p:%g q:%g  ~p+~q:%g ~q+~p:%g q+p:%g\n", p, q, 
		   zFloatwiseScoreAdd(p,q), zFloatwiseScoreAdd(q,p), realAdd(p,q));
	p = -1; q =-3;
	printf("p:%g q:%g  ~p+~q:%g ~q+~p:%g q+p:%g\n", p, q, 
		   zFloatwiseScoreAdd(p,q), zFloatwiseScoreAdd(q,p), realAdd(p,q));
	
}

