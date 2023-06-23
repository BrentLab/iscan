#include <stdio.h>
#include <stdlib.h>
#include "ZOE.h"


size_t ASCII_SET = 256;

int main (int argc, char *argv[]) {
	FILE* file;
	zFastaFile seq;
	double entropy = 0;
	double d;
	float f;
	int file_count = 0;
	int total_length = 0;
	int this_arg;
	coor_t i;
	int count[256];
	for(i=0;i<ASCII_SET;i++) count[i] = 0;
	
	/* set the program name */
	zSetProgramName(argv[0]);
	
	if (argc == 1) {
		(void)fprintf(stderr, "usage: %s <fasta files>\n", argv[0]);
		exit(1);
	}
	
	for (this_arg=1;this_arg<argc;this_arg++) {
	
		/* open the file and get the sequences */
		file = fopen(argv[this_arg], "r");
		while (zReadFastaFile(file, &seq)) {
			file_count++;
			
			/* count the letters */
			for (i=0;i<seq.length;i++) count[(int) seq.seq[i]]++;
			total_length += seq.length;
			
			zFreeFastaFile(&seq);
		}
		(void)fclose(file);
	}
	
	/* output some stats */
	(void)printf("%d files\n", file_count);
	(void)printf("%d total letters\n", total_length);
	f = total_length/file_count;
	(void)printf("%g average\n", f);
	
	for (i=0;i<ASCII_SET;i++) {
		if (count[i] != 0) {
			f = (double)count[i]/total_length;
			d = - ((double)count[i]/total_length)
				* zLog2((double)count[i]/total_length);
			(void)printf("%c\t%d\t%f\t%g\n", i, count[i], f, d);
			entropy += d;
		}
	}
	(void)printf("entropy = %g bits\n", entropy);
	return 0;
}
