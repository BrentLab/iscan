###############################################################################
# Makefile for zoe
###############################################################################

###########
# LIBRARY #
###########

LIBRARY = libzoe.a

OBJECTS = \
	zTools.o\
	zTBTree.o\
	bntree/bntree.o\
	bntree/cache.o\
	bntree/celllistitem.o\
	bntree/intlistitem.o\
	bntree/node.o\
	zAlignment.o\
	zAlnFeature.o\
	zConseq.o\
	zDistribution.o\
	zDNA.o\
	zDuration.o\
	zEstseq.o\
	zEstseqSearch.o\
	zFastaFile.o\
	zFeatureFactory.o\
	zGTF.o\
	zHardCoding.o\
	zHMM.o\
	zHMM_State.o\
	zMath.o\
	zMathTables.o\
	zModel.o\
	zPhasePref.o\
	zProtein.o\
	zScanner.o\
	zSeedUtils.o\
	zSequence.o\
	zSfeature.o\
	zStopSeq.o\
	zTransition.o\
	zTrellis.o\
	zViterbi.o\
	zPairTransition.o\
	zPairTrellis.o\
	zPairViterbi.o\

#################
# TEST PROGRAMS #
#################

LFLAGS = -L. -lzoe -lm `pkg-config --libs glib-2.0`

EXE0 = tooltest
SRC0 = tooltest.c
OBJ0 = $(SRC0:.c=.o)

EXE1 = seqstats
SRC1 = seqstats.c
OBJ1 = $(SRC1:.c=.o)

EXE2 = seqedit
SRC2 = seqedit.c
OBJ2 = $(SRC2:.c=.o)

EXE3 = mathtest
SRC3 = mathtest.c
OBJ3 = $(SRC3:.c=.o)

EXE4 = modeltest
SRC4 = modeltest.c
OBJ4 = $(SRC4:.c=.o)

EXE5 = scannertest
SRC5 = scannertest.c
OBJ5 = $(SRC5:.c=.o)

EXE6 = sftest
SRC6 = sftest.c
OBJ6 = $(SRC6:.c=.o)

EXE7 = factorytest
SRC7 = factorytest.c
OBJ7 = $(SRC7:.c=.o)

EXE8 = gtftest
SRC8 = gtftest.c
OBJ8 = $(SRC8:.c=.o)

EXE9 = conseqtest
SRC9 = conseqtest.c
OBJ9 = $(SRC9:.c=.o)

EXE11 = conseqtest2
SRC11 = conseqtest2.c
OBJ11 = $(SRC11:.c=.o)

EXE12 = iscan
SRC12 = iscan.c
OBJ12 = $(SRC12:.c=.o)

EXE13 = zoe2gtf
SRC13 = zoe2gtf.c
OBJ13 = $(SRC13:.c=.o)

EXE14 = gtf2zoe
SRC14 = gtf2zoe.c
OBJ14 = $(SRC14:.c=.o)

EXE15 = pairagon
SRC15 = pairagon.c
OBJ15 = $(SRC15:.c=.o)

EXE16 = pairagon2estgen
SRC16 = pairagon2estgen.c
OBJ16 = $(SRC16:.c=.o)

EXE17 = pairameter_estimate
SRC17 = pairameter_estimate.c
OBJ17 = $(SRC17:.c=.o)

EXECUTABLES = $(EXE12) $(EXE13)
PAIRAGON_EXECUTABLES = $(EXE15) $(EXE16) $(EXE17)
#EXECUTABLES = $(EXE12) $(EXE0) $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE7) $(EXE8) $(EXE9) $(EXE10) $(EXE11) $(EXE12) $(EXE13) $(EXE14)

TABLE = zMathTables.c
TABLEGEN = zFloatwiseScoreAdd

########################
# Program Options      #
########################

BUILD = 20060622MA
VERSION = 3.5
PROGRAM = Twinscan
PAIRAGON_VERSION = 0.99
PAIRAGON_PROGRAM = Pairagon
GCC_FLAGS = "-W -Wall -Werror -ansi -pedantic -Wwrite-strings `pkg-config --cflags glib-2.0`"
DIST_FLAGS = "-DBUILD=\\\"$(PROGRAM)\ version\ $(VERSION)\ build\ $(BUILD)\\\" -DNDEBUG -O2"
PAIRAGON_DIST_FLAGS = "-DBUILD=\\\"$(PAIRAGON_PROGRAM)\ version\ $(PAIRAGON_VERSION)\ build\ $(BUILD)\\\" -DNDEBUG -O2"

############################
# End Program Options      #
############################

################################################################################
# TARGETS 
################################################################################

include Makefile.target

default:
	echo "recalled_target: default" > Makefile.target
	make gcc-development

debug:
	echo "recalled_target: debug" > Makefile.target
	make gcc-debug

product:
	echo "recalled_target: product" > Makefile.target
	make gcc-distribution

product64:
	echo "recalled_target: product64" > Makefile.target
	make gcc-distribution-64

distro:
	echo "recalled_target: distro" > Makefile.target
	make gcc-distro

profile:
	echo "recalled_target: profile" > Makefile.target
	make gcc-profile

insure:
	echo "recalled_target: insure" > Makefile.target
	make all CC="insure" CFLAGS+=-DDEBUG CFLAGS+=$(GCC_FLAGS) CFLAGS+=-g 

insure_link:
	echo "recalled_target: insure_link" > Makefile.target
	make CFLAGS+=$(GCC_FLAGS) CFLAGS+=-g gcc-debug
	rm -f cp_iscan iscan snp_iscan pin_iscan pin_iscan2
	make all CC="insure" CFLAGS+=-DDEBUG CFLAGS+=$(GCC_FLAGS) CFLAGS+=-g 

$(LIBRARY): $(OBJECTS)
	rm -f $(LIBRARY)
	ar cr $(LIBRARY) $(OBJECTS)
	ranlib $(LIBRARY)

all: $(LIBRARY) $(EXECUTABLES)

pairagon-exe: $(LIBRARY) $(PAIRAGON_EXECUTABLES)

clean:
	rm -f *.o bntree/*.o $(LIBRARY) $(EXECUTABLES) core *~ *.bck $(TABLE) $(TABLEGEN)

include:
	gcc `pkg-config --cflags glib-2.0` -MM *.c > Makefile.include

test:
	perl test.pl

zoe:
	make clean
	make default
	make test

Makefile.target:
	touch Makefile.target

gcc-development:
	make all CC="gcc" CFLAGS+=$(GCC_FLAGS) CFLAGS+=-DTRACE=1 

gcc-debug:
	make all CC="gcc" CFLAGS+="-DTRACE=2" CFLAGS+=-DDEBUG CFLAGS+=$(GCC_FLAGS) CFLAGS+=-g CFLAGS+=-pg CFLAGS+=-O0

gcc-distribution:
	make all CC="gcc" CFLAGS+=$(GCC_FLAGS) CFLAGS+=-DNDEBUG CFLAGS+=-O2 CFLAGS+=$(DIST_FLAGS)

gcc-distribution-64:
	make all CC="gcc" CFLAGS+=$(GCC_FLAGS) CFLAGS+=-march=opteron CFLAGS+=-DNDEBUG CFLAGS+=-O2 CFLAGS+=$(DIST_FLAGS)

gcc-distro:
	make all CC="gcc" CFLAGS+=$(GCC_FLAGS) CFLAGS+=-static CFLAGS+=-DNDEBUG CFLAGS+=-O2 CFLAGS+=$(DIST_FLAGS)

gcc-profile:
	make all CC="gcc" CFLAGS+=$(GCC_FLAGS) CFLAGS+=-pg

pairagon-debug:
	make pairagon-exe CC="gcc" CFLAGS+="-DTRACE=2" CFLAGS+=-DDEBUG CFLAGS+=$(GCC_FLAGS) CFLAGS+=-g CFLAGS+=-pg

pairagon-linux:
	make pairagon-exe CC="gcc" CFLAGS+=$(GCC_FLAGS) CFLAGS+=$(PAIRAGON_DIST_FLAGS)

pairagon-sparc:
	make pairagon-exe CC="cc" CFLAGS+=$(SPARC_FLAGS) CFLAGS+=$(PAIRAGON_DIST_FLAGS)

pairagon-distro:
	make pairagon-exe CC="gcc" CFLAGS+=$(GCC_FLAGS) CFLAGS+=-march=i686 CFLAGS+=-static CFLAGS+=$(PAIRAGON_DIST_FLAGS)

pairagon-insure:
	make pairagon-exe CC="insure" CFLAGS+="-DTRACE=2" CFLAGS+=$(GCC_FLAGS) CFLAGS+=-g

########################
# Architecture Section #
########################

ALPHA_FLAGS = "-c99 -msg_enable level5"

alpha:
	make all CC="cc" CFLAGS+=$(ALPHA_FLAGS) CFLAGS+=$(DIST_FLAGS)

MACOSX_FLAGS+="-W -Wall -Werror -ansi -pedantic -g -Wwrite-strings -DMAC_OS_X"

macosx:
	make all CC="gcc" CFLAGS+=$(MACOSX_FLAGS) CFLAGS+=$(DIST_FLAGS)

SPARC_FLAGS = ""

sparc:
	make all CC="cc" CFLAGS+=$(SPARC_FLAGS) CFLAGS+=$(DIST_FLAGS)

#######################
# Excecutable Section #
#######################

$(EXE0): $(OBJ0) $(LIBRARY)
	$(CC) -o $(EXE0) $(CFLAGS) $(OBJ0) $(LFLAGS)

$(EXE1): $(OBJ1) $(LIBRARY)
	$(CC) -o $(EXE1) $(CFLAGS) $(OBJ1) $(LFLAGS)

$(EXE2): $(OBJ2) $(LIBRARY)
	$(CC) -o $(EXE2) $(CFLAGS) $(OBJ2) $(LFLAGS)

$(EXE3): $(OBJ3) $(LIBRARY)
	$(CC) -o $(EXE3) $(CFLAGS) $(OBJ3) $(LFLAGS)

$(EXE4): $(OBJ4) $(LIBRARY)
	$(CC) -o $(EXE4) $(CFLAGS) $(OBJ4) $(LFLAGS)

$(EXE5): $(OBJ5) $(LIBRARY)
	$(CC) -o $(EXE5) $(CFLAGS) $(OBJ5) $(LFLAGS)

$(EXE6): $(OBJ6) $(LIBRARY)
	$(CC) -o $(EXE6) $(CFLAGS) $(OBJ6) $(LFLAGS)

$(EXE7): $(OBJ7) $(LIBRARY)
	$(CC) -o $(EXE7) $(CFLAGS) $(OBJ7) $(LFLAGS)

$(EXE8): $(OBJ8) $(LIBRARY)
	$(CC) -o $(EXE8) $(CFLAGS) $(OBJ8) $(LFLAGS)

$(EXE9): $(OBJ9) $(LIBRARY)
	$(CC) -o $(EXE9) $(CFLAGS) $(OBJ9) $(LFLAGS)

$(EXE10): $(OBJ10) $(LIBRARY)
	$(CC) -o $(EXE10) $(CFLAGS) $(OBJ10) $(LFLAGS)

$(EXE11): $(OBJ11) $(LIBRARY)
	$(CC) -o $(EXE11) $(CFLAGS) $(OBJ11) $(LFLAGS)

$(EXE12): $(OBJ12) $(LIBRARY)
	$(CC) -o $(EXE12) $(CFLAGS) $(OBJ12) $(LFLAGS)

$(EXE13): $(OBJ13) $(LIBRARY)
	$(CC) -o $(EXE13) $(CFLAGS) $(OBJ13) $(LFLAGS)

$(EXE14): $(OBJ14) $(LIBRARY)
	$(CC) -o $(EXE14) $(CFLAGS) $(OBJ14) $(LFLAGS)

$(EXE15): $(OBJ15) $(LIBRARY)
	$(CC) -o $(EXE15) $(CFLAGS) $(OBJ15) $(LFLAGS)

$(EXE16): $(OBJ16) $(LIBRARY)
	$(CC) -o $(EXE16) $(CFLAGS) $(OBJ16) $(LFLAGS)

$(EXE17): $(OBJ17) $(LIBRARY)
	$(CC) -o $(EXE17) $(CFLAGS) $(OBJ17) $(LFLAGS)

$(EXE18): $(OBJ18) $(LIBRARY)
	$(CC) -o $(EXE18) $(CFLAGS) $(OBJ18) $(LFLAGS)

###################
# Inference Rules #
###################

.SUFFIXES : .c .o

.c.o:
	$(CC) $(CFLAGS) -c -o $@ $<


#######################
# Source Dependencies #
#######################

$(TABLEGEN): zFloatwiseScoreAdd.c
	$(CC) $(CFLAGS) -o $(TABLEGEN) zFloatwiseScoreAdd.c -lm

$(TABLE): $(TABLEGEN)
	./$(TABLEGEN) > $(TABLE)

include Makefile.include
