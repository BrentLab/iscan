TWINSCAN/N-SCAN 3.0 DISTRIBUTION
--------------------------------

CONTENTS:

This README file is split into the following parts:

  I. TWINSCAN/N-SCAN Documentation
 II. Twinscan specific Documentation
III. N-SCAN specific Documentation

Address questions or comments to nscan@mblab.wustl.edu.


I. TWINSCAN/N-SCAN Documentation
--------------------------------

No panic.
This distribution comes with example inputs and outputs, and two pipeline scripts:
one for Twinscan and one for N-SCAN. These scripts generate most of the inputs to
the actual Twinscan/N-SCAN executable, and runs it, too. 
For an example of their use, see the Quickstart Guides under II and III.

1. Twinscan/N-SCAN Executable
==============================

Twinscan/N-SCAN 4.0 is written in C and the source code, if included, is found in 
the src/ directory.  

If you downloaded the N-SCAN_(build nr)_src.tar.gz, you must create an executable.
To create the executable for linux, type 

make linux

from the main directory.  The Twinscan/N-SCAN executable will be placed in the bin/ 
directory. 
If you wish to build for a different architecture/OS than Linux Intel, several
special Makefile targets have been made and tested for this: alpha, macosx, 
macosxg5, and sparc.  

There are several pre-built architectures for N-SCAN available through
http://mblab.wustl.edu/software/Twinscan_Nscan/

To build Pairagon, use one of the following make targets: 
pairagon-linux, pairagon-sparc, pairagon-macosx or pairagon-macosxg5.

You may wish to add the Twinscan/N-SCAN environment variable to your login script.
For example, under Linux with a bash shell, type

export TWINSCAN="/home/bart/Twinscan"

If this is a binary only (pre-built) distribution, the Twinscan/N-SCAN executable will already
be in the bin/ directory.  To test the executable, type  

./test-executable   

in the main directory.

Twinscan/N-SCAN has been tested under Linux with gcc version 3.4.6
HOWEVER THERE IS NO GUARANTEE THAT TWINSCAN/N-SCAN WILL COMPILE OR RUN WITH YOUR 
PARTICULAR COMBINATION.  

2. Twinscan/N-SCAN Parameter Files
==================================

Each Twinscan/N-SCAN HMM parameter file is specific for a particular target-informant pair
such as human-mouse. The reason for this is that the Twinscan/N-SCAN parameters are
sensitive to evolutionary distance.  Development of new parameter files is an 
ongoing project. If your favorite pairs of genomes are not present, you may be 
able to use something similar and still get satisfactory results. Alternatively, 
you may download iParameterEstimation from http://mblab.wustl.edu/software/ and create
your own parameters.

3.  Version History
===================
4.0        Package updated with parameterfiles, and reorganized. Nscan_driver script was added.
3.5
3.0  build 20051110RB, Nov. 2005
           Merged Twinscan/N-SCAN in a single distribution. Documentation updated to reflect the merge.
2.03 build 20051004CW, Oct. 2005
           Code optimization and several bugfixes.
           Added HMM parameter files for mammalian (human/mouse/rat) gene prediction with EST evidence.
        Dec 2004
           Added HMM and blast parameter files for mammalian (human/mouse/rat) gene prediction.
           The mammalian HMM parameter file does not include the EST mode.
2.02 build 20041011CW, Oct. 2004
           Added gene prediction with EST support.
           Added the ability to use WWAMs for INITCONS and TERMCONS models.
           Fixed a bug that counted signal scores (start, stop codon and 
           splice sites) twice for the conservation sequence model.
           Added parameter file for Arabidopsis gene prediction.
2.01 build 20040819MA, Aug. 2004
           Patch release for version 2.01. August 2004.
2.01, July 2004
           Explicit intron length model was implemented. New parameter
           file for C. elegans with explicit intron length was added.
2.0 beta, Feb. 2004
           Parameter file for Cryptococcus added February 2004.
2.0 beta, Dec. 2003  
           Released version 2.0 beta.



II. TWINSCAN SPECIFIC DOCUMENTATION
-----------------------------------

This part of the file contains the following sections:

A.  Quick Start Guide
B.  Twinscan Overview
C.  Running Twinscan - Basic Instructions
D.  Running Twinscan using the runTwinscan2 script
E.  Known Limitations


A.  Quick Start Guide
=====================

An example script (described in detail below) for the Twinscan analysis pipeline is included.
To access, go to the /examples directory and run

../bin/runTwinscan2 -r ../parameters/twinscan_parameters/human_twinscan.zhmm -d output -B ../parameters/blast_params/Hsapiens.blast.param example.fa.masked informant.fa

After running you can find output files in the newly created /output directory.

Several programs must be installed on your system to run runTwinscan
You may need to change runTwinscan to point it to these programs. To do so, open the
script in a text editor and look for the following:
my $REPEATMASKER        = "RepeatMasker";       # Format for local environment
my $BLASTN              = "blastn";             # Format for local environment
my $BLAT                = "blat";               # Format for local environment
my $XDFORMAT            = "xdformat";           # Format for local environment
my $PRESSDB             = "pressdb";            # Format for local environment

Example: If RepeatMasker is in "/bin/john/RepeatMasker", put this in place of "RepeatMasker".
Alternatively, you can add these programs to your path.

B.  Twinscan Overiew
=====================

Twinscan finds genes in a "target" genomic sequence by simultaneously
maximizing the probability of the gene structure in the target and the
evolutionary conservation dervied from "informant" genomic sequences.

The target sequence (i.e. the sequence to be annotated) should generally be
of draft or finished quality.  The informant can range from a single sequence 
to a whole genome in any condition from raw shotgun reads to finished assembly.  
Details about how the quality of the informant database effects predictive 
accuracy can be found in Flicek, et. al. 

Information complementary to this file can be found in the following:

P. Hu and M.R. Brent.  Using Twinscan to predict gene structures in
genomic DNA sequence.  Current Protocols in Bioinformatics (in press).

If you use Twinscan in your research, please cite the following
references:

P. Flicek, E. Keibler, P. Hu, I. Korf, M.R. Brent. Leverging the mouse
genome for gene prediction in human: from whole genome shotgun reads
to a global synteny map.  Genome Research 13. 46-54.

I. Korf, P. Flicek, D. Duan, M.R. Brent. 2001.  Integrating genomic
homology into gene-structure prediction.  Bioinformatics 17. S140-S148.

In order to run Twinscan you will need the following components:

  (1) Nscan 4.0 executable
  (2) Twinscan parameter file
  (3) DNA sequence
  (4) Conservation sequence
  (5) EST sequence (optional)


(1) Twinscan Executable
-----------------------

See Section I.1

(2) Twinscan Parameter File
---------------------------

The parameterfiles can be found in /parameters/twinscan_parameters. Each filename contains the
name of the target organism that was used to create it, eg maize_twinscan.zhmm. Twinscan results
will be optimal for this species, but it may be possible to use it for a related organism.

(3) DNA Sequence
----------------

The target sequence must be in FASTA format, must be longer that 500 bp, and should
have the repetitive elements masked.  While masking is not required to run Twinscan, 
it will improve performance by reducing false-positive predictions.  

We normally mask with RepeatMasker, go to (http://www.repeatmasker.org/RMDownload.html).
The RepeatMasker program is not included with this distribution.
Note: RepeatMasker will also mask low complexity and simple repeats. We recommend 
switching this off by using the -nolow flag. Real genes sometimes contain such repeats, 
and we find that gene prediction works better this way.

(4) Conservation Sequence
-------------------------

Conservation sequence is a symbolic representation of the the best alignments
between the target and informant sequences. The format of the conservation
sequence file is very simple: a definition line that includes the BLAST
database name and a second line of conservation symbols (which are
just numbers). For an example, see examples/example.conseq.

To create this conservation sequence, you need a BLAST program. 
We generally use WU-BLAST (http://blast.wustl.edu) to create the BLAST 
report.  NCBI BLAST works with our software, but the input parameters
need to be changed. Parameters for WU-BLAST can be found in 
examples/example_blast_parameters.txt. 

The choice of BLAST parameters is an important consideration and will
affect both the time required for the Twinscan analysis pipeline and the 
performance of the gene-prediction algorithm.  See Flicek et. al. for
the BLAST parameters we chose to annotate the human genome. 

WU-BLAST comes with the xdformat program, which formats the informant sequences
to create a blast database. 
After running BLAST, the output must be formatted with conseq, which is 
included in this package.

Example:
xdformat -n informant.fa
Blast M=1 N=-1 Q=5 R=1 B=10000 V=100 -cpus=1 -warnings -lcfilter filter=seg filter=dust topcomboN=1 informant.fa target.fa > blast.out
conseq target.fa blast.out > conseq.fa

Note: The runTwinscan2 script will run these steps without user intervention (see below).

(5) EST Sequence
-------------------------

EST sequence is a symbolic representation of evidence from ESTs that align to
the target sequence. The format is similar to the Conservation sequence, but 
the possible values for each position are 1, 2, 0 (to represent Exon, Intron and
not known). The estseq script included in the distribution creates
EST sequence when given a DNA sequence and a (set of) BLAT reports of the
the ESTs aligned to the target.
For downloading BLAT, go to http://genome.ucsc.edu/FAQ/FAQblat.html#blat3
and follow the instructions.

C.  Running Twinscan - Basic instructions
===========================================

Twinscan takes a number of command-line parameters.  One parameter
file (e.g. human_twinscan.zhmm) and two sequence 
files (the target sequence and the conservation sequence) are required.

Twinscan's output is in GTF2 format ( see http://mblab.wustl.edu/GTF2.html).

When all files described above are present, twinscan can be run like so:

twinscan <parameter file> <masked sequence file> -c=<conseq file> [-e=estseq_file] > <outputfile>

example:
twinscan human_twinscan.zhmm mySequence.masked.fa -c=conseq.fa > mySequence.gtf

Notes:
Twinscan may be run in "Genscan-compatible" mode by skipping the "-c=<conseq>" 
option.  In this case only the zoe HMM parameters
and the target sequence are required.  

In practice, Twinscan's memory requirements are approximately linear 
with the length of the target sequence.  A rough guideline is 1 GB of memory
for 1 Mb of input sequence.  


D.  Running Twinscan using the runTwinscan2 script
======================================================

In summary, there are 5 steps required to run Twinscan:
Step 1: Mask target sequence with RepeatMasker
Step 2: Create informant BLAST database
Step 3: Run BLAST
Step 4: Create conservation sequence
(Step 4b: Create EST sequence)
Step 5: Run Twinscan 

These five steps are all contained in the example script runTwinscan2, which 
comes with this distribution. You may have to tweak this script. 
for your particular environment (See Quickstart guide).

The default BLAST parameters used by runTwinscan2 are those for C.elegans 
(see parameters/blast_params/Celegans.blast.param).  This can and should be changed 
for any other species  with the -B option to the runTwinscan2 script. We have 
also included  files specific to Cryptococcus, Arabidopsis rice, maize and 
human annotation in the parameters directory. 

The file example.output in the /examples directory contains the output 
from runTwincan2 using the BLAST parameters found in the script.  


E.  Known Limitations
======================

Genscan-compatible mode does not produce predictions that are identical
to Genscan predictions.  Specifically promoters are often predicted in 
different places and exons may be slightly different near very long introns.


III. N-SCAN SPECIFIC DOCUMENTATION
----------------------------------

This part of the file contains the following sections:   

A.  Quick Start Guide
B.  N-SCAN Overview
C.  Running N-SCAN - Basic Instructions       
D.  Running N-SCAN using Nscan_driver

A.  Quick Start Guide
=====================

An example script (described in detail below) for the N-SCAN analysis pipeline is included.
To access, go to the /examples directory and run

../bin/Nscan_driver -d nscanOutput example.fa nscandriver.config

After running you can find output files in the newly created /nscanOutput directory.

B.  N-SCAN Overiew
=====================

N-SCAN performs gene prediction on a "target" genome using information from DNA 
sequence modeling and from single or multiple genome alignments to the target. 

The target sequence (i.e. the sequence to be annotated) should generally be
of draft or finished quality.  The informant can range from a single sequence 
to a whole genome in any condition from raw shotgun reads to finished assembly.  

Information complementary to this file can be found in the following:

Gross SS, Brent MR. Using multiple alignments to improve gene prediction. J Comput Biol. 2006 Mar;13(2):379-93.

In order to run N-SCAN you will need the following components:

  (1) Nscan 4.0 executable
  (2) N-SCAN parameter file
  (3) DNA sequence
  (4) Alignment sequence     
  (5) EST sequence (optional)


(1) N_SCAN Executable
-----------------------

See Section I.1

(2) N_SCAN parameter file
-----------------------

The parameterfiles can be found in /parameters/nscan_parameters. Each filename contains the
name of the target organism that was used to create it, eg mouse_nscan.zhmm. N-SCAN results
will be optimal for this species, but it may be possible to use it for a related organism.

(3) DNA sequence
-----------------------

The target sequence must be in FASTA format, must be longer that 500 bp, and should
have the repetitive elements masked.
See the DNA sequence section under Twinscan overview (II.B.3) for masking information.

(4) Alignment sequence
-----------------------

The informant-alignment fragment consists 
of a FASTA header line and one line for each informant (note that the DNA sequence is 
present in the informant-alignment file, but not the informant-alignment fragment). The 
length of each informant line in the informant-alignment fragment is equal to the length 
of the DNA sequence fragment to which it corresponds. For an example, see
examples/example.fa.masked.align

To create this alignment sequence, two programs are needed: blastz and lav2maf. Both can 
be downloaded from http://www.bx.psu.edu/miller_lab/. Note that lav2maf is found in the
Multiz distribution. A third program, maf_to_align, is included in this package.

Once the programs are installed, the following commands must be run:
blastz ########
lav2maf ##########
maf_to_align ########

(5) EST sequence (optional)
-----------------------

See the EST sequence section under Twinscan overview (II.B.5) for more information.


C.  Running N-SCAN - Basic Instructions       
========================================

N-SCAN takes a number of command-line parameters.  One parameter
file (e.g. human_nscan.zhmm) and two sequence files (the target 
sequence and the alignment sequence) are required.

N-SCAN's output is in GTF2 format ( see http://mblab.wustl.edu/GTF2.html).

When all files described above are present, N-SCAN can be run like so:

nscan <parameter file> <masked sequence file> -a=<align file> [-e=estseq_file] > <outputfile>

example:
nscan human_nscan.zhmm mySequence.masked.fa -a=align.fa > mySequence.gtf

Notes:
When large sequences are used, nscan can be run with memory optimization. To do this, 
add '-o' to the input:
nscan -o human_nscan.zhmm mySequence.masked.fa -a=align.fa > mySequence.gtf


D.  Running N-SCAN using Nscan_driver
========================================

In summary, there are 6 steps required to run N-SCAN:
Step 1: Mask target sequence with RepeatMasker
Step 2: Create informant BLAST database
Step 3: Run Blastz
Step 4: Convert blastz output to maf output
Step 5: Convert maf output to align output
(Step 5b: Create EST sequence)
Step 6: Run Twinscan

These six steps (without 5b) are all contained in the example script Nscan_driver,
which comes with this distribution. 
Nscan_driver needs a configuration file that contains the paths to all input files.
To create an example configuration file, run
Nscan_driver --config > config.file

Open the file in a text editor and change all the paths according to your system.

