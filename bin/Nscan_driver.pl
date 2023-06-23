#! /usr/bin/perl -w 

use strict;
use FAlite;
use File::Copy;

use Getopt::Long;
my ($opt_d, $opt_m, $opt_c, $opt_h, $opt_lcmask,
    $opt_blastz, $opt_no_interspersed, $opt_nomask, $opt_add_lowcomplex, $opt_config);
GetOptions(
           "d=s" => \$opt_d,
           "m=s" => \$opt_m,
           "h" => \$opt_h,
           "c" => \$opt_c,
           "no_interspersed" => \$opt_no_interspersed,
           "nomask" => \$opt_nomask,
           "lcmask" => \$opt_lcmask,
           "add_lowcomplex" => \$opt_add_lowcomplex,
           "config" => \$opt_config,
           "blastz=s" => \$opt_blastz);

my$usage = "$0 [options] <target sequence> <configuration file> 

Driver to run Nscan.

The configuration file should contain the locations of the informant genome, the parameter file and all programs called by the driver.
To see an example configuration file, run $0 --config
Program starts with masking sequence, then runs Blastz, converts the blastz output and runs nscan.
The output is a gtf file on STDOUT, be sure to redirect this to a file.

The output is a gtf file.
Options:
  -h                                   Print this usage statement
  -m  <string> RepeatMasker species    Species can be human, mouse, rattus, mammal, carnivore, rodentia, rat, cow, pig, cat, dog, chicken, fugu,
                                       danio, ciona intestinalis drosophila, anopheles, elegans, diatoaea, artiodactyl, arabidopsis, rice, wheat, and maize
  --lcmask                             Treat lowercase as masked.
  --nomask                             Sequence has already been N masked
  --add_lowcomplex                     Mask the low-complexity DNA or simple repeats
  --no_interspersed                    Do not mask the sequence for interspersed repeats
  --config                             Show example configuration file 
  --blastz <align file>	               Do not run blastz, instead use this input file	
  -d                                   Specify Directory for output file [default is current directory]
  -c                                   Just output command list, do not run

";

# do not run blastz (and if so, do not check program!)
# Check if outputfiles already exist 

if($opt_h){
	print $usage;
	exit(0);
}

if($opt_config){
	print_config();
	exit(0);
}

my($seq, $conf) = @ARGV;
die $usage unless $conf;

our@commands = ();
our@troubles = ();
our%conf = ();

# get all file and program locations from the configuration file
parse_config_file($conf);

my$tmpdir = make_tmpdir($conf{'tmpdir'});

# check outputdir - create if it doesn't exist
my$outdir = check_outputdir($opt_d);

# check the input sequence and put a masked copy in the temporary directory
my ($base, $tname, $targetseq) = check_sequence_file ($seq, $tmpdir, $opt_lcmask, $opt_c);

# create the commands for running the pipeline
# repeatMasker
my $maskseq = run_masking ($targetseq, $tmpdir, $opt_add_lowcomplex, $opt_no_interspersed, $opt_nomask, $opt_m);

# blastz
my $alignfile = run_blastz ($maskseq, $conf{'informant'}, $opt_blastz, $tmpdir);

# Nscan
my $gtf = run_Nscan ($conf{'nscan'}, $alignfile, $maskseq, $tmpdir);

# test if output files already exist
test_outputfiles($maskseq, $alignfile, $gtf, $opt_blastz, $opt_nomask, $outdir);

if ($troubles[0]){
	print STDERR "\n Cannot run program, please correct errors\n";
	print STDERR @troubles;
	print STDERR "\n";
	unlink glob ("$tmpdir/*");
	rmdir "$tmpdir";
	exit(0);
}

# opt_c is do not run, just list
if ($opt_c){
	for my$program(@commands){
		print STDERR "Running\n$program\n\n";
	}
	rmdir "$tmpdir";
	exit(0);
}

# now run the actual commands
for my$program(@commands){
	print STDERR "Running\n$program\n\n";
	if($opt_c){
		exit();
	}
        my $retvalue = system $program;
        if ($retvalue) {
		print STDERR "ERROR: The following command died with exit code ($retvalue):\n$program\n\n";
		print STDERR "You can find your intermediate files in $tmpdir.\nDon't forget to delete them!\n";
		exit(1);
        }
}

# move files to output directory
	move_files($outdir, $gtf, $base, $tname, $maskseq, $alignfile);

sub move_files{
	my($outdir, $gtf, $base, $tname, $mask, $align) = (@_);

# the gtf file must be cleaned: replace the filename with the (truncated)
# definition line of the input fasta
	$tname =~ s/^.*\///;	# remove directory name, if present
	open(GTF, "$gtf");
	open (OUT, ">$outdir/$base.masked.gtf") || die "cannot open $outdir/$base.masked.gtf: $!";
	while(<GTF>){
		$_ =~ s/($base.*?.masked)/$tname/g;
		print OUT $_;
	}
	close OUT;
	close GTF;
	copy("$mask", "$outdir")  || die "ERROR, cannot copy $mask to $outdir: $!";
	unless($opt_blastz){
		copy("$align", "$outdir")  || die "ERROR, cannot copy $align to $outdir: $!";
	}
}

# clean up
unlink glob ("$tmpdir/*");
rmdir "$tmpdir";

exit(0);

############## SUBROUTINES ########################

sub test_outputfiles{
	my($mask, $align, $gtf, $opt_blastz, $opt_nomask, $outdir) = (@_);
	# all files contain the path to the temporary directory
	$mask =~ s/.*\///;
	$align =~ s/.*\///;
	$gtf =~ s/.*\///;
	if( -e "$outdir/$mask" && !$opt_nomask){
		push(@troubles, "Mask file $outdir/$mask already exists, remove or set --nomask flag to run Nscan\n"); 
	}
	if( -e "$outdir/$align" && !$opt_blastz){
		push(@troubles, "Align file $outdir/$align already exists, remove or set --blastz flag to run Nscan\n"); 
	}
	if( -e "$outdir/$gtf"){
		push(@troubles, "$outdir already contains a Nscan outputfile $gtf, remove it to rerun  Nscan\n"); 
	}
}


sub check_outputdir{
	(my$dir) = (@_);
	unless ($dir){
		return './';
	}
	if ( -w $opt_d ){
		return $opt_d;
	}else{
		mkdir $opt_d || die "cannot create $opt_d: $!\n";
		return $opt_d;
	}
}

sub config_error{
	my($missing) = (@_);
	push(@troubles, "Cannot find $missing in your config file $conf. Please add the location of $missing (including the filename) to the file\n");
}

sub error_log{
	my($complaint) = (@_);
	push(@troubles, $complaint);

}

# make_tmpdir: create a temporary directory under tmpdir
# seed rand (for making unique filename)

sub make_tmpdir{
	my($tmpdir) = (@_);
	die "no tmpdir in inputfile" unless $tmpdir;
	srand(time() ^ ($$ + ($$ << 15)));
	my$nr = rand;
	$nr =~ s/^0\.//;
	$tmpdir .= "tmp.$nr";
	die ("Unique dir $tmpdir exists") if ( -d "$tmpdir");
	mkdir "$tmpdir" || error_log "cannot create temporary directory $tmpdir: $!\n";
	return $tmpdir;
}



sub run_Nscan{
	my ($nscan_program, $alignfile, $targetseq, $tmpdir) = (@_);
	my $parameterfile = $conf{'parameter_file'};
	config_error("nscan") unless (-f "$nscan_program");
	config_error("parameter_file") unless ( -f "$parameterfile");
	my $nscan_cmd = "$nscan_program -o $parameterfile $targetseq -a=$alignfile > $targetseq.gtf";
	push(@commands, $nscan_cmd);
	return "$targetseq.gtf";
}


sub run_blastz{
	my ($targetseq, $informantseq, $opt_blastz, $tmpdir) = (@_);
	config_error("informant") unless (-f "$informantseq");

        open (FAFILE, "$informantseq");
        my $fasta = new FAlite(\*FAFILE);
	my$info_name;
        while(my $entry = $fasta->nextEntry) {
                $info_name = $entry->def;
		last;
	}

	if($opt_blastz){
		# check if it exists
		unless ( -e $opt_blastz ){
			push (@troubles, "Cannot find Blastz align file $opt_blastz, please add file or remove --blatz option");
		}
		return $opt_blastz;	# the value of the option is the blastz file
	}

	my $blastz_program = $conf{'blastz'};
	config_error("blastz") unless ( -f "$blastz_program");
	my $blastz_command = "$blastz_program $targetseq $informantseq K=2200 > $targetseq.lav";

	my $lav2maf = $conf{'lav2maf'};
	config_error("lav2maf") unless ( -f "$lav2maf");
	my $lav2maf_command = "$lav2maf $targetseq.lav $targetseq $informantseq > $targetseq.maf";

	my $alignfile = "$informantseq.align";
	my $maf2align = $conf{'maf2align'};
	config_error("maf2align") unless ( -f "$maf2align");
	my $maf2align_command = "$maf2align $tmpdir $targetseq.maf A $targetseq $informantseq  > $targetseq.align";

	push (@commands, $blastz_command, $lav2maf_command, $maf2align_command);
	return "$targetseq.align";
}

# check_sequence_file checks if the input file is a fasta sequence
# it then runs RepeatMasker based on the input options

sub check_sequence_file{
	my($targetseq, $tmpdir, $lcmask, $opt_c) = (@_);
	my($def, $seq);
	(my $base = $targetseq) =~ s/.*\///;
	my $copy_sequence = "$tmpdir/$base";
	unless (open (FAFILE, "$targetseq")){
		error_log ("cannot open sequence file $targetseq: $!\n");
		return (undef, undef, $copy_sequence);
	}
	my $fasta = new FAlite(\*FAFILE);
	while(my $entry = $fasta->nextEntry) {
		$def = $entry->def;
		unless ($def) {
			push (@troubles, "There does not seem to be a definition line (starting with > ) in the sequence file. Please correct.\n");
		}
		$seq = $entry->seq;
		unless ($seq) {
			push (@troubles, "There does not seem to be a sequence in the sequence file. Please correct.\n");
			$seq="";
		}
		if ($seq && $lcmask && !$opt_c){
			print STDERR "lowercase masking sequence...\n";
			$seq =~ tr/actg/NNNN/;
		}
	}
	close FAFILE;
# keep target name for output gtf, but truncate if necessary
	$def =~ s/>//;
	$def =~ s/\s.*//;
	$def =~ s/(.{10})(.*)/$1/;
# definition line should be the same as file name
	unless($opt_c){
		open (WORKSEQ, ">$copy_sequence");
		print WORKSEQ ">$targetseq\n$seq\n";
		close WORKSEQ;
	}
	return ($base, $def, $copy_sequence);
}

sub run_masking{
	my($targetseq, $tmpdir, $add_lowcomplex, $no_interspersed, $nomask, $species) = (@_);
	if($nomask){
		print STDERR "not repeat masking sequence\n";
		return($targetseq);		
	}else{
		my $repMask_program = $conf{'repmask'};
		config_error("repmask") unless ( -f "$repMask_program");

                my $repMask_command = "$repMask_program -dir $tmpdir";
		unless ($add_lowcomplex){
			$repMask_command .= " -nolow";
		}
		if ($no_interspersed){
			$repMask_command .= " -noint";
		}
		if ($species){
			$repMask_command .= " -species $species";
		}else{
			$repMask_command .= " -species human";
		}
		$repMask_command .= " $targetseq";
# If no repeats are found, a file is not created. In that case, there should be a link to the original sequence
		(my $base = $targetseq) =~ s/.*\///;
		my$checkFile_command = "if [ ! -f $targetseq.masked ]; then ln -s $base $targetseq.masked; fi";

		push(@commands, $repMask_command, $checkFile_command);
	}
	
	return "$targetseq.masked";

}

sub parse_config_file{
	my($file) = (@_);
	open(FILE, "$file") || die "cannot open config file $file: $!\n";
	while(<FILE>){
		if($_ =~ /^>/){
			print "Looks like you switched the configuration and sequence files!\n";
			exit;
		}
		next if ($_ =~ /^#/);
		next unless ($_ =~ /=/);
		chomp $_;
		my($program, $location) = split('=', $_);
		unless ( -e "$location"){
			push (@troubles, "$location in $file does not seem to exist\n");
		}
		$conf{$program}=$location;
	}
	close FILE;
}


sub print_config{

my$text = <<TOEND

# This file should contain the locations of all programs needed by Nscan_driver.pl
# You can redirect this output to a file and use it as a config file.
# Make sure you set all locations according to your system!

# directory for storage of temporary files (must be large enough to store informant genome)
tmpdir=/tmp/

# Nscan_executable
# N-SCAN can be downloaded from http://mblab.wustl.edu/software/
nscan=/home/myscripts/iscan

# RepeatMasker program
# RepeatMasker can be downloaded from http://www.repeatmasker.org/RMDownload.html
# If your sequence is already masked or you do not want to mask, 
# you can run Nscan_driver.pl with the --nomask option
# the repmask path will be ignored
repmask=/bio/bin/RepeatMasker

# Blastz program
# Blastz and multiz can be downloaded from http://www.bx.psu.edu/miller_lab/
# If you already have an ALIGN format file, 
# you can run Nscan_driver.pl with the --blastz=[align file] option
# in that case, you can remove the following three paths, or leave them:
# they will be ignored.
blastz=/usr/bin/blastz

# Path to blastz output formatter. This program comes with Multiz
lav2maf=/bio/bin/lav2maf

# the maf2align script is included in the N-SCAN download
maf2align=/home/myscripts/maf_to_align.pl

# Informant sequence file. 
# This file should be in (multiple) fasta format with Blastz headers: 
# name, chromosome, startpos, strand, and endpos separated by colons
# Example:
# >mm8:chr1:1:+:197069962

informant=/bio/myFasta/mouseGenome.fa

# parameterfile
parameter_file=/home/myParameters/param.zhmm

TOEND
;
	print $text;
	return();
}
