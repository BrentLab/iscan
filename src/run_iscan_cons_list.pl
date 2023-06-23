#!/usr/bin/perl
use strict;
use Getopt::Std;
use vars qw($opt_d $opt_c $opt_i $opt_n $opt_o  $opt_p $opt_q $opt_Q $opt_s);

my $usage = "run_iscan_cons_list.pl [options] <file_list>
run run_iscan_cons.pl for a list of files

     options: 
        -i <iscan_cons>   iscan_cons. default is \"./iscan_cons\".
        -p <para_dir>     directory contains the parameter files. 
                          default is the current directory. 
        -s <para_suffix>  suffix for the parameter files. default is \"intron_zhmm\".

        -d <seq_dir>      directory contains the sequence file. default is the current directory. 
        -c <conseq_dir>   directory contains the conseq file. default is the current directory.
        -o <out_dir>      directory contains the output result. default is the current directory.

        -q <que_file>    default is \"iscan_worm.que\";
        -Q <queue name>  default is research queue.
        -n <job_num>     divided into <n> jobs run in the queue system.
";

getopts('d:c:i:n:o:p:q:Q:s:');

die $usage if @ARGV < 1;
my $PWD = $ENV{PWD};
my $iscan = "$PWD/iscan_cons";
if ($opt_i) { $iscan = $opt_i; }

my $suffix = "intron_zhmm";
if ($opt_s) { $suffix = $opt_s; }

my $para_dir = $ENV{PWD};
if ($opt_p) { $para_dir = $opt_p; }

my $data_dir = $ENV{PWD};
if ($opt_d) { $data_dir = $opt_d; }

my $conseq_dir = $ENV{PWD};
if ($opt_c) { $conseq_dir = $opt_c; }

my $out_dir = $ENV{PWD};
if ($opt_o) { $out_dir = $opt_o; }

my $que = "iscan_worm.que";
if ($opt_q) {
    $que = $opt_q;
    open(QUE, ">$que") or die "can't open file $que:$!\n"; 
}

my $file = shift;
open(DATA, $file ) or die "can't open file $file:$!\n";
while (<DATA>) {
    chomp;
    if (/^\s*$/ or /^\$\$/) {
	next;
    }
    else {
	my $seq = $_;
	if ($seq =~ /^(\S+)/) {
	    $seq = $1;
	}
	if ($seq =~ /\/([^\/]+)\s*$/ ) {
	    $seq = $1;
	}
	my $name  = $seq;
	print STDERR "seq: $seq\n";
	if ($seq =~ /^(\S+)\.seq/ ){
	    $name = $1;
	}
	elsif ($seq =~ /^(\S+)\.fasta/ ){
	    $name = $1;
	}

	my $op = "perl /bio/tmp/wei/iscan/intron_length_acceptor_pseudo_para/run_iscan_cons.pl ".
	    "-p $para_dir -i $iscan -s $suffix ";
	$op .= " -q " if $opt_q;

	$op.=" cons  $data_dir/$seq $conseq_dir/$name.icon > $out_dir/$name.iscan_gtf ";
	print STDERR "=========$op\n";
	if ($opt_q) {
	    print QUE $op, "\n";
	}
	else {
	    system($op);
	}
    }
}


if ($opt_q) {

    close(QUE);
    my $op = "enqueue.pl $que iscan_worm.jobq $out_dir ";

    if ($opt_n) {
	 $op = "enqueue.pl -n $opt_n $que iscan_worm.jobq $out_dir ";
    }

    if ($opt_Q) {
	$op.= " -q $opt_Q ";
    }
    else {
	$op .= "  -q production ";
    }

    system($op);
}




