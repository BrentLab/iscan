#! /usr/bin/perl -w
use strict;
use Getopt::Std;
use vars qw($opt_p $opt_i $opt_q $opt_s);
getopts('i:p:qs:');
my $usage = "run_iscan_cons.pl [options]<mode> <seq_file> <conseq_file> [<verb>] 

run iscan_cons for <seq_file> with <conseq_file>.

options: 
    -p <p_dir> parameter file directory. 
               default is \"/bio/tmp/wei/iscan/Ce_4_7_2003\".
    -i <iscan> iscan_cons. 
               default is \"/bio/tmp/wei/iscan/intron_length/iscan_cons\".
    -s <para_suffix>  suffix for the parameter file. default is \"intron_zhmm\".
    -q         generate a check file after finishing iscan. 
";

die $usage if @ARGV < 3;

my($model, $seq_file, $conseq_file, $verb) = @ARGV;

#my $iscan = "/bio/tmp/wei/iscan/Ce_4_7_2003/data/iscan_cons";
my $iscan = "/bio/tmp/wei/iscan/intron_length/iscan_cons";
#my $iscan = "/bio/tmp/wei/iscan/old_iscan/iscan_cons";
if ($opt_i) {
    $iscan = $opt_i;
}

my $para_dir = "/bio/tmp/wei/iscan/Ce_4_7_2003";
#my $para_dir = "/bio/tmp/wei/iscan/ws82_single_gene";
#my $para_dir = "/bio/tmp/wei/iscan/intron_length";
#my $para_dir = "/bio/tmp/wei/iscan/intron_length_acceptor_pseudo_para";

if ($opt_p) { $para_dir = $opt_p; }

my $suffix = "intron_zhmm";
if ($opt_s) { $suffix = $opt_s;}

if(!defined($conseq_file)) {
    $conseq_file = "";
}

if(!defined($verb)) {
    $verb = "";
}

my %seq_hash;

read_fasta($seq_file, \%seq_hash);

my $seq_name = (keys %seq_hash)[0];

my $seq = $seq_hash{$seq_name};

my %hash;

for(my $i=0;$i<length($seq);$i++) {
    my $base = substr($seq, $i, 1);
    $hash{$base}++;

}

#foreach my $base (keys %hash) {
#    print $base, "\t", $hash{$base}, "\n";
#}

my $gc = $hash{"C"} + $hash{"G"};

my $at = $hash{"A"} + $hash{"T"} + $hash{"C"} + $hash{"G"};

my $gc_percent = 100*($gc/$at);

print STDERR "GC percent = ", $gc_percent, "%\n";

#my @result = ();
my $op="";

if($gc_percent < 32) {
    print STDERR "using isochore 1\n";
    #@result = `$iscan $verb -m=$model $para_dir/worm_6phase_iso1.$suffix $seq_file $conseq_file`;
    $op = "$iscan $verb -m=$model $para_dir/worm_6phase_iso1.$suffix $seq_file $conseq_file";
}
elsif($gc_percent < 35) {
    print STDERR "using isochore 2\n";

    #@result = `$iscan  $verb -m=$model $para_dir/worm_6phase_iso2.$suffix $seq_file $conseq_file`;
    $op = "$iscan  $verb -m=$model $para_dir/worm_6phase_iso2.$suffix $seq_file $conseq_file";

}
elsif($gc_percent < 39) {
    print STDERR "using isochore 3\n";
    #@result = `$iscan $verb -m=$model $para_dir/worm_6phase_iso3.$suffix $seq_file $conseq_file`;
    $op = "$iscan $verb -m=$model $para_dir/worm_6phase_iso3.$suffix $seq_file $conseq_file";

}
else {
    print STDERR "using isochore 4\n";
    #@result =  `$iscan  $verb -m=$model $para_dir/worm_6phase_iso4.$suffix $seq_file $conseq_file`;
    $op = "$iscan  $verb -m=$model $para_dir/worm_6phase_iso4.$suffix $seq_file $conseq_file";
}

print STDERR $op, "\n";
system($op);

my $seq_file_name = $seq_file;
if ($seq_file =~ /([^\/]+)\.[^\.]+$/) {
    $seq_file_name = $1;
}

$op = "touch $seq_file_name.check";

print STDERR "$op\n";
system("$op");

sub read_fasta {
    my($filename, $hash_ref) = @_;
    open(FASTA, $filename) || die "couldnt open fasta file $filename\n";
    chomp(my @fasta = <FASTA>);

    my $read_name;
    foreach my $line (@fasta) {
	if($line =~ /^>(.*)/) {
	    $read_name = $1;
	}
	else {
	    $hash_ref->{$read_name} .= $line;
	}
    }   
}









