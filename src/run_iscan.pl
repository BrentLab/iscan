#! /usr/bin/perl -w
use strict;

my($seq_file) = @ARGV;

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

my @result = ();

if($gc_percent < 43) {
    print STDERR "using isochore 1\n";
    @result = `iscan temp_iso1.zhmm $seq_file`;
}
elsif($gc_percent < 51) {
    print STDERR "using isochore 2\n";
    @result = `iscan temp_iso2.zhmm $seq_file`;
}
elsif($gc_percent < 57) {
    print STDERR "using isochore 3\n";
    @result = `iscan temp_iso3.zhmm $seq_file`;
}
else {
    print STDERR "using isochore 4\n";
    @result = `iscan temp_iso4.zhmm $seq_file`;
}

foreach my $line (@result) {
    print $line;
}

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
