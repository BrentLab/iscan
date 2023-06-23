#!/usr/bin/perl -w
$| = 1;
use strict;

my $DNA = "T/test.dna";

my @test = (
	    {cmd => "./tooltest -int T/some_int",     std => "tooltest.1"},
	    {cmd => "./tooltest -float T/some_float", std => "tooltest.2"},
	    {cmd => "./tooltest -text T/some_text",   std => "tooltest.3"},
	    {cmd => "./tooltest -table T/some_table", std => "tooltest.4"},
	    {cmd => "./tooltest -intern T/some_text", std => "tooltest.5"},
	    {cmd => "./seqstats $DNA", std => "seqstats.1"},
	    {cmd => "./seqedit -s=521 -l=120 $DNA",    std => "seqedit.1"},
	    {cmd => "./seqedit -s=521 -l=120 -x $DNA", std => "seqedit.2"},
	    {cmd => "./seqedit -s=521 -l=120 $DNA | ./seqedit -a -",  std => "seqedit.3"},
	    {cmd => "./seqedit -s=521 -l=120 $DNA | ./seqedit -a -x -", std => "seqedit.4"},
	    {cmd => "./gctest -w=5000000 T/gctest1.fa", std => "gctest.1"},
	    {cmd => "./gctest -w=50 T/gctest1.fa",      std => "gctest.2"},
	    {cmd => "./mathtest func",                  std => "mathtest.1"},
	    {cmd => "./mathtest dist T/dist.constant",  std => "mathtest.2"},
	    {cmd => "./mathtest dist T/dist.defined",   std => "mathtest.3"},
	    {cmd => "./mathtest dist T/dist.geometric", std => "mathtest.4"},
	    {cmd => "./mathtest dist T/dist.poisson",   std => "mathtest.5"},
	    {cmd => "./mathtest dura T/duration",       std => "mathtest.6"},
	    {cmd => "./mathtest add",                   std => "mathtest.8"},
	    {cmd => "./modeltest T/model.start",         std => "modeltest.1"},
	    {cmd => "./modeltest T/model.cds",           std => "modeltest.2"},
	    {cmd => "./modeltest T/model.acceptor.wwam", std => "modeltest.3"},
	    {cmd => "./scannertest T/model.start T/test.dna 519 522",      std => "scannertest.1"},
	    {cmd => "./scannertest T/model.stop T/test.dna 4243 4247",     std => "scannertest.2"},
	    {cmd => "./scannertest T/model.donor T/test.dna 640 644",      std => "scannertest.3"},
	    {cmd => "./scannertest T/model.acceptor T/test.dna 1063 1067", std => "scannertest.4"},
	    {cmd => "./scannertest T/model.iso.acceptor T/test.dna 1063 1067", std => "scannertest.5"},
	    {cmd => "./scannertest T/model.iso.acceptor T/test.enriched.dna 1063 1067", std => "scannertest.4"},
	    {cmd => "./scannertest T/model.iso.windowed T/gctest1.fa 97 111", std => "scannertest.6"},
	    {cmd => "./sftest T/test.sf", std => "sftest.1"},
	    {cmd => "./factorytest T/test.dna", std => "factorytest.1"},
	    {cmd => "./conseqtest T/GO_10000.seq T/GO_10000.conseq", std => "conseqtest.1"},
	    {cmd => "./conseqscannertest T/model.conseqscannertest T/GO_10000.conseq 300 400", std => "conseqscannertest.1"},
	    {cmd => "./gtftest -s T/alt.gtf", std => "gtftest.1"},

	    {cmd => "./iscan human.zhmm T/alt.dna -c=T/alt.icon -f=Esngl", std => "iscantest.alt.esngl"},
	    {cmd => "./iscan human.zhmm T/alt.dna -c=T/alt.icon -f=Eterm", std => "iscantest.alt.eterm"},
	    {cmd => "./iscan human.zhmm T/alt.dna -c=T/alt.icon -f=Exon", std => "iscantest.alt.exon"},
	    {cmd => "./iscan human.zhmm T/alt.dna -c=T/alt.icon -f=Einit", std => "iscantest.alt.einit"},
	    {cmd => "./iscan human.zhmm T/alt.dna -c=T/alt.icon -f=PolyA", std => "iscantest.alt.polya"},
	    {cmd => "./iscan human.zhmm T/alt.dna -c=T/alt.icon -f=Prom", std => "iscantest.alt.prom"},
	    {cmd => "./iscan human.zhmm T/alt.dna -c=T/alt.icon -i|tail +3", std => "iscantest.alt.conseq"},
	    {cmd => "./iscan human.zhmm T/alt.dna -i|tail +3", std => "iscantest.alt"},
	    {cmd => "./iscan human.zhmm T/alt.dna -t=T/partialsfeature", std => "iscantest.alt"},

	    {cmd => "./iscan human.zhmm T/GO_10000.seq -c=T/GO_10000.conseq -f=Esngl", std => "iscantest.GO_10000.esngl"},
	    {cmd => "./iscan human.zhmm T/GO_10000.seq -c=T/GO_10000.conseq -f=Eterm", std => "iscantest.GO_10000.eterm"},
	    {cmd => "./iscan human.zhmm T/GO_10000.seq -c=T/GO_10000.conseq -f=Exon",  std => "iscantest.GO_10000.exon" },
	    {cmd => "./iscan human.zhmm T/GO_10000.seq -c=T/GO_10000.conseq -f=Einit", std => "iscantest.GO_10000.einit"},
	    {cmd => "./iscan human.zhmm T/GO_10000.seq -c=T/GO_10000.conseq -f=PolyA", std => "iscantest.GO_10000.polya"},
	    {cmd => "./iscan human.zhmm T/GO_10000.seq -c=T/GO_10000.conseq -f=Prom",  std => "iscantest.GO_10000.prom" },
	    {cmd => "./iscan human.zhmm T/GO_10000.seq -c=T/GO_10000.conseq -i|tail +3",  std => "iscantest.GO_10000" }
);

my $fail = 0;
foreach my $test (@test) {
	print "$test->{cmd}  ";
	my $error = system("$test->{cmd} > zoe_test_out") | 
	    system("diff zoe_test_out T/$test->{std} > zoe_test_diff");
	print "." x (60 - length($test->{cmd}));
	if ($error || -s "zoe_test_diff" || ! -e "T/$test->{std}") {
		print " failed\n";
		$fail++;
	}
	else {
		print " passed\n";
	}
	unlink "zoe_test_out";
	unlink "zoe_test_diff";
}

print "Passed ", scalar(@test) - $fail, " / ", scalar(@test), "\n";
if ($fail > 0) {exit(1);}
