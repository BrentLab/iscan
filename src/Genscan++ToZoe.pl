#! /usr/bin/perl -w
use strict;

use Getopt::Std;
use vars qw($opt_h $opt_I $opt_o $opt_s);
getopts('hI:o:s:');
my $usage  = "Genscan++ToZoe_6phase.pl [option] <seq_para> <conseq_para> <info_file>

Convert Genscan++ parameter files into iscan paramerter files. 
The result parameter files will be like 
   \${base_name}_iso\${iso_counter}.\$suffix

base_name is defined in the first line of the <info_file> 
suffix can be specified by option -s, default is \"zhmm\".

options: 
   -h     hack verion: UTR use the intergenic length distribution and 
                       intron use its own length distribution.
          default is both UTR and intron use the intergenic length distribution. 
          with option -h, the name of UTR changes from Internal to GInternal.

   -I  <intron_lengh_dir_file> intron length distribution. 
                        change the nanme of intron from Internal to Explicit.
       change the geometric length distribution to the length distribution 
       according to this file.

   -o <out_dir>  output directory. default is the current directory. 
 
   -s <suffix>   suffix for the result parameter files. default is \"zhmm\".
";



# takes our paramater file (genscan.smat) and turns it into
# a Zoe style .zhmm file
die $usage if @ARGV < 3;

#my ($infile, $outfile) = @ARGV;
my ($infile, $conseq_file, $isochore_info) = @ARGV;
my $suffix = "zhmm";
if ($opt_s) { $suffix = $opt_s; }
my $out_dir = $ENV{PWD};
if ($opt_o) { $out_dir =  $opt_o; }

open(IN, $infile) or die "can't open file $infile:$!\n";
open(CONSEQ, $conseq_file) or die "can't open file $conseq_file:$!\n";
open(INFO, $isochore_info) or die "can't open file $isochore_info:$!\n";

chomp(my @infile = <IN>);
chomp(my @conseq_file = <CONSEQ>);
chomp(my @info = <INFO>);


my $output_test = 1;

my %isochore_info;
my %InitialProbibility;       #  $InitialProbibility{$feature}{$isochore} = $value
my %TransistionMatrix;        #  $TransistionMatrix{$feature1}{$feature2}{$isochore} = $value
my %IntronLength;

my %ExonDistLength;
my %ExonLength;               #ExonLength{$table} = @list
my %OtherTables;               #OtherTables{$table} = @list

my %CodingRegion;
my %DonorSite;

my @section = ();

my @initial_probability_feature_list = qw(F+ T+ N I0+ I1+ I2+);
my @intron_length_feature_list = qw(F+ T+ N I0+ I1+ I2+);
my @transistion_matrix_feature_list = qw(Einit+,I0+ Einit+,I1+ Einit+,I2+ E0+,I0+ E0+,I1+ E0+,I2+ E1+,I0+ E1+,I1+ E1+,I2+ E2+,I0+ E2+,I1+ E2+,I2+ F+,Esngl+ I0+,Eterm+ I0-,E0- I0-,E1- I0-,E2- I1-,E0- I1-,E1- I1-,E2- I2-,E0- I2-,E1- I2-,E2- F+,Einit+ N,P+ N,A- T+,A+ Eterm-,I0- Eterm-,I1- Eterm-,I2- T-,Eterm- T-,Esngl- F+,Esngl+ I0+,Eterm+ I1+,Eterm+ I2+,Eterm+ I0-,Einit- I1-,Einit- I2-,Einit- I0+,E0+ I1+,E1+ I2+,E2+ );

my $donor_counter = 0;

foreach my $line (@infile) {

    # celan up the line
    $line =~ s/\s+$//;
    $line =~ s/^\s+//;
    if($line =~ /^<(.*)>$/) { # changing sections call update_section
	update_section($line, \@section);
	if($line =~ /ExonLength\s+(\d+)\s+(\S+)/) {
#	    print $1, "\t", $2, "\n";
	    $ExonDistLength{$2} = $1;
	}
    }
    else { # not changing sections want to keep these lines somewhere
	my $section_id;
	my @section_list;
	if(@section > 0) {
	    @section_list = split(/\s+/, $section[0]);
	    $section_id = $section_list[0];
	}

        # initial probability section
	if((@section == 2) && ($section_id eq "InitialProbability") && ($line ne "")) {
	    my @isochore_list = split(/\s+/, $section[1]);
	    my $isochore = $isochore_list[2];
	    my ($feature1,$value);
	    if($line =~ /(\S+)\s+(\S+)/) {
		($feature1,$value) = ($1,$2);
	    }
	    foreach my $feature (@initial_probability_feature_list) {
		if($feature1 eq $feature) {
		    $InitialProbibility{$feature1}{$isochore} = $value;
		}
	    }
	}
	# non-isochore specific transistion matrix section
	elsif((@section == 1) && ($section_id eq "TransitionMatrix") && ($line ne "")) {
	    my ($feature1,$feature2,$value) = ("", "", "");
	    if($line =~ /(\S+)\s+\((\S+)\s+(\S+)\)/) {
		($feature1,$feature2,$value) = ($1,$3,$2);
	    }
	    foreach my $feature (@transistion_matrix_feature_list) {
		my ($first, $second) = split(/,/, $feature);
		if(($first eq $feature1) && ($second eq $feature2)) {
		    $TransistionMatrix{$feature1}{$feature2}{"none"} = $value;
		}
	    }
	}
	# isochore specific transistion matrix section
	elsif((@section == 2) && ($section_id eq "TransitionMatrix") && ($line ne "")) {

	    my @isochore_list = split(/\s+/, $section[1]);
	    my $isochore = $isochore_list[2];

	    my ($feature1,$feature2,$value) = ("", "", "");
	    if($line =~ /(\S+)\s+\((\S+)\s+(\S+)\)/) {
		($feature1,$feature2,$value) = ($1,$3,$2);
	    }
	    foreach my $feature (@transistion_matrix_feature_list) {
		my ($first, $second) = split(/,/, $feature);
		if(($first eq $feature1) && ($second eq $feature2)) {
		      $TransistionMatrix{$feature1}{$feature2}{$isochore} = $value;
		}
	    }
	}
	# isochore specific intron length section
	elsif((@section == 2) && ($section_id eq "IntronLength") && ($line ne "")) {
	    my @isochore_list = split(/\s+/, $section[1]);
	    my $isochore = $isochore_list[2];
	    my ($feature1,$value);
	    if($line =~ /(\S+)\s+(\S+)/) {
		($feature1,$value) = ($1,$2);
	    }
	    foreach my $feature (@intron_length_feature_list) {
		if($feature1 eq $feature) {
		    $IntronLength{$feature1}{$isochore} = $value;
		}
	    }
	}
	# non-isochore specific intron length section
	elsif((@section == 1) && ($section_id eq "IntronLength") && ($line ne "")) {
	    my ($feature1,$value);
	    if($line =~ /(\S+)\s+(\S+)/) {
		($feature1,$value) = ($1,$2);
	    }
	    foreach my $feature (@intron_length_feature_list) {
		if($feature1 eq $feature) {
		    $IntronLength{$feature1}{"none"} = $value;
		}
	    }
	}
	# exon length sections
	elsif((@section == 1) && ($section_id eq "ExonLength") && ($section_list[2] eq "Einit+") && ($line ne "")) {
	    push(@{$ExonLength{"Einit+"}}, $line);
	}
	elsif((@section == 1) && ($section_id eq "ExonLength") && ($section_list[2] eq "E0+") && ($line ne "")) {
	    push(@{$ExonLength{"E0+"}}, $line);
	}
	elsif((@section == 1) && ($section_id eq "ExonLength") && ($section_list[2] eq "Eterm+") && ($line ne "")) {
	    push(@{$ExonLength{"Eterm+"}}, $line);
	}
	elsif((@section == 1) && ($section_id eq "ExonLength") && ($section_list[2] eq "Esngl+") && ($line ne "")) {
	    push(@{$ExonLength{"Esngl+"}}, $line);
	}
	# other tables sections 
	elsif((@section == 1) && ($section_id eq "PolyASignal") && ($line ne "")) {
	    	    push(@{$OtherTables{"PolyASignal"}}, $line);
	}
	elsif((@section == 1) && ($section_id eq "StartCodon") && ($line ne "")) {
	    	    push(@{$OtherTables{"StartCodon"}}, $line);
	}
	elsif((@section == 1) && ($section_id eq "StopCodon") && ($line ne "")) {
	    	    push(@{$OtherTables{"StopCodon"}}, $line);
	}
	elsif((@section == 1) && ($section_id eq "TATA") && ($line ne "")) {
	    	    push(@{$OtherTables{"TATA"}}, $line);
	}
	elsif((@section == 1) && ($section_id eq "PB_CAP") && ($line ne "")) {
	    	    push(@{$OtherTables{"PB_CAP"}}, $line);
	}
	elsif((@section == 1) && ($section_id eq "SigPept") && ($line ne "")) {
	    	    push(@{$OtherTables{"SigPept"}}, $line);
	}
	elsif((@section == 1) && ($section_id eq "AcceptorRegion") && ($line ne "")) {
	    	    push(@{$OtherTables{"AcceptorRegion"}}, $line);
	}
        # coding region section
	elsif((@section == 2) && ($section_id eq "CodingRegion") && ($line ne "")) {
	    my @isochore_list = split(/\s+/, $section[1]);
	    my $isochore = $isochore_list[2];
	    push(@{$CodingRegion{$isochore}}, $line);
	}
        # donor site section
	elsif((@section == 1) && ($section_id eq "DonorSite")) {
	    if($line eq "") {
		$donor_counter++;
	    }
	    else {
		push(@{$DonorSite{$donor_counter}}, $line);
#		print $donor_counter, "\t\t", $line, "\n";
	    }
	}
    }
    
}

my %conseq_params;
my $conseq_id;
foreach my $line (@conseq_file) {

    # celan up the line
    $line =~ s/\s+$//;
    $line =~ s/^\s+//;

    if($line =~ /[a-zA-Z]+/) {
	$conseq_id = $line;
#	print $line, "\n";
    }
    else {
	push(@{$conseq_params{$conseq_id}}, $line);
    }

}


################ print out zoe .zhmm file  ####################3

my $iso_counter = 0;

my $base_name = $info[0];

for(my $i=1;$i<@info;$i++) {

    my $pair = $info[$i];

#foreach my $pair (@info) {

    $iso_counter++;

    my($isochore, $coding_isochore) = split(/\s+/, $pair);

    open(OUT, ">$out_dir/${base_name}_iso${iso_counter}.$suffix");

#my $isochore = "43.0";
#my $coding_isochore = "100.0";

my %TransistionMatrixMunge;

# these need to be rechecked carefully before testing
$TransistionMatrixMunge{'N'}{'A-'}{'none'} = score2prob(prob2score($TransistionMatrix{'N'}{'A-'}{'none'}) + 5);
$TransistionMatrixMunge{'T+'}{'A+'}{'none'} = score2prob(prob2score($TransistionMatrix{'T+'}{'A+'}{'none'}) - 5);

$TransistionMatrixMunge{'Eterm-'}{'I0-'}{'none'} = score2prob(0.5*prob2score($TransistionMatrix{'Eterm-'}{'I0-'}{'none'}) + 
                                                                  prob2score($TransistionMatrix{'I0-'}{'Einit-'}{$isochore}));


$TransistionMatrixMunge{'Eterm-'}{'I1-'}{'none'} = score2prob(0.5*prob2score($TransistionMatrix{'Eterm-'}{'I1-'}{'none'}) + 
                                                                  prob2score($TransistionMatrix{'I1-'}{'Einit-'}{$isochore}));

$TransistionMatrixMunge{'Eterm-'}{'I2-'}{'none'} = score2prob(0.5*prob2score($TransistionMatrix{'Eterm-'}{'I2-'}{'none'}) + 
                                                                  prob2score($TransistionMatrix{'I2-'}{'Einit-'}{$isochore}));

# this one kind of makes sense
$TransistionMatrixMunge{'Einit+'}{'I0+'}{'none'} = score2prob(0.5*prob2score($TransistionMatrix{'Einit+'}{'I0+'}{'none'}) + 
                                                                  prob2score($TransistionMatrix{'F+'}{'Einit+'}{$isochore}));
$TransistionMatrixMunge{'Einit+'}{'I1+'}{'none'} = score2prob(0.5*prob2score($TransistionMatrix{'Einit+'}{'I1+'}{'none'}) + 
                                                                  prob2score($TransistionMatrix{'F+'}{'Einit+'}{$isochore}));
$TransistionMatrixMunge{'Einit+'}{'I2+'}{'none'} = score2prob(0.5*prob2score($TransistionMatrix{'Einit+'}{'I2+'}{'none'}) + 
                                                                  prob2score($TransistionMatrix{'F+'}{'Einit+'}{$isochore}));


# this one makes no sense at all.  the only alternative is to use Einit+ to I+ transistions instead of the negative strand
# but that makes even less sense!
$TransistionMatrixMunge{'I0+'}{'Eterm+'}{$isochore} = score2prob( 0.5*prob2score($TransistionMatrix{'Eterm-'}{'I0-'}{'none'}) + 
                                                                      prob2score($TransistionMatrix{'I0+'}{'Eterm+'}{$isochore}));

$TransistionMatrixMunge{'I1+'}{'Eterm+'}{$isochore} = score2prob( 0.5*prob2score($TransistionMatrix{'Eterm-'}{'I1-'}{'none'}) + 
                                                                      prob2score($TransistionMatrix{'I1+'}{'Eterm+'}{$isochore}));

$TransistionMatrixMunge{'I2+'}{'Eterm+'}{$isochore} = score2prob( 0.5*prob2score($TransistionMatrix{'Eterm-'}{'I2-'}{'none'}) + 
                                                                      prob2score($TransistionMatrix{'I2+'}{'Eterm+'}{$isochore}));


# I have given up trying to understand why twinscan uses these numbers, it looks like black magic but it works.
$TransistionMatrixMunge{'I0-'}{'Einit-'}{$isochore} = score2prob( 0.5*prob2score($TransistionMatrix{'Einit+'}{'I0+'}{'none'}) + 
                                                                      prob2score($TransistionMatrix{'I0+'}{'E0+'}{$isochore}));
$TransistionMatrixMunge{'I1-'}{'Einit-'}{$isochore} = score2prob( 0.5*prob2score($TransistionMatrix{'Einit+'}{'I1+'}{'none'}) + 
                                                                      prob2score($TransistionMatrix{'I1+'}{'E1+'}{$isochore}));
$TransistionMatrixMunge{'I2-'}{'Einit-'}{$isochore} = score2prob( 0.5*prob2score($TransistionMatrix{'Einit+'}{'I2+'}{'none'}) + 
                                                                      prob2score($TransistionMatrix{'I2+'}{'E2+'}{$isochore}));

#print $TransistionMatrix{'I0+'}{'Eterm+'}{$isochore}, "\t", prob2score($TransistionMatrix{'I0+'}{'Eterm+'}{$isochore}), "\n";
#print $TransistionMatrix{'I1+'}{'Eterm+'}{$isochore}, "\t", prob2score($TransistionMatrix{'I1+'}{'Eterm+'}{$isochore}), "\n";
#print $TransistionMatrix{'I2+'}{'Eterm+'}{$isochore}, "\t", prob2score($TransistionMatrix{'I2+'}{'Eterm+'}{$isochore}), "\n\n";


#print "Einit+ --> I+ transistions\n";
#print $TransistionMatrixMunge{'Einit+'}{'I0+'}{'none'}, "\t", prob2score($TransistionMatrixMunge{'Einit+'}{'I0+'}{'none'}), "\n";
#print $TransistionMatrixMunge{'Einit+'}{'I1+'}{'none'}, "\t", prob2score($TransistionMatrixMunge{'Einit+'}{'I1+'}{'none'}), "\n";
#print $TransistionMatrixMunge{'Einit+'}{'I2+'}{'none'}, "\t", prob2score($TransistionMatrixMunge{'Einit+'}{'I2+'}{'none'}), "\n\n";

#print "Eterm- --> I- transistions\n";
#print $TransistionMatrixMunge{'Eterm-'}{'I0-'}{'none'}, "\t", prob2score($TransistionMatrixMunge{'Eterm-'}{'I0-'}{'none'}), "\n";
#print $TransistionMatrixMunge{'Eterm-'}{'I1-'}{'none'}, "\t", prob2score($TransistionMatrixMunge{'Eterm-'}{'I1-'}{'none'}), "\n";
#print $TransistionMatrixMunge{'Eterm-'}{'I2-'}{'none'}, "\t", prob2score($TransistionMatrixMunge{'Eterm-'}{'I2-'}{'none'}), "\n\n";


#print "I+ --> Eterm+ transistions\n";
#print $TransistionMatrixMunge{'I0+'}{'Eterm+'}{$isochore}, "\t", prob2score($TransistionMatrixMunge{'I0+'}{'Eterm+'}{$isochore}), "\n";
#print $TransistionMatrixMunge{'I1+'}{'Eterm+'}{$isochore}, "\t", prob2score($TransistionMatrixMunge{'I1+'}{'Eterm+'}{$isochore}), "\n";
#print $TransistionMatrixMunge{'I2+'}{'Eterm+'}{$isochore}, "\t", prob2score($TransistionMatrixMunge{'I2+'}{'Eterm+'}{$isochore}), "\n\n";

#print "I- --> Einit- transistions\n";
#print $TransistionMatrixMunge{'I0-'}{'Einit-'}{$isochore}, "\t", prob2score($TransistionMatrixMunge{'I0-'}{'Einit-'}{$isochore}), "\n";
#print $TransistionMatrixMunge{'I1-'}{'Einit-'}{$isochore}, "\t", prob2score($TransistionMatrixMunge{'I1-'}{'Einit-'}{$isochore}), "\n";
#print $TransistionMatrixMunge{'I2-'}{'Einit-'}{$isochore}, "\t", prob2score($TransistionMatrixMunge{'I2-'}{'Einit-'}{$isochore}), "\n\n";

#print $TransistionMatrix{'F+'}{'Einit+'}{$isochore}, "\n";
#print $TransistionMatrixMunge{'N'}{'A-'}{'none'}, "\n";
#print $TransistionMatrixMunge{'T+'}{'A+'}{'none'}, "\n";

if($output_test) {
    


    print OUT "zHMM i2e3 49 134 10 20\n\n";
    
    
    
    print OUT "<STATES>\n\n";
    
# need to put internal states here, need their initial probs
    print OUT "Inter		Internal + 0   Inter   Inter  	$InitialProbibility{'N'}{$isochore}";

    if ($opt_I) { #explicit intron length distribution
	print OUT "
Intron0		Explicit + 0   Intron  Intron 	$InitialProbibility{'I0+'}{$isochore}
Intron0-	Explicit - 0   Intron  Intron 	$InitialProbibility{'I0+'}{$isochore}
Intron1		Explicit + 1   Intron  Intron   $InitialProbibility{'I1+'}{$isochore}
Intron1T	Explicit + 1T  Intron  Intron   $InitialProbibility{'I1+'}{$isochore}
Intron1-	Explicit - 1   Intron  Intron   $InitialProbibility{'I1+'}{$isochore}
Intron1T-	Explicit - 1T  Intron  Intron   $InitialProbibility{'I1+'}{$isochore}
Intron2		Explicit + 2   Intron  Intron   $InitialProbibility{'I2+'}{$isochore}
Intron2TA	Explicit + 2TA Intron  Intron   $InitialProbibility{'I2+'}{$isochore}
Intron2TG	Explicit + 2TG Intron  Intron   $InitialProbibility{'I2+'}{$isochore}
Intron2-	Explicit - 2   Intron  Intron   $InitialProbibility{'I2+'}{$isochore}
Intron2TA-	Explicit - 2TA Intron  Intron   $InitialProbibility{'I2+'}{$isochore}
Intron2TG-	Explicit - 2TG Intron  Intron   $InitialProbibility{'I2+'}{$isochore}";
    }
    else {
	print OUT "
Intron0		Internal + 0   Intron  Intron 	$InitialProbibility{'I0+'}{$isochore}
Intron0-	Internal - 0   Intron  Intron 	$InitialProbibility{'I0+'}{$isochore}
Intron1		Internal + 1   Intron  Intron   $InitialProbibility{'I1+'}{$isochore}
Intron1T	Internal + 1T  Intron  Intron   $InitialProbibility{'I1+'}{$isochore}
Intron1-	Internal - 1   Intron  Intron   $InitialProbibility{'I1+'}{$isochore}
Intron1T-	Internal - 1T  Intron  Intron   $InitialProbibility{'I1+'}{$isochore}
Intron2		Internal + 2   Intron  Intron   $InitialProbibility{'I2+'}{$isochore}
Intron2TA	Internal + 2TA Intron  Intron   $InitialProbibility{'I2+'}{$isochore}
Intron2TG	Internal + 2TG Intron  Intron   $InitialProbibility{'I2+'}{$isochore}
Intron2-	Internal - 2   Intron  Intron   $InitialProbibility{'I2+'}{$isochore}
Intron2TA-	Internal - 2TA Intron  Intron   $InitialProbibility{'I2+'}{$isochore}
Intron2TG-	Internal - 2TG Intron  Intron   $InitialProbibility{'I2+'}{$isochore}";
    }
    if ($opt_h) { #hack verion, UTR use intergenic length distribution
	print OUT "
Utr3		GInternal + 0   Utr3    Utr3     $InitialProbibility{'T+'}{$isochore}
Utr3-		GInternal - 0   Utr3    Utr3     $InitialProbibility{'T+'}{$isochore}
Utr5		GInternal + 0   Utr5    Utr5   	$InitialProbibility{'F+'}{$isochore}
Utr5-		GInternal - 0   Utr5    Utr5     $InitialProbibility{'F+'}{$isochore}";
    }
    else{

	print OUT "
Utr3		Internal + 0   Utr3    Utr3     $InitialProbibility{'T+'}{$isochore}
Utr3-		Internal - 0   Utr3    Utr3     $InitialProbibility{'T+'}{$isochore}
Utr5		Internal + 0   Utr5    Utr5   	$InitialProbibility{'F+'}{$isochore}
Utr5-		Internal - 0   Utr5    Utr5     $InitialProbibility{'F+'}{$isochore}";
    }
print OUT "
Einit0		External + 0   Einit   Einit    0.0
Einit1		External + 1   Einit   Einit    0.0
Einit1T		External + 1T  Einit   Einit    0.0
Einit2    	External + 2   Einit   Einit    0.0
Einit2TA	External + 2TA Einit   Einit    0.0
Einit2TG	External + 2TG Einit   Einit    0.0
Eterm		External + 0   Eterm   Eterm    0.0
Exon0		External + 0   Exon    Exon     0.0
Exon1		External + 1   Exon    Exon     0.0
Exon1T		External + 1T  Exon    Exon     0.0
Exon2		External + 2   Exon    Exon     0.0
Exon2TA		External + 2TA Exon    Exon     0.0
Exon2TG		External + 2TG Exon    Exon     0.0
Einit-		External - 0   Einit   Einit    0.0
Eterm0-		External - 0   Eterm   Eterm    0.0
Eterm1-		External - 1   Eterm   Eterm    0.0
Eterm1T-	External - 1T  Eterm   Eterm    0.0
Eterm2-		External - 2   Eterm   Eterm    0.0
Eterm2TA-	External - 2TA Eterm   Eterm    0.0
Eterm2TG-	External - 2TG Eterm   Eterm    0.0
Exon0-		External - 0   Exon    Exon     0.0
Exon1-		External - 1   Exon    Exon     0.0
Exon1T-		External - 1T  Exon    Exon     0.0
Exon2-		External - 2   Exon    Exon     0.0
Exon2TA-	External - 2TA Exon    Exon     0.0
Exon2TG-	External - 2TG Exon    Exon     0.0
PolyA		External + 0   PolyA   PolyA    0.0
PolyA-		External - 0   PolyA   PolyA    0.0
Prom		External + 0   Prom    Prom     0.0
Prom-		External - 0   Prom    Prom     0.0
Esngl		External + 0   Esngl   Esngl    0.0
Esngl-		External - 0   Esngl   Esngl    0.0\n";
    
    print OUT "\n\n<STATE_TRANSITIONS>\n\n";
    
    print OUT "
Utr3-   Eterm0-               $TransistionMatrix{'T-'}{'Eterm-'}{$isochore}
Utr3-   Eterm1-               $TransistionMatrix{'T-'}{'Eterm-'}{$isochore}
Utr3-   Eterm1T-              $TransistionMatrix{'T-'}{'Eterm-'}{$isochore}
Utr3-   Eterm2-               $TransistionMatrix{'T-'}{'Eterm-'}{$isochore}
Utr3-   Eterm2TA-             $TransistionMatrix{'T-'}{'Eterm-'}{$isochore}
Utr3-   Eterm2TG-             $TransistionMatrix{'T-'}{'Eterm-'}{$isochore}
Eterm0-  Intron0-             $TransistionMatrixMunge{'Eterm-'}{'I0-'}{'none'}
Eterm1-  Intron1-             $TransistionMatrixMunge{'Eterm-'}{'I1-'}{'none'}
Eterm1T- Intron1T-            $TransistionMatrixMunge{'Eterm-'}{'I1-'}{'none'}
Eterm2-  Intron2-             $TransistionMatrixMunge{'Eterm-'}{'I2-'}{'none'}
Eterm2TA-  Intron2TA-         $TransistionMatrixMunge{'Eterm-'}{'I2-'}{'none'}
Eterm2TG-  Intron2TG-         $TransistionMatrixMunge{'Eterm-'}{'I2-'}{'none'}
Inter   PolyA-                $TransistionMatrixMunge{'N'}{'A-'}{'none'}
PolyA-	Utr3-                 1.0
Intron0- Exon0-               $TransistionMatrix{'I0-'}{'E0-'}{$isochore}
Intron0- Exon1-               $TransistionMatrix{'I1-'}{'E0-'}{$isochore}
Intron0- Exon1T-              $TransistionMatrix{'I1-'}{'E0-'}{$isochore}
Intron0- Exon2-               $TransistionMatrix{'I2-'}{'E0-'}{$isochore}
Intron0- Exon2TA-             $TransistionMatrix{'I2-'}{'E0-'}{$isochore}
Intron0- Exon2TG-             $TransistionMatrix{'I2-'}{'E0-'}{$isochore}
Intron1- Exon0-               $TransistionMatrix{'I0-'}{'E1-'}{$isochore}
Intron1T- Exon0-              $TransistionMatrix{'I0-'}{'E1-'}{$isochore}
Intron1- Exon1-               $TransistionMatrix{'I1-'}{'E1-'}{$isochore}
Intron1- Exon1T-              $TransistionMatrix{'I1-'}{'E1-'}{$isochore}
Intron1T- Exon1-              $TransistionMatrix{'I1-'}{'E1-'}{$isochore}
Intron1T- Exon1T-             $TransistionMatrix{'I1-'}{'E1-'}{$isochore}
Intron1- Exon2-               $TransistionMatrix{'I2-'}{'E1-'}{$isochore}
Intron1- Exon2TA-             $TransistionMatrix{'I2-'}{'E1-'}{$isochore}
Intron1- Exon2TG-             $TransistionMatrix{'I2-'}{'E1-'}{$isochore}
Intron1T- Exon2-              $TransistionMatrix{'I2-'}{'E1-'}{$isochore}
Intron1T- Exon2TA-            $TransistionMatrix{'I2-'}{'E1-'}{$isochore}
Intron1T- Exon2TG-            $TransistionMatrix{'I2-'}{'E1-'}{$isochore}
Intron2- Exon0-               $TransistionMatrix{'I0-'}{'E2-'}{$isochore}
Intron2TA- Exon0-             $TransistionMatrix{'I0-'}{'E2-'}{$isochore}
Intron2TG- Exon0-             $TransistionMatrix{'I0-'}{'E2-'}{$isochore}
Intron2- Exon1-               $TransistionMatrix{'I1-'}{'E2-'}{$isochore}
Intron2- Exon1T-              $TransistionMatrix{'I1-'}{'E2-'}{$isochore}
Intron2TA- Exon1-             $TransistionMatrix{'I1-'}{'E2-'}{$isochore}
Intron2TA- Exon1T-            $TransistionMatrix{'I1-'}{'E2-'}{$isochore}
Intron2TG- Exon1-             $TransistionMatrix{'I1-'}{'E2-'}{$isochore}
Intron2TG- Exon1T-            $TransistionMatrix{'I1-'}{'E2-'}{$isochore}
Intron2- Exon2-               $TransistionMatrix{'I2-'}{'E2-'}{$isochore}
Intron2- Exon2TA-             $TransistionMatrix{'I2-'}{'E2-'}{$isochore}
Intron2- Exon2TG-             $TransistionMatrix{'I2-'}{'E2-'}{$isochore}
Intron2TA- Exon2-             $TransistionMatrix{'I2-'}{'E2-'}{$isochore}
Intron2TA- Exon2TA-           $TransistionMatrix{'I2-'}{'E2-'}{$isochore}
Intron2TA- Exon2TG-           $TransistionMatrix{'I2-'}{'E2-'}{$isochore}
Intron2TG- Exon2-             $TransistionMatrix{'I2-'}{'E2-'}{$isochore}
Intron2TG- Exon2TA-           $TransistionMatrix{'I2-'}{'E2-'}{$isochore}
Intron2TG- Exon2TG-           $TransistionMatrix{'I2-'}{'E2-'}{$isochore}
Intron0- Einit-               $TransistionMatrixMunge{'I0-'}{'Einit-'}{$isochore}
Intron1- Einit-               $TransistionMatrixMunge{'I1-'}{'Einit-'}{$isochore}
Intron1T- Einit-              $TransistionMatrixMunge{'I1-'}{'Einit-'}{$isochore}
Intron2- Einit-               $TransistionMatrixMunge{'I2-'}{'Einit-'}{$isochore}
Intron2TA- Einit-             $TransistionMatrixMunge{'I2-'}{'Einit-'}{$isochore}
Intron2TG- Einit-             $TransistionMatrixMunge{'I2-'}{'Einit-'}{$isochore}
Exon0-   Intron0-             1.0
Exon1-   Intron1-             1.0
Exon1T-   Intron1T-           1.0
Exon2-   Intron2-             1.0
Exon2TA-   Intron2TA-         1.0
Exon2TG-   Intron2TG-         1.0
Einit-	Utr5-	              1.0
Utr5-	Prom-	              1.0
Prom-	Inter	              1.0
Inter	Prom	              $TransistionMatrix{'N'}{'P+'}{'none'}
Prom	Utr5	              1.0
Utr5	Einit0	              $TransistionMatrix{'F+'}{'Einit+'}{$isochore}
Utr5	Einit1	              $TransistionMatrix{'F+'}{'Einit+'}{$isochore}
Utr5	Einit1T	              $TransistionMatrix{'F+'}{'Einit+'}{$isochore}
Utr5	Einit2	              $TransistionMatrix{'F+'}{'Einit+'}{$isochore}
Utr5	Einit2TA              $TransistionMatrix{'F+'}{'Einit+'}{$isochore}
Utr5	Einit2TG	      $TransistionMatrix{'F+'}{'Einit+'}{$isochore}
Intron0	Exon0	              $TransistionMatrix{'I0-'}{'E0-'}{$isochore}
Intron0	Exon1	              $TransistionMatrix{'I0-'}{'E1-'}{$isochore}
Intron0	Exon1T	              $TransistionMatrix{'I0-'}{'E1-'}{$isochore}
Intron0	Exon2	              $TransistionMatrix{'I0-'}{'E2-'}{$isochore}
Intron0	Exon2TA	              $TransistionMatrix{'I0-'}{'E2-'}{$isochore}
Intron0	Exon2TG	              $TransistionMatrix{'I0-'}{'E2-'}{$isochore}
Intron0	Eterm	              $TransistionMatrixMunge{'I0+'}{'Eterm+'}{$isochore}
Intron1	Exon0	              $TransistionMatrix{'I1-'}{'E0-'}{$isochore}
Intron1	Exon1	              $TransistionMatrix{'I1-'}{'E1-'}{$isochore}
Intron1	Exon1T	              $TransistionMatrix{'I1-'}{'E1-'}{$isochore}
Intron1	Exon2	              $TransistionMatrix{'I1-'}{'E2-'}{$isochore}
Intron1	Exon2TA	              $TransistionMatrix{'I1-'}{'E2-'}{$isochore}
Intron1	Exon2TG	              $TransistionMatrix{'I1-'}{'E2-'}{$isochore}
Intron1	Eterm	              $TransistionMatrixMunge{'I1+'}{'Eterm+'}{$isochore}
Intron1T	Exon0	      $TransistionMatrix{'I1-'}{'E0-'}{$isochore}
Intron1T	Exon1	      $TransistionMatrix{'I1-'}{'E1-'}{$isochore}
Intron1T	Exon1T	      $TransistionMatrix{'I1-'}{'E1-'}{$isochore}
Intron1T	Exon2	      $TransistionMatrix{'I1-'}{'E2-'}{$isochore}
Intron1T	Exon2TA	      $TransistionMatrix{'I1-'}{'E2-'}{$isochore}
Intron1T	Exon2TG	      $TransistionMatrix{'I1-'}{'E2-'}{$isochore}
Intron1T	Eterm	      $TransistionMatrixMunge{'I1+'}{'Eterm+'}{$isochore}
Intron2	Exon0	              $TransistionMatrix{'I2-'}{'E0-'}{$isochore}
Intron2	Exon1	              $TransistionMatrix{'I2-'}{'E1-'}{$isochore}
Intron2	Exon1T	              $TransistionMatrix{'I2-'}{'E1-'}{$isochore}
Intron2	Exon2	              $TransistionMatrix{'I2-'}{'E2-'}{$isochore}
Intron2	Exon2TA	              $TransistionMatrix{'I2-'}{'E2-'}{$isochore}
Intron2	Exon2TG	              $TransistionMatrix{'I2-'}{'E2-'}{$isochore}
Intron2	Eterm	              $TransistionMatrixMunge{'I2+'}{'Eterm+'}{$isochore}
Intron2TA	Exon0	      $TransistionMatrix{'I2-'}{'E0-'}{$isochore}
Intron2TA	Exon1	      $TransistionMatrix{'I2-'}{'E1-'}{$isochore}
Intron2TA	Exon1T	      $TransistionMatrix{'I2-'}{'E1-'}{$isochore}
Intron2TA	Exon2	      $TransistionMatrix{'I2-'}{'E2-'}{$isochore}
Intron2TA	Exon2TA	      $TransistionMatrix{'I2-'}{'E2-'}{$isochore}
Intron2TA	Exon2TG	      $TransistionMatrix{'I2-'}{'E2-'}{$isochore}
Intron2TA	Eterm	      $TransistionMatrixMunge{'I2+'}{'Eterm+'}{$isochore}
Intron2TG	Exon0	      $TransistionMatrix{'I2-'}{'E0-'}{$isochore}
Intron2TG	Exon1	      $TransistionMatrix{'I2-'}{'E1-'}{$isochore}
Intron2TG	Exon1T	      $TransistionMatrix{'I2-'}{'E1-'}{$isochore}
Intron2TG	Exon2	      $TransistionMatrix{'I2-'}{'E2-'}{$isochore}
Intron2TG	Exon2TA	      $TransistionMatrix{'I2-'}{'E2-'}{$isochore}
Intron2TG	Exon2TG	      $TransistionMatrix{'I2-'}{'E2-'}{$isochore}
Intron2TG	Eterm	      $TransistionMatrixMunge{'I2+'}{'Eterm+'}{$isochore}
Einit0	Intron0	              $TransistionMatrixMunge{'Einit+'}{'I0+'}{'none'}
Einit1	Intron1	              $TransistionMatrixMunge{'Einit+'}{'I1+'}{'none'}
Einit1T	Intron1T	      $TransistionMatrixMunge{'Einit+'}{'I1+'}{'none'}
Einit2	Intron2	              $TransistionMatrixMunge{'Einit+'}{'I2+'}{'none'}
Einit2TA	Intron2TA     $TransistionMatrixMunge{'Einit+'}{'I2+'}{'none'}
Einit2TG	Intron2TG     $TransistionMatrixMunge{'Einit+'}{'I2+'}{'none'}
Exon0	Intron0	              1.0
Exon1	Intron1	              1.0
Exon1T	Intron1T	      1.0
Exon2	Intron2	              1.0
Exon2TA	Intron2TA	      1.0
Exon2TG	Intron2TG	      1.0
Eterm	Utr3	              1.0
Utr3	PolyA	              $TransistionMatrixMunge{'T+'}{'A+'}{'none'}
PolyA	Inter	              1.0
Utr5	Esngl	              $TransistionMatrix{'F+'}{'Esngl+'}{$isochore}
Esngl	Utr3	              1.0
Utr3-	Esngl-	              $TransistionMatrix{'T-'}{'Esngl-'}{$isochore}
Esngl-	Utr5-	              1.0\n";

    print OUT "\n\n<PHASE_PREFERENCES>\n\n";

    for(my $i=0;$i<18;$i++) {
	print OUT "\t1.0\n";
    }
    
    print OUT "\n\n<STATE_DURATIONS>\n\n";
    
    
    print OUT "
Einit 2
	DEFINED 0 ", 3*$ExonDistLength{"Einit+"}-4, "\n";  ##need to finish doing these!!!!!
    
    my @einit_score = convert_exon_probs(@{$ExonLength{"Einit+"}});
    print OUT "-100\t";
    for(my $i=3;$i<@einit_score - 4;$i+=4) {
	print OUT $einit_score[$i], "\t", $einit_score[$i+1], "\t", $einit_score[$i+2], "\t", $einit_score[$i+3], "\n";
    }
    
    print OUT "
	CONSTANT ", 3*$ExonDistLength{"Einit+"}-3, " -1\n";	
    print OUT"		-300	\n";
    
    print OUT "
Eterm 2
	DEFINED 0 ", 3*$ExonDistLength{"Einit+"}-4, "\n";
    
    my @eterm_score = convert_exon_probs(@{$ExonLength{"Eterm+"}});
    print OUT "-100\t";
    for(my $i=3;$i<@eterm_score - 4;$i+=4) {
	print OUT $eterm_score[$i], "\t", $eterm_score[$i+1], "\t", $eterm_score[$i+2], "\t", $eterm_score[$i+3], "\n";
    }
    
    print OUT "
	CONSTANT ", 3*$ExonDistLength{"Einit+"}-3, " -1\n";	
    print OUT "		-300	\n";
    
    print OUT "
Exon 2
	DEFINED 0 ", 3*$ExonDistLength{"Einit+"}-4, "\n";
    
    my @exon_score = convert_exon_probs(@{$ExonLength{"E0+"}});
    print OUT "-100\t";
    for(my $i=3;$i<@exon_score - 4;$i+=4) {
	print OUT $exon_score[$i], "\t", $exon_score[$i+1], "\t", $exon_score[$i+2], "\t", $exon_score[$i+3], "\n";
    }
    
    print OUT "
	CONSTANT ", 3*$ExonDistLength{"Einit+"}-3, " -1\n";	
    print OUT "		-300	\n";
    
    print OUT "
Esngl 2
	DEFINED 0 ", 3*$ExonDistLength{"Einit+"}-4, "\n";
    
    my @esngl_score = convert_exon_probs(@{$ExonLength{"Esngl+"}});
    print OUT "-100\t";
    for(my $i=3;$i<@esngl_score - 4;$i+=4) {
	print OUT $esngl_score[$i], "\t", $esngl_score[$i+1], "\t", $esngl_score[$i+2], "\t", $esngl_score[$i+3], "\n";
    }
    
    print OUT "
	CONSTANT ", 3*$ExonDistLength{"Einit+"}-3, " -1\n";	
    print OUT "		-300	\n";

    if (!$opt_I) {
	print OUT "
Intron 1
	GEOMETRIC 0 -1	
                $IntronLength{'I0+'}{$isochore}";

    }
    else {
	my $f  = $opt_I;
	open(INTRON_LEN, $f) or die "can't open file $f:$!\n";
	my @lines = <INTRON_LEN>;
	for (my $i = 0; $i < @lines; $i ++) {
	    print OUT $lines[$i];
	}
	close(INTRON_LEN);
    }


    print OUT "

Inter 1
	GEOMETRIC 0 -1
                $IntronLength{'N'}{$isochore}

Utr3 1
	GEOMETRIC 0 -1
                $IntronLength{'T+'}{'none'}
Utr5 1
	GEOMETRIC 0 -1
                $IntronLength{'F+'}{'none'}

PolyA 2
	DEFINED 0 6
-300	-300	-300	-300	-300	-300	0

	CONSTANT 7 -1	
		-300	

Prom 2
	DEFINED 0 40

-300	-300	-300	-300	-300
-300	-300	-300	-300	-300
-300	-300	-300	-300	-300
-300	-300	-300	-300	-300
-300	-300	-300	-300	-300
-300	-300	-300	-300	-300
-300	-300	-300	-300	-300
-300	-300	-300	-300	-300	0

	CONSTANT 41 -1
		-300
\n";


    print OUT "\n\n<SEQUENCE_MODELS>\n\n";
    
    print OUT "
Acceptor SDT DNA 2 1 4 2
	AG WWAM DNA 43 42 4 0 3\n";
    
    for(my $i=0;$i<@{$OtherTables{"AcceptorRegion"}};$i++) {
	print OUT  ${$OtherTables{"AcceptorRegion"}}[$i], "\n";
    }
    
    print OUT "
	NN WMM DNA 1 0 4 0
		.	.	.	.\n";
    
    print OUT "
Donor SDT DNA 9 3 4 14
	NNNGTNNAN WMM DNA 9 3 4 0\n";
    
    for(my $i=0;$i<@{$DonorSite{"0"}};$i++) {
	print OUT ${$DonorSite{"0"}}[$i], "\n";
    }
    
    print OUT "	NNNGTNNCN WMM DNA 9 3 4 0\n";
    
    for(my $i=0;$i<@{$DonorSite{"0"}};$i++) {
	print OUT  ${$DonorSite{"0"}}[$i], "\n";
    }
    
    print OUT "	NNNGTNNTN WMM DNA 9 3 4 0\n";
    
    for(my $i=0;$i<@{$DonorSite{"0"}};$i++) {
	print OUT ${$DonorSite{"0"}}[$i], "\n";
    }
    
    print OUT "	NNAGTNNGN WMM DNA 9 3 4 0\n";
    
    for(my $i=0;$i<@{$DonorSite{"1"}};$i++) {
	print OUT  ${$DonorSite{"1"}}[$i], "\n";
    }
    
    print OUT "	NNCGTNNGN WMM DNA 9 3 4 0\n";
    
    for(my $i=0;$i<@{$DonorSite{"1"}};$i++) {
	print OUT ${$DonorSite{"1"}}[$i], "\n";
    }
    
    print OUT "	NNTGTNNGN WMM DNA 9 3 4 0\n";
    
    for(my $i=0;$i<@{$DonorSite{"1"}};$i++) {
	print OUT  ${$DonorSite{"1"}}[$i], "\n";
    }
    
    print OUT "	NCGGTNNGN WMM DNA 9 3 4 0\n";
    
    for(my $i=0;$i<@{$DonorSite{"2"}};$i++) {
	print OUT ${$DonorSite{"2"}}[$i], "\n";
    }
    
    print OUT "	NGGGTNNGN WMM DNA 9 3 4 0\n";
    
    for(my $i=0;$i<@{$DonorSite{"2"}};$i++) {
	print OUT  ${$DonorSite{"2"}}[$i], "\n";
    }
    
    print OUT "	NTGGTNNGN WMM DNA 9 3 4 0\n";
    
    for(my $i=0;$i<@{$DonorSite{"2"}};$i++) {
	print OUT  ${$DonorSite{"2"}}[$i], "\n";
    }
    
    print OUT "	NAGGTNNGA WMM DNA 9 3 4 0\n";
    
    for(my $i=0;$i<@{$DonorSite{"3"}};$i++) {
	print OUT  ${$DonorSite{"3"}}[$i], "\n";
    }
    
    print OUT "	NAGGTNNGC WMM DNA 9 3 4 0\n";
    
    for(my $i=0;$i<@{$DonorSite{"3"}};$i++) {
	print OUT ${$DonorSite{"3"}}[$i], "\n";
    }
    
    print OUT "	NAGGTNNGG WMM DNA 9 3 4 0\n";
    
    for(my $i=0;$i<@{$DonorSite{"3"}};$i++) {
	print OUT ${$DonorSite{"3"}}[$i], "\n";
    }
    
    print OUT "	NAGGTNNGT WMM DNA 9 3 4 0\n";
    
    for(my $i=0;$i<@{$DonorSite{"4"}};$i++) {
	print OUT  ${$DonorSite{"4"}}[$i], "\n";
    }
    
    print OUT "
	NNNNNNNNN WMM DNA 1 0 4 0
		.	.	.	.\n\n\n";
    
    
    print OUT "PB_CAP WMM DNA 8 3 4 0\n";
    
    for(my $i=0;$i<@{$OtherTables{"PB_CAP"}};$i++) {
	print OUT  ${$OtherTables{"PB_CAP"}}[$i], "\n";
    }
    
    print OUT "\n\n";
    
    print OUT "TATA WMM DNA 15 3 4 0\n";
    
    for(my $i=0;$i<@{$OtherTables{"TATA"}};$i++) {
	print OUT  ${$OtherTables{"TATA"}}[$i], "\n";
    }
    
    print OUT "\n\n";
    
    print OUT "SIG_PEP SIG DNA 3 2 4 0\n";
    
    
    for(my $i=0;$i<@{$OtherTables{"SigPept"}};$i++) {
	print OUT ${$OtherTables{"SigPept"}}[$i], "\n";
    }
    
    print OUT "\n\n";
    
    print OUT "PolyA WMM DNA 6 3 4 0\n";
    
    for(my $i=0;$i<@{$OtherTables{"PolyASignal"}};$i++) {
	print OUT ${$OtherTables{"PolyASignal"}}[$i], "\n";
    }
    
    print OUT "\n\n";
    
    print OUT "
Start SDT DNA 3 0 4 2
	ATG WMM DNA 12 6 4 0\n";
    
    for(my $i=0;$i<@{$OtherTables{"StartCodon"}};$i++) {
	print OUT ${$OtherTables{"StartCodon"}}[$i], "\n";
    }
    
    print OUT "
NNN WMM DNA 1 0 4 0
		.	.	.	.\n";
    
    print OUT "
Stop SDT DNA 3 0 4 4
	TAA WMM DNA 6 0 4 0\n";
    
    
    for(my $i=0;$i<@{$OtherTables{"StopCodon"}};$i++) {
	print OUT "\t\t", ${$OtherTables{"StopCodon"}}[$i], "\n";
    }
    
    
    print OUT "	TAG WMM DNA 6 0 4 0\n";
    
    for(my $i=0;$i<@{$OtherTables{"StopCodon"}};$i++) {
	print OUT "\t\t", ${$OtherTables{"StopCodon"}}[$i], "\n";
    }
    
    print OUT "	TGA WMM DNA 6 0 4 0\n";
    
    for(my $i=0;$i<@{$OtherTables{"StopCodon"}};$i++) {
	print OUT "\t\t", ${$OtherTables{"StopCodon"}}[$i], "\n";
    }
    
    print OUT "
	NNN WMM DNA 1 0 4 0
		.	.	.	.\n";
    print OUT "
Coding CDS DNA 3 2 4 3
	frame0 LUT DNA 6 5 4 0\n";
    
    for(my $i=2048;$i<3072;$i++) {
	print OUT  ${$CodingRegion{$coding_isochore}}[$i], "\n";
    }
    
    print OUT "	frame1 LUT DNA 6 5 4 0\n";
    
    for(my $i=0;$i<1024;$i++) {
	print OUT  ${$CodingRegion{$coding_isochore}}[$i], "\n";
    }
    
    print OUT "	frame2 LUT DNA 6 5 4 0\n";
    
    for(my $i=1024;$i<2048;$i++) {
	print OUT  ${$CodingRegion{$coding_isochore}}[$i], "\n";
    }
    
    print OUT "Intron LUT DNA 6 5 4 0\n";
    
    for(my $i=0;$i<1024;$i++) {
	print OUT "0\t0\t0\t0\n";
    }
    
    print OUT "\n\n";
    
    print OUT "Inter LUT DNA 6 5 4 0\n";
    
    for(my $i=0;$i<1024;$i++) {
	print OUT "0\t0\t0\t0\n";
    }
    
    print OUT "\n\n";
    
    print OUT "Utr3 LUT DNA 6 5 4 0\n";
    
    for(my $i=0;$i<1024;$i++) {
	print OUT "0\t0\t0\t0\n";
    }
    
    print OUT "\n\n";
    
    print OUT "Utr5 LUT DNA 6 5 4 0\n";
    
    for(my $i=0;$i<1024;$i++) {
	print OUT "0\t0\t0\t0\n";
    }
    
    print OUT "\n\n";
    
    print OUT "CDSCONS LUT CONSEQ 6 5 3 0\n";
    
    for(my $i=0;$i<@{$conseq_params{"CdsConsScore"}};$i++) {
	print OUT ${conseq_params{"CdsConsScore"}}[$i], "\n";
    }
    
    print OUT "\n\n";
    
    print OUT "INTRONCONS LUT CONSEQ 6 5 3 0\n";
    
    for(my $i=0;$i<@{$conseq_params{"IntronConsScore"}};$i++) {
	print OUT ${conseq_params{"IntronConsScore"}}[$i], "\n";
    }print OUT "\n\n";
    
    print OUT "INITCONS LUT CONSEQ 6 5 3 0\n";
    
    for(my $i=0;$i<@{$conseq_params{"InitConsScore"}};$i++) {
	print OUT ${conseq_params{"InitConsScore"}}[$i], "\n";
    }print OUT "\n\n";
    
    print OUT "TERMCONS LUT CONSEQ 6 5 3 0\n";
    
    for(my $i=0;$i<@{$conseq_params{"TermConsScore"}};$i++) {
	print OUT ${conseq_params{"TermConsScore"}}[$i], "\n";
    }print OUT "\n\n";
    
    print OUT "ACCCONS WWAM CONSEQ 43 41 3 0 2\n";
    
    for(my $i=0;$i<@{$conseq_params{"AccConsScore"}};$i++) {
	print OUT ${conseq_params{"AccConsScore"}}[$i], "\n";
    }print OUT "\n\n";
    
    print OUT "DONCONS WWAM CONSEQ 9 5 3 0 2\n";
    
    for(my $i=0;$i<@{$conseq_params{"DonorConsScore"}};$i++) {
	print OUT ${conseq_params{"DonorConsScore"}}[$i], "\n";
    }print OUT "\n\n";
    
    print OUT "UTRCONS LUT CONSEQ 6 5 3 0\n";
    
    for(my $i=0;$i<@{$conseq_params{"UtrConsScore"}};$i++) {
	print OUT ${conseq_params{"UtrConsScore"}}[$i], "\n";
    }
    
}

}

sub update_section {
    my ($line, $section_ref) = @_;
    
    my @line_list = split(/\s+/, $1);
    if($line_list[0] ne "End") {
	push(@section, $1);
    }
    else {
	pop(@section);
    }
}

sub convert_exon_probs {

    my @file = @_;
    
    my @score_list = ();
    
    foreach my $line (@file) {
	my @probs = split(/\s+/, $line);
	foreach my $prob (@probs) {
	    
	    my $score;
	    
	    if ($prob != 0) {
		$score = 10 * log($prob)/log(2);
	    }
	    else {
		$score = -1000.1;
	    }
	    
	    push(@score_list, $score);
	    push(@score_list, $score);
	    push(@score_list, $score);
	}
	
    }
 return @score_list;   
}

sub  prob2score {

    my($num)=@_;

    if($num != 0) {
	return 10*log($num)/log(2);	
    }
    else {
	return (-300);
    }
}

sub score2prob {

    my($num) = @_;
    
    return 2**($num/10);
    
}
