#!/usr/bin/perl -w
use strict;

my $TRUE = 1;
my $FALSE = 0;

my $START = 0;
my $STOP = 1;
my $SEQ = 2;
my $SCORE = 3;

my $count_alignment_blocks = 0;

my ($i,  $k);

my $unaligned_seq;
my $maf_file;

my @informant_seq;

my $unaligned_char = ".";
my $gap_char = "_";

##### input command line options RANDY

my $usage = "$0
        <tempdir>
        <MAF file>
        <A (ascending include overlapping alignments) | D (descending do not include lower scoring overlaps)>
        <space separated list of target and informant sequence files in order>
";

@ARGV >= 5 || die $usage;

my $TMP_PATH =  $ARGV[0];
my $input_maf_file =  $ARGV[1];
my $sort_code = $ARGV[2];
# delete rest of string
$sort_code = substr($sort_code, 0, 1);
$sort_code =~ tr/ad/AD/; 
my $seq_file =  $ARGV[3];

# check files
die "ERROR, $input_maf_file should be a valid filename\n" unless ( -f "$input_maf_file");
die "ERROR, $seq_file should be a valid filename\n" unless ( -f "$seq_file");


my @seq_names;
my %name_map;
for ($i = 3; $i < @ARGV; $i++) {
    my$id=$ARGV[$i];
    unless( -f "$id" ){
	die "ERROR, $id should be a valid filename\n";
    }
    $id =~ s/.*\///;		# remove path
    $seq_names [$i - 3] = $id;
    $name_map {$seq_names[$i - 3]} = $i - 3;
    if ($i > 3) {
        $informant_seq [$i - 1] = "";
    }
}
my $target_name = $seq_names [0];

##### read in the target sequence and translate to upper case

my $target_seq = "";
open (SEQ_FILE, "$seq_file");
<SEQ_FILE>;  #skip header
my $line_num = 0;
while (<SEQ_FILE>) {
  my $seq_line = $_;
  chomp ($seq_line);
  $target_seq .= uc($seq_line);
  $line_num++;
}
close (SEQ_FILE);
my $seq_length = length($target_seq);

my @blocks;
my $block_num = 0;

##### read in the MAF file and count lines
open (MAF_FILE, "< $input_maf_file") || die "cannot open $input_maf_file: $!\n";
my @maf_lines = <MAF_FILE>;
close MAF_FILE;

##### process maf lines
#####	check for recognized genome build names
#####	store target sequence length and start
#####	make sure target sequence matches target block 
##### 	build target for block
#####	build informants with gaps and unaligned for block

my $j = 0;
while ($j<@maf_lines) {
    if ($maf_lines[$j] =~ /^a/) {  #alignment block
        $maf_lines[$j] =~ /^a\s+score=(\S+)/;
        my $score = $1;
        $count_alignment_blocks++;
        my $start_coor;
        my $block_length;
        my @this_block;
        $j++;
        while ($maf_lines[$j] =~ /^s\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
    	    my $seq_name = $1;
	    my $start = $2;
	    my $size = $3;
	    my $strand = $4;
	    my $source_size = $5;
	    my $text = $6;

            $seq_name =~ s/.*\///;		# remove path

	    if (! exists($name_map{$seq_name})) {
	        die "Alignment block starting at $start: unrecognized sequence name $seq_name\n";
	    }

	    my @text_arr = split(//, uc($text));
	    $this_block[$name_map{$seq_name}] = \@text_arr;

	    if ($seq_name eq $target_name) {
	        $block_length = length($text);	
	        $start_coor = $start;
	        $blocks[$block_num][$START] = $start_coor;
	        $blocks[$block_num][$SCORE] = $score;
	    }
	    $j++;
        }
	   
        my @target_arr = split(//, substr($target_seq, $start_coor, $block_length));
	    	    
        for (my $u=0; $u<@seq_names; $u++) {
    	    $blocks[$block_num][$SEQ][$u] = "";
        }

        my $seq_coor = 0;
        for (my $v=0; $v<$block_length; $v++) {
            ##### ignore gaps in target sequence
	    if ($this_block[0][$v] ne "-") { 
	        ##### first check to make sure target sequence matches this block
	        if ($target_arr[$seq_coor] ne $this_block[0][$v] &&
		    $target_arr[$seq_coor] ne "N") {
		    printf STDERR "Mismatch between sequence and alignment at block $seq_coor, offset $v %s %s\n",
                      $target_arr[$seq_coor], $this_block[0][$v];
	        }
	        for (my $u=0; $u<@seq_names; $u++) {
		    ##### is there an informant alignment for this block?
		    if (defined $this_block[$u]) {
                        ##### does the informant have a gap at this position?
		        if ($this_block[$u][$v] eq "-") {
			    $blocks[$block_num][$SEQ][$u] .= $gap_char;
                        ##### take informant character from MAF and copy to alignment block
		        } else {
			    if ($u == 0) {
			        #include N's from target sequence in output file
			        $blocks[$block_num][$SEQ][$u] .= $target_arr[$seq_coor];
			    } else {
			        $blocks[$block_num][$SEQ][$u] .= $this_block[$u][$v];
			    }
		        }
                    ##### if no informant alignment, then unaligned
		    } else {  #no alignment for sequence $u in this block
		        $blocks[$block_num][$SEQ][$u] .= $unaligned_char;
		    }
 	        }
	        $seq_coor++;
	    }
        }
        $blocks[$block_num][$STOP] = $blocks[$block_num][$START] + 
          length($blocks[$block_num][$SEQ][0]) - 1;
        $block_num++;
        if ($block_num % 5000 == 0) {
	    printf STDERR "Processed $block_num blocks\n";
        }
    }
    $j++;
}

##### sort blocks by ascending stop coordinate

my @unsorted = @blocks;
if ($sort_code eq "A") {
    @blocks = sort by_score_ascending @unsorted; 
} elsif ($sort_code eq "D") {
    @blocks = sort by_score_descending @unsorted; 
} else {
    die "ERROR: illegal sort code on block sort $sort_code\n";
}
    

##### now write out the data
##### first print the header

printf ">";
for (my $i=1; $i<@seq_names; $i++) {
    printf "$seq_names[$i] ";
}
printf "\n";

##### do not print the target sequence

$target_seq =~ tr/ACGTacgtBD-FH-SU-Zbd-fh-su-z/ACGTACGTN/;
#printf "$target_seq\n";
#my $target_length = length($target_seq);

##### print the informant sequences
$unaligned_seq = $target_seq;
$unaligned_seq  =~ tr/A-Za-z/\./;
for (my $i=1; $i<@seq_names; $i++) {
    ##### initialize informant sequence to unaligned
    $informant_seq [$i ] = $target_seq;
    $informant_seq [$i ] =~ tr/A-Za-z/\./;
    my $length = 0;
    my $current_coor = 0;
    for (my $block_num = 0; $block_num < @blocks; $block_num++) {
        my $block_len = $blocks [$block_num][$STOP] - $blocks [$block_num][$START] + 1;

	if ($sort_code eq "A") {
	    substr ($informant_seq[$i], $blocks [$block_num][$START], $block_len) =
              $blocks [$block_num][$SEQ][$i];
	} else {
            if (substr ($unaligned_seq, 0, $block_len) eq substr($informant_seq[$i], $blocks[$block_num][$START], $block_len) ) {
	    substr ($informant_seq[$i], $blocks [$block_num][$START], $block_len) =
              $blocks [$block_num][$SEQ][$i];
        }
    }
	    






    }
}

##### print informants to align file
for ($i = 1; $i < @seq_names; $i++) {
    $informant_seq[$i] =~ tr/ACGTacgt_.BD-FH-SU-Zbd-fh-su-z/ACGTACGT_./;
    printf "$informant_seq[$i]\n";
}

sub by_stop_coor {
    my @a_arr = @$a;
    my @b_arr = @$b;
    return $a_arr[$STOP] <=> $b_arr[$STOP];
}

sub by_score_ascending {
    my @a_arr = @$a;
    my @b_arr = @$b;
    return $a_arr[$SCORE] <=> $b_arr[$SCORE];
}

sub by_score_descending {
    my @a_arr = @$a;
    my @b_arr = @$b;
    return $b_arr[$SCORE] <=> $a_arr[$SCORE];
}
