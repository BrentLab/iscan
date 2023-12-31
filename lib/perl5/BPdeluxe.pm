package BPdeluxe;
use strict;
use overload '""' => '_overload';
###############################################################################
# BPdeluxe
###############################################################################
sub new {
	my ($class, $fh) = @_;
	if (ref $fh !~ /GLOB/)
		{die "BPdeluxe error: new expects a GLOB reference not $fh\n"}
	my $this = bless {};
	$this->{FH} = $fh;
	$this->{LASTLINE} = "";
	if ($this->_parseHeader) {
		# there are alignments
		$this->{REPORT_DONE} = 0;
		die if not defined $this->{QUERY};
		die if not defined $this->{DATABASE};
	}
	else {
		# empty report
		$this->{REPORT_DONE} = 1; 
	}
	

	return $this;
}
sub query    {shift->{QUERY}}
sub queryLength {shift->{QUERY_LENGTH}}
sub database {shift->{DATABASE}}
sub gi       {
	my ($this) = @_;
	my $gi;
	if ($this->query =~ /gi\|(\d+)/) { $gi = $1;}
	elsif ($this->query =~ /^>\s*(\S+)/) { $gi = $1;}
	else { $gi = '0';}
}

sub length   {  #added jg
	my ($this) = @_;
	my $len    = ($this->query =~ /(\d+,*\d*,*\d*,*\d*)\s+letters/) ? $1 : '1';
	$len =~ s/,//g;
	return $len;

}
sub nextSbjct {
	my ($this) = @_;
	$this->_fastForward or return 0;
	
	#######################
	# get all sbjct lines #
	#######################
	my $def = $this->{LASTLINE};
	my $FH = $this->{FH};
	while(<$FH>) {
		if    ($_ !~ /\w/)            {next}
		elsif ($_ =~ /Strand HSP/)    {next} # WU-BLAST non-data
		elsif ($_ =~ /^\s{0,2}Score/) {$this->{LASTLINE} = $_; last}
		else                          {$def .= $_}
	}
	return 0 unless $def =~ /^>/;
	$def =~ s/\s+/ /g;
	$def =~ s/\s+$//g;
	my ($sbjct_length) = $def =~ /Length = ([\d,]+)$/;
	$sbjct_length =~ s/,//g;
	$def =~ s/Length = [\d,]+$//g;
	
	####################
	# the Sbjct object #
	####################
	my $sbjct = BPdeluxe::Sbjct::new($def, $sbjct_length,
		$this->{FH}, $this->{LASTLINE}, $this);
	return $sbjct;
}
sub _parseHeader {
	my ($this) = @_;
	my $FH = $this->{FH};
	while(<$FH>) {
		if ($_ =~ /^Query=\s+(.+)/)    {
			my $query = $1;
			while(<$FH>) {
				last if $_ !~ /\S/;
				$query .= $_;
			}
			$query =~ s/\s+/ /g;
			$this->{QUERY} = $query;
			($this->{QUERY_LENGTH}) = $query =~ /([\d,]+) letters/;
			$this->{QUERY_LENGTH} =~ s/\D//g;
		}
		elsif ($_ =~ /(\d+) letters/) {$this->{LENGTH} = $1}  #jg added
		elsif ($_ =~ /^Database:\s+(.+)/) {$this->{DATABASE} = $1}
		elsif ($_ =~ /^>/)                {$this->{LASTLINE} = $_; return 1}
		elsif ($_ =~ /^Parameters|^\s+Database:/) {
			$this->{LASTLINE} = $_;
			return 0; # there's nothing in the report
		}
	}
}
sub _fastForward {
	my ($this) = @_;
	return 0 if $this->{REPORT_DONE};
	return 1 if $this->{LASTLINE} =~ /^>/;
	return 0 if $this->{LASTLINE} =~ /^Parameters|^\s+Database:/;
	my $FH = $this->{FH};
	while(<$FH>) {
		if ($_ =~ /^>|^Parameters|^\s+Database:/) {
			$this->{LASTLINE} = $_;
			return 1;
		}
	}
	warn "Possible parse error in _fastForward in BPdeluxe.pm\n";
}
sub _overload {
	my ($this) = @_;
	return $this->{QUERY} . " vs. " . $this->{DATABASE};
}

###############################################################################
# BPdeluxe::Sbjct
###############################################################################
package BPdeluxe::Sbjct;
#use overload '""' => 'name';
sub new {
	my $sbjct = bless {};
	($sbjct->{NAME}, $sbjct->{LENGTH}, $sbjct->{FH},$sbjct->{LASTLINE},
	$sbjct->{PARENT}) = @_; $sbjct->{HSP_ALL_PARSED} = 0;
	return $sbjct;
}
sub name {shift->{NAME}}
sub length {shift->{LENGTH}}
sub nextHSP {
	my ($sbjct) = @_;
	return 0 if $sbjct->{HSP_ALL_PARSED};
	
	############################
	# get and parse scorelines #
	############################
	my $scoreline = $sbjct->{LASTLINE};
	my $FH = $sbjct->{FH};
	my $nextline = <$FH>;
	return undef if not defined $nextline;
	$scoreline .= $nextline;
	my ($score, $bits);
	if ($scoreline =~ /\d bits\)/) {
		($score, $bits) = $scoreline =~
			/Score = (\d+) \((\S+) bits\)/; # WU-BLAST
	}
	else {
		($bits, $score) = $scoreline =~
			/Score =\s+(\S+) bits \((\d+)/; # NCBI-BLAST
	}
	
	my ($match, $length) = $scoreline =~ /Identities = (\d+)\/(\d+)/;
	my ($positive) = $scoreline =~ /Positives = (\d+)/;
	$positive = $match if not defined $positive;
	my ($p,$group) = $scoreline =~ /[Sum ]*P[\(\d+\)]* = (\S+), Group = (\d+)/;  	## Deluxe Code ##
	
	my ($e) = $scoreline =~ /Expect = (\S+)/; ##miao's code ##	
	my ($strand) = $scoreline=~ /Strand = (\S+) \/ \S+/; ##miao's code ##
	
	if (not defined $p)     {$p     = $scoreline =~ /Expect = (\S+)/}
	if (not defined $group) {$group = '0';} ##miao's code ##
	$p =~ s/,//g;
	$e =~ s/,//g;
	### Deluxe code JG###
	my $sP; #JG
	if ($scoreline =~ /Sum P[\(\d+\)] = \S+/) { #JG
		$sP = $1; #JG
	} else { $sP = 1;} #JG
	
	die "parse error $scoreline\n" if not defined $score;

	#######################
	# get alignment lines #
	#######################
	my @hspline;
	while(<$FH>) {
		if ($_ =~ /^WARNING:|^NOTE:/) {
			while(<$FH>) {last if $_ !~ /\S/}
		}
		elsif ($_ !~ /\S/)            {next}
		elsif ($_ =~ /Strand HSP/)    {next} # WU-BLAST non-data
		elsif ($_ =~ /^\s*Frame/)     {next} # NCBI-BLAST non-data
		elsif ($_ =~ /^\s*Strand/)    {next} # NCBI-BLAST non-data
		elsif ($_ =~ /^\s*Score/)     {$sbjct->{LASTLINE} = $_; last}
		elsif ($_ =~ /^>|^Parameters|^\s+Database:/)   {
			$sbjct->{LASTLINE} = $_;
			$sbjct->{PARENT}->{LASTLINE} = $_;
			$sbjct->{HSP_ALL_PARSED} = 1;
			last;
		}
		else {
			push @hspline, $_;           # store the query line
			my $l1 = <$FH>;              # either alignment line or sbjct line
			if ($l1 =~ /^Sbjct/) {
				push @hspline, "";  # dummy line, this is a -noseq option
				push @hspline, $l1; # so store a fake alignment and real sbjct
				next;
			}
			push @hspline, $l1;                 # grab/store the alignment line
			my $l2 = <$FH>; push @hspline, $l2; # grab/store the sbjct line
		}
	}
	
	#########################
	# parse alignment lines #
	#########################
	my ($ql, $sl, $as) = ("", "", "");
	my ($qb, $qe, $sb, $se) = (0,0,0,0);
	my (@QL, @SL, @AS); # for better memory management
			
	for(my $i=0;$i<@hspline;$i+=3) {
		#warn $hspline[$i], $hspline[$i+2];
		$hspline[$i]   =~ /^Query:\s+(\d+)\s+(\S+)\s+(\d+)/;
		$ql = $2; $qb = $1 unless $qb; $qe = $3;
		
		my $offset = index($hspline[$i], $ql);
		$as = substr($hspline[$i+1], $offset, CORE::length($ql))
			if $hspline[$i+1];
		
		$hspline[$i+2] =~ /^Sbjct:\s+(\d+)\s+(\S+)\s+(\d+)/;
		$sl = $2; $sb = $1 unless $sb; $se = $3;
		
		push @QL, $ql; push @SL, $sl; push @AS, $as;
	}

	##################
	# the HSP object #
	##################
	$ql = join("", @QL);
	$sl = join("", @SL);
	$as = join("", @AS);
	my $qgaps = $ql =~ tr/-/-/;
	my $sgaps = $sl =~ tr/-/-/;
	my $hsp = BPdeluxe::HSP::new( $score, $bits, $match, $positive, $length, $p, $e, $strand,
		$qb, $qe, $sb, $se, $ql, $sl, $as, $qgaps, $sgaps, $group, $sP);  	## Deluxe Code ##
	return $hsp;
}

sub group_list { ### DELUXE CODE jg ###
	#purpose: to return a very informative string for a subject found in a BLAST query
	my $sbjct = shift;
	my $query = @_ ? shift : '0';  #this allows you to hand it the query object if you like
	my $gi;
	
#	print "QL = ", $query->length , "\n";
	if ($sbjct->name =~ /gi\|(\d+)/) { $gi = $1;}
	elsif ($sbjct->name =~ /^>\s*(\S+)/) { $gi = $1;}
	else { die "no identifier in defline $!";}
	
	my $string;
	my @groups = [];
	my @TOT_COV = ();
	my @TOT_MAT = ();
	
	while(my $hsp = $sbjct->nextHSP) {  # divy up according to group number
		my $tmp_string = $hsp->sb."-".$hsp->se.":".$hsp->qb."-".$hsp->qe."(".$hsp->percent.")"."; ";
		push @{$groups[($hsp->group - 1)]}, $tmp_string;
		$TOT_COV[$hsp->group - 1] += $hsp->length;  #NOTE: took out the -1 for array
		$TOT_MAT[$hsp->group - 1] += $hsp->match;   #NOTE: took out the -1 for array
	}

	#foreach my $array (@groups) {   #for each group number there are an array of coords
	for (my $i = 0; $i <= $#groups; $i++) {
		@{$groups[$i]} = sort {$a cmp $b} @{$groups[$i]};
		my $min = ${$groups[$i]}[0];
		my $max = ${$groups[$i]}[$#{$groups[$i]}];
		my $cov = sprintf("%.1f", (100*($TOT_COV[$i]/$sbjct->length)));
		my $id  = sprintf("%.1f", (100*($TOT_MAT[$i]/$TOT_COV[$i])));
		my $q_gi = $query ? $query : '0';
		$string .= "$gi\t$q_gi\t".$cov."\t".$id."\t@{$groups[$i]}\n";   #a line has to made for each line	
	}

	return $string;   #return, each line corresponds to a group
}


sub group_list_long { ### DELUXE CODE jg ###
	#purpose: to return a very informative string for a subject found in a BLAST query
	my $sbjct = shift;
	my $query = @_ ? shift : '0';  #this allows you to hand it the query object if you like
	my $gi;
	
#	print "QL = ", $query->length , "\n";
	if ($sbjct->name =~ /gi\|(\d+)/) { $gi = $1;}
	elsif ($sbjct->name =~ /^>\s*(\S+)/) { $gi = $1;}
	else { die "no identifier in defline $!";}
	
	my $string;
	my @groups = [];
	my @TOT_COV = ();
	my @TOT_MAT = ();
	my @TOT_SCO = ();
	
	
	while(my $hsp = $sbjct->nextHSP) {  # divy up according to group number
		my $tmp_string = $hsp->sb."-".$hsp->se.":".$hsp->qb."-".$hsp->qe."(".$hsp->percent.")"."; ";
		push @{$groups[($hsp->group - 1)]}, $tmp_string;
		$TOT_COV[$hsp->group - 1] += $hsp->length;
		$TOT_MAT[$hsp->group - 1] += $hsp->match;
		$TOT_SCO[$hsp->group - 1] += $hsp->bits;
	}

	#foreach my $array (@groups) {   #for each group number there are an array of coords
	for (my $i = 0; $i <= $#groups; $i++) {
		@{$groups[$i]} = sort {$a cmp $b} @{$groups[$i]};
		my $min = ${$groups[$i]}[0];
		my $max = ${$groups[$i]}[$#{$groups[$i]}];
		my $cov = sprintf("%.1f", (100*($TOT_COV[$i]/$sbjct->length)));
		my $qcov = sprintf("%.1f", (100*($TOT_COV[$i]/$sbjct->{PARENT}->length)));
		my $id  = sprintf("%.1f", (100*($TOT_MAT[$i]/$TOT_COV[$i])));
		my $score = $TOT_SCO[$i];
		my $q_gi = $query ? $query : '#';
		$string .= "$gi\t$q_gi\t".$cov."\t".$qcov."\t".$id."\t".$score."\t@{$groups[$i]}\n";   #a line has to made for each line	
	}

	return $string;   #return, each line corresponds to a group
}




###############################################################################
# BPdeluxe::HSP
###############################################################################
package BPdeluxe::HSP;
use overload '""' => '_overload';
sub new {
	my $hsp = bless {};
	($hsp->{SCORE}, $hsp->{BITS},
		$hsp->{MATCH}, $hsp->{POSITIVE}, $hsp->{LENGTH},$hsp->{P}, $hsp->{E}, $hsp->{STRAND},
		$hsp->{QB}, $hsp->{QE}, $hsp->{SB}, $hsp->{SE}, $hsp->{QL},
		$hsp->{SL}, $hsp->{AS}, $hsp->{QG}, $hsp->{SG}, $hsp->{GRP}, $hsp->{SP}) = @_;  	## Deluxe Code ##
	   $hsp->{PERCENT} = int(1000 * $hsp->{MATCH}/$hsp->{LENGTH})/10;
	return $hsp;
}
sub _overload {
	my $hsp = shift;
	return $hsp->queryBegin."..".$hsp->queryEnd." ".$hsp->bits;
}
sub score           {shift->{SCORE}}
sub bits            {shift->{BITS}}
sub percent         {shift->{PERCENT}}
sub match           {shift->{MATCH}}
sub positive        {shift->{POSITIVE}}
sub length          {shift->{LENGTH}}
sub P               {shift->{P}}
sub E					  {shift->{E}} ##miao's code## 
sub strand			  {shift->{STRAND}}; ##miao's code##
sub queryBegin      {shift->{QB}}
sub queryEnd        {shift->{QE}}
sub sbjctBegin      {shift->{SB}}
sub sbjctEnd        {shift->{SE}}
sub queryAlignment  {shift->{QL}}
sub sbjctAlignment  {shift->{SL}}
sub alignmentString {shift->{AS}}
sub queryGaps       {shift->{QG}}
sub sbjctGaps       {shift->{SG}}
sub qb              {shift->{QB}}
sub qe              {shift->{QE}}
sub sb              {shift->{SB}}
sub se              {shift->{SE}}
sub qa              {shift->{QL}}
sub sa              {shift->{SL}}
sub as              {shift->{AS}}
sub qg              {shift->{QG}}
sub sg              {shift->{SG}}
sub group			{shift->{GRP}}		## Deluxe Code ##
sub sumP			{shift->{SP}}		## Deluxe Code ##

###############################################################################
# BPdeluxe::Multi
###############################################################################
package BPdeluxe::Multi;
sub new {
	my ($class, $fh) = @_;
	if (ref $fh !~ /GLOB/)
		{die "BPdeluxe error: new expects a GLOB reference not $fh\n"}
	my $this = bless {};
	$this->{FH} = $fh;
	return $this;
}
sub nextReport {
	my ($this) = @_;
	$this->_fastForward or return 0;
	my $blast = new BPdeluxe($this->{FH});
	return $blast;
}
sub _fastForward {
	my ($this) = @_;
	my $FILEHANDLE = $this->{FH};
	while(<$FILEHANDLE>) {
		if ($_ =~ /^T?BLAST[NPX]/) {
			$this->{LASTLINE} = $_;
			return 1;
		}
	}
	return 0;
}



1;
__END__

=head1 NAME

BPdeluxe - Lightweight BLAST parser

=head1 SYNOPSIS

 use BPdeluxe;
 
 # single BLAST report
 
 my $report = new BPdeluxe(\*STDIN);
 $report->query;
 $report->database;
 while(my $sbjct = $report->nextSbjct) {
     $sbjct->name;
	 $sbjct->length;
     $sbjct->group_list   #jg
	 while (my $hsp = $sbjct->nextHSP) {
         $hsp->score;
         $hsp->bits;
         $hsp->percent;
         $hsp->P;        
         $hsp->match;
         $hsp->positive;
         $hsp->length;
         $hsp->queryBegin;
         $hsp->queryEnd;
         $hsp->sbjctBegin;
         $hsp->sbjctEnd;
         $hsp->queryAlignment;
         $hsp->sbjctAlignment;
         $hsp->alignmentString;
		 $hsp->queryGaps;
		 $hsp->sbjctGaps;
     	 $hsp->group   #ry
		 $hsp->sumP    #jg
	 }
 }
 
 # multiple (concatenated) BLAST reports
 
 my $multiple_report = new BPdeluxe::Multi(\*STDIN);
 while(my $blast = $multiple_report->nextReport) {
     while(my $sbjct = $blast->nextSbjct) {
         while(my $hsp = $sbjct->nextHSP) {
         }
     }
 }


=head1 DESCRIPTION

BPdeluxe is a package for parsing BLAST reports. The BLAST programs are a family
of widely used algorithms for sequence database searches. The reports are
non-trivial to parse, and there are differences in the formats of the various
flavors of BLAST. BPdeluxe parses BLASTN, BLASTP, BLASTX, TBLASTN, and TBLASTX
reports from both the high performance WU-BLAST, and the more generic
NCBI-BLAST.

Many people have developed BLAST parsers (I myself have made at least three).
BPdeluxe is for those people who would rather not have a giant object
specification, but rather a simple handle to a BLAST report that works well
in pipes.

=head2 Object

BPdeluxe has three kinds of objects, the report, the subject, and the HSP. To
create a new report, you pass a filehandle reference to the BPdeluxe constructor.

 my $report = new BPdeluxe(\*STDIN); # or any other filehandle

The report has two attributes (query and database), and one method (nextSbjct).

 $report->query;     # access to the query name
 $report->database;  # access to the database name
 $report->nextSbjct; # gets the next subject
 while(my $sbjct = $report->nextSbjct) {
     # canonical form of use is in a while loop
 }

A subject is a BLAST hit, which should not be confused with an HSP (below). A
BLAST hit may have several alignments associated with it. A useful way of
thinking about it is that a subject is a gene and HSPs are the exons. Subjects
have one attribute (name) and one method (nextHSP).

 $sbjct->name;    # access to the subject name
 "$sbjct";        # overloaded to return name
 $sbjct->nextHSP; # gets the next HSP from the sbjct
 while(my $hsp = $sbjct->nextHSP) {
     # canonical form is again a while loop
 }

An HSP is a high scoring pair, or simply an alignment. HSP objects do not have
any methods, just attributes (score, bits, percent, P, match, positive, length,
queryBegin, queryEnd, sbjctBegin, sbjctEnd, queryAlignment, sbjctAlignment)
that should be familiar to anyone who has seen a blast report. For
lazy/efficient coders, two-letter abbreviations are available for the
attributes with long names (qb, qe, sb, se, qa, sa).

 $hsp->score;
 $hsp->bits;
 $hsp->percent;
 $hsp->P;
 $hsp->match;
 $hsp->positive;
 $hsp->length;
 $hsp->queryBegin;      $hsp->qb;
 $hsp->queryEnd;        $hsp->qe;
 $hsp->sbjctBegin;      $hsp->sb;
 $hsp->sbjctEnd;        $hsp->se;
 $hsp->queryAlignment;  $hsp->qa;
 $hsp->sbjctAlignment;  $hsp->sa;
 $hsp->alignmentString; $hsp->as;
 $hsp->queryGaps;       $hsp->qg;
 $hsp->sbjctGaps;       $hsp->sg;
 "$hsp"; # overloaded for qb..qe bits

I've included a little bit of overloading for double quote variable
interpolation convenience. A report will return the query and database, a
subject will return its name and an HSP will return its queryBegin, queryEnd,
and bits in the alignment. Feel free to modify this to whatever is most
frequently used by you.

So a very simple look into a BLAST report might look like this.

 my $report = new BPdeluxe(\*STDIN);
 while(my $sbjct = $report->nextSbjct) {
     print "$sbjct\n";
     while(my $hsp = $sbjct->nextHSP) {
	 	print "\t$hsp\n";
     }
 }

The output of such code might look like this:

 >foo
     100..155 29.5
     268..300 20.1
 >bar
     100..153 28.5
     265..290 22.1


=head2 Concatenated BLAST reports

You can step through multiple BLAST reports if you have a file of concatenated
BLAST reports using the following construct.

 my $multiple_report = new BPdeluxe::Multi(\*STDIN);
 while(my $blast = $multiple_report->nextReport) {
     while(my $sbjct = $blast->nextSbjct) {
         while(my $hsp = $sbjct->nextHSP) {
         }
     }
 }

=head2 Additional Features of BPdeluxe (over BPlite):
 $hsp->group			returns the group #		## Deluxe Code RY##
 $hsp->sumP			returns # of HSPS		## Deluxe Code JG##
 $sbjct->group_list  returns synop. of hits  ## Deluxe Code JG##


 So a very simple look at the group_list method may be as follows.

 my $report = new BPdeluxe(\*STDIN);
 while(my $sbjct = $report->nextSbjct) {
     print $sbjct->group_list
 }

 The output of such code might look like this.
 NOTE: each line corresponds to a different HSP group.

 338422  0       100.0   99.4    3-44:1885532-1885573(100);  44-680:1886286-1886922(99.3);
 338422  0       100.0   99.4    3-44:1843226-1843185(100);  44-680:1842472-1841836(99.3);
 338422  0       93.5    92.8    44-676:1805375-1804742(92.7);
 338422  0       93.4    92.4    44-676:1819479-1818848(92.4); 



=head1 AUTHOR

Ian Korf (ikorf@sapiens.wustl.edu, http://sapiens.wustl.edu/~ikorf)

=head1 CONTRIBUTING AUTHORS

Raymond Yeh      (ryeh@sapiens.wustl.edu)
Jarret Glasscock (jglassco@sapiens.wustl.edu)

=head1 ACKNOWLEDGEMENTS

This software was developed at the Genome Sequencing Center at Washington
Univeristy, St. Louis, MO.

=head1 COPYRIGHT

Copyright (C) 1999 Ian Korf. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut










