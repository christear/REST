#!/usr/bin/perl -w
use strict;
=use
search the motif in the upstream of cleavage site ....
....
=cut

=test
my $test = "agctgggg";
$test =~ tr/[agct]/[AGCT]/;
print "$test\n";
=cut

my @MOTIFS=('AATAAA', 'ATTAAA', 'TATAAA', 'AGTAAA', 'AAGAAA', 'AATATA', 'AATACA', 'CATAAA', 'GATAAA','AATGAA', 'TTTAAA','ACTAAA','AATAGA');
#print "$MOTIFS[0]\t$MOTIFS[$#MOTIFS]\n";
my $start = 0;
my $end = 104;
sub searchMotif {
	my ($fa,$output,$start,$end) = @_;
	open IN,"$fa" or die $!;
	my $id;
	my $strand;
	my $seq;
	my $num = 0;
	open OUT,">$output" or die $!;
	my %Sequence;
	while(<IN>){
		chomp;
		if(/>/){
			$id = substr($_,1);
			$num++;
			$seq = "";
			#print OUT "$id\t";
		}else{
			my $subseq = $_;
			$subseq =~ tr/[agct]/[AGCT]/;
			$seq = $seq.$subseq;
			$Sequence{$id} = $seq;
		}
	}
	close IN;
	foreach my $id(keys %Sequence){
		print OUT "$id\t";
		my $seq = $Sequence{$id};
		my $sub = substr($seq,$start,$end-$start);
		print "$id\t$seq\t$sub\n";
		my @Pos;
		foreach my $m(@MOTIFS){
			my $p = -1;
			if($sub =~ /$m/){
				$p = index($sub,$m);
			}
			## convert to 1-based position 
			$p++;
			push @Pos,$p;
		}
		my $posl = join("\t",@Pos);
		print OUT "$posl\n";
	}
	close OUT;
} 

if(@ARGV < 2){
	print "usage: fasta_input\toutput_tab\n";
	exit(0);
}else{
	if(@ARGV > 2){
		$start = $ARGV[2];
	}
	if(@ARGV > 3){
		$end = $ARGV[3];
	}
	print "search sequence from $start to $end\n";
	&searchMotif($ARGV[0],$ARGV[1],$start,$end);
}
