#!/usr/bin/perl
#
# This script is used to compare 2 CNV profiles and compute a distance metric
#
# Files describing the CNV profiles should be TAB-delimited
# and include the following information:
# chromosome, start-end, LRR, copy-number (optionnal)
#
# example:
#
# chr2	90397715-91794601	-0.704403	1
# chr5	71904-46115086	0.377794	3
# chr5	69238677-70587018	-0.450079	1
# chr6	259318-293615	-0.962523	1
# chr9	71129824-141018984	0.328202	3
#
# Author: Stephane Wenric s.wenric@ulg.ac.be
#

use strict;
use warnings;
use v5.10;

use Math::Round;

my %CNVProfile1;
my %CNVProfile2;

my $precision = 40000;	

# feel free to adjust this value. Default here is 40 kb.
# Keep in mind that this is a somewhat naive implementation, and
# if you want to use a lower precision, this algorithm 
# should be adapted to use something like Set::IntSpan
# instead of iterating on each chromosomal position.
#
# Several tests have shown no change in scores when adjusting
# the precision at values lower than 50 kb.

my $usage = "Usage: perl $0 <file 1> <file 2>\n";

# file format:
# chr <TAB> start-end <TAB> LRR <TAB> copy-number (optionnal)
# eg. chr1	16335799-16379884	0.348422	3

my $file1 = shift || die($usage);
my $file2 = shift || die($usage);

my $scoreLRR = 0;
my $scoreCN = 0;
my $nbSegmentsWAlterations = 0;
my $nbSegmentsTotal = 0;

# chromosome lengths values for hg19
# adjust accordingly for other organisms
my %chrLengths = (
	"chr1"	=>	249250621,
	"chr2"	=>	243199373,
	"chr3"	=>	198022430,
	"chr4"	=>	191154276,
	"chr5"	=>	180915260,
	"chr6"	=>	171115067,
	"chr7"	=>	159138663,
	"chrX"	=>	155270560,
	"chr8"	=>	146364022,
	"chr9"	=>	141213431,
	"chr10"	=>	135534747,
	"chr11"	=>	135006516,
	"chr12"	=>	133851895,
	"chr13"	=>	115169878,
	"chr14"	=>	107349540,
	"chr15"	=>	102531392,
	"chr16"	=>	90354753,
	"chr17"	=>	81195210,
	"chr18"	=>	78077248,
	"chr20"	=>	63025520,
	"chrY"	=>	59373566,
	"chr19"	=>	59128983,
	"chr22"	=>	51304566,
	"chr21"	=>	48129895,	
);

# set the default LRR to 0 for all chromosomes for both profiles
for my $chr (sort sortChr keys %chrLengths) {
	my $length = round($chrLengths{$chr} / $precision);
	for (my $pos = 0; $pos <= $length; $pos++) {
		$CNVProfile1{$chr}{$pos} = 0;
		$CNVProfile2{$chr}{$pos} = 0;
		$nbSegmentsTotal++;
	}
}

# load LRRs for profile 1
open(my $in, "<", $file1) || die($usage);
while (my $line = <$in>) {
	chomp($line);
	my ($chr, $startEnd, $lrr, $copyNumber) = split(/\t/, $line);
	my ($start, $end) = split(/-/, $startEnd);
	if ($end - $start < $precision) { next; }
	$start = round($start/$precision);
	$end = round($end/$precision);
	for (my $pos = $start; $pos <= $end; $pos++) {
		$CNVProfile1{$chr}{$pos} = $lrr;
		$nbSegmentsWAlterations++;
	}
}
close($in);

# load LRRs for profile 2
open($in, "<", $file2) || die($usage);
while (my $line = <$in>) {
	chomp($line);
	my ($chr, $startEnd, $lrr, $copyNumber) = split(/\t/, $line);
	my ($start, $end) = split(/-/, $startEnd);
	if ($end - $start < $precision) { next; }
	$start = round($start/$precision);
	$end = round($end/$precision);
	for (my $pos = $start; $pos <= $end; $pos++) {
		$CNVProfile2{$chr}{$pos} = $lrr;
		$nbSegmentsWAlterations++;
	}
}
close($in);

for my $chr (sort sortChr keys %CNVProfile1) {
	for my $pos (sort keys %{$CNVProfile1{$chr}}) {
			my $lrr1 = $CNVProfile1{$chr}{$pos};
			my $lrr2 = $CNVProfile2{$chr}{$pos};

			my $cn1 = 2*(2**$CNVProfile1{$chr}{$pos});
			my $cn2 = 2*(2**$CNVProfile2{$chr}{$pos});

			$scoreLRR += abs($lrr2-$lrr1);
			$scoreCN += abs($cn2-$cn1);
	}
}

say "Score: ".(round(1000*$scoreCN/$nbSegmentsTotal)/100)." * 10^-1";


sub sortChr {
	my $chr1;
	my $chr2;
	if ($a =~ m/chr/i) { $chr1 = substr($a, 3); } else { $chr1 = $a; }
	if ($b =~ m/chr/i) { $chr2 = substr($b, 3); } else { $chr2 = $b; }
	if ($chr1 =~ /\d+/ && $chr2 =~ /\D+/) { return -1; }
	if ($chr1 =~ /\D+/ && $chr2 =~ /\d+/) { return 1; }
	if ($chr1 =~ /\d+/ && $chr2 =~ /\d+/) { return $chr1 <=> $chr2; }
	if ($chr1 =~ /\D+/ && $chr2 =~ /\D+/) { return $chr1 cmp $chr2; }
}

