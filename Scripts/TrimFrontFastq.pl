#!/usr/bin/perl
use warnings;
use strict;

my $fastq=$ARGV[0];
my $trim=$ARGV[2];
my $output=$ARGV[1];
die "usage: perl $0 <fastq input file (MUST end in .fastq)> <trim # from 5' end>\n\n" unless $#ARGV==2;
open(OUT,">",$output) || die "cannot open $output!\n$!\nexiting....\n";
open(FASTQ,"<",$fastq) || die "cannot open $fastq!\n$!\nexiting...\n";
until(eof(FASTQ)){
	my $h1=<FASTQ>;
	my $s1=<FASTQ>;
	my $h2=<FASTQ>;
	my $s2=<FASTQ>;
	chomp $h1;
	chomp $s1;
	chomp $h2;
	chomp $s2;
	my $junk;
	$junk=substr($s1,0,$trim,"");
	$junk=substr($s2,0,$trim,"");
	print OUT $h1."\n".$s1."\n".$h2."\n".$s2."\n";
}
close FASTQ;
close OUT;
