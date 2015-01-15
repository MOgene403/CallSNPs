#!/usr/bin/perl

use warnings;
use strict;

use lib '/home/hpriest/Scripts/Library';

my $usage = "
perl $0 <fastq file 1> <fastq file 2> <output root>\n\n";

die $usage unless $#ARGV==2;

my $FQ1=$ARGV[0];
my $FQ2=$ARGV[1];
my $OR =$ARGV[2];
my $outOne=$OR.".R1.fastq";
my $outTwo=$OR.".R2.fastq";
my $outOrphan=$OR.".orphan.fastq";


open(FQ1,"<",$FQ1) || die "cannot open $FQ1!\n$!\nexiting...\n";
open(FQ2,"<",$FQ2) || die "cannot open $FQ2!\n$!\nexiting...\n";
open(O1,">",$outOne) || die "cannot open $outOne!\n$!\nexiting...\n";
open(O2,">",$outTwo) || die "cannot open $outTwo!\n$!\nexiting...\n";
open(OO,">",$outOrphan) || die "cannot open $outOrphan!\n$!\nexiting...\n";
my %D;
my $i=0;
while(1){
	unless(eof(FQ1)){
		my $H1s=<FQ1>;
		my $S1 =<FQ1>;
		my $H1q=<FQ1>;
		my $Q1 =<FQ1>;
		chomp $H1s;
		chomp $S1;
		chomp $H1q;
		chomp $Q1;
		my ($R,$NR)=split(/\s/,$H1s);
#		$NR=~m/\d\:[A-Z]\:\d\:([ACGT]\{4,12\})$/;
		if(defined($D{$R})){
			print O1 $H1s."\n".$S1."\n".$H1q."\n".$Q1."\n";
			print O2 $D{$R};
			delete($D{$R});
		}else{
			$D{$R}=$H1s."\n".$S1."\n".$H1q."\n".$Q1."\n";
		}
	}
	unless(eof(FQ2)){
		my $H1s=<FQ2>;
		my $S1 =<FQ2>;
		my $H1q=<FQ2>;
		my $Q1 =<FQ2>;
		chomp $H1s;
		chomp $S1;
		chomp $H1q;
		chomp $Q1;
		my ($R,$NR)=split(/\s/,$H1s);
		if(defined($D{$R})){
			print O2 $H1s."\n".$S1."\n".$H1q."\n".$Q1."\n";
			print O1 $D{$R};
			delete($D{$R});
		}else{
			$D{$R}=$H1s."\n".$S1."\n".$H1q."\n".$Q1."\n";;
		}	
	}
	last if ((eof(FQ1))&&(eof(FQ2)));
	$i++;
}

foreach my $entry (keys %D){
	print OO $D{$entry};
}

close O1;
close O2;
close OO;
close FQ1;
close FQ2;

