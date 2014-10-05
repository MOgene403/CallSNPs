#!/usr/bin/perl
use warnings;
use strict;
use FindBin;
use lib "$FindBin::Bin/Lib";

use Tools;
use Configuration;

use threads;
use threads::shared;
use Thread::Queue;

my $q = Thread::Queue->new();
die "usage: perl $0 <Config file>\n\n" unless $#ARGV==0;
my $Config = Configuration->new($ARGV[0]);

my $nThreads = $Config->get("OPTIONS","BWAThreads");

warn "Recognizing $nThreads as max threading...\n";

my $ref=$Config->get("PATHS","References");

warn "Finding Vectors...\n";
my @LineNo = $Config->getAll("GROUPS");

foreach my $i (@LineNo){
      $q->enqueue($i);
}
for(my$i=0;$i<1;$i++){
      my $thr=threads->create(\&worker);
}
while(threads->list()>0){
      my @thr=threads->list();
      $thr[0]->join();
}


sub worker {
	my $TID=threads->tid() -1 ;
	while(my$j=$q->dequeue_nb()){
		my @Refs=split(/\,/,$Config->get("GROUPS",$j));
		my $P1=$Config->get("DIRECTORIES","Filtered")."/".$j.".R1.fastq";
		my $P2=$Config->get("DIRECTORIES","Filtered")."/".$j.".R2.fastq";
		
		my $outputDir = $Config->get("DIRECTORIES","Output")."/".$j;
		mkdir $outputDir unless -e $outputDir;
		my $base = $j;
		my $samtools = $Config->get("PATHS","samtools");
		foreach my $ref (@Refs){
			my $bwaRef=$Config->get("DIRECTORIES","References")."/".$ref;

			my $bwaRoot=$outputDir."/$base.vs.$ref.Alignments";
		
			my $bwaAln=$bwaRoot.".bam";

			my $cmd=$Config->get("PATHS","bwa")." mem -t $nThreads $bwaRef $P1 $P2 | $samtools view -bS - > $bwaAln";
			warn $cmd."\n";
			`$cmd`;
			my $sorted=$bwaRoot.".sorted";
			$cmd = $samtools." sort $bwaAln $sorted";
			`$cmd`;
			$cmd = $samtools." index ".$sorted.".bam";
			`$cmd`;
		}
	}
}

# /home/ec2-user/Store1/bin/delly  -t TRA -o TRA.vcf -q 20 -g TwoChrom.fasta pGC1_Raw.sorted.bam

exit(0);


