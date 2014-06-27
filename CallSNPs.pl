#!/usr/bin/perl
use warnings;
use strict;
use threads;
use Thread::Queue;
use FindBin;
use lib "$FindBin::Bin/Lib";
use Configuration;
use Tools;

my $configFile=$ARGV[0];

die "usage : perl $0 <config file governing all alignment>\n\n" unless $#ARGV==0;

my $q = Thread::Queue->new();
my $config = Configuration->new($configFile);
my $threads = $config->get("OPTIONS","ManagerThreads");
my @Groups = $config->getAll("GROUPS");

for(my $i=0;$i<=$#Groups;$i++){
	warn "enqueuing $i ($Groups[$i])\n";
	$q->enqueue($Groups[$i]);
}

for(my$i=0;$i<=$threads;$i++){
	my $thr=threads->create(\&workerThread);
}
while(threads->list()>0){
	my @thr=threads->list();
	$thr[0]->join();
}


sub workerThread{
	while(my $work=$q->dequeue_nb()){
		my $grp		= $work;
		my $DataDir 	= $config->get("DIRECTORIES","Data");
		my $OutDir  	= $config->get("DIRECTORIES","Output");
		my $RefDir	= $config->get("DIRECTORIES","References");
		my $workThreads = $config->get("OPTIONS","BWAThreads");
		my $bwa		= $config->get("PATHS","bwa");
		my $samtools	= $config->get("PATHS","samtools");

		my @CurrentSourcePaths;
		my @GarbageCollector;

		my $file1=$DataDir."/".$grp.".R1.fastq";
		my $file2=$DataDir."/".$grp.".R2.fastq";
	
		die "Cannot find read 1 for group: $grp\nFile missing: $file1\nexiting...\n" unless -e $file1;
		die "Cannot find read 2 for group: $grp\nFile missing: $file2\nexiting...\n" unless -e $file2;
		
		my @Indicies 	= split(",",$config->get("GROUPS",$grp));
		foreach my $index (@Indicies){
			my $IndexPath=$RefDir."/".$index;
			my $alias=$index;
			my $baseOutput = $OutDir."/".$grp."_vs_".$alias;
			$alias =~ s/\..+//;
			my $command = "$bwa mem -t $workThreads $IndexPath $file1 $file2 > $baseOutput.sam";
			warn $command."\n";
			`$command`;
			$command = "$samtools view -bS $baseOutput.sam > $baseOutput.bam";
			warn $command."\n";
			`$command`;
			$command = "$samtools sort -\@ $workThreads $baseOutput.bam $baseOutput.sorted";
			warn $command."\n";
			`$command`;
			$command = "$samtools index $baseOutput.sorted.bam";
			warn $command."\n";
			`$command`;
			push @GarbageCollector, $baseOutput.".sam";
		}
		collectTheGarbage(@GarbageCollector);
	}
}

sub collectTheGarbage {
	my @files = @_;
	foreach my $file (@files){
		my $command="rm -rf $file";
		warn $command."\n";
		`$command`;
	}
	return 1;
}

sub prepFinal {
	my $finalDir = shift @_;
	my @files = @_;
	foreach my $file (@files){
		my $sPath=$file;
		my $oPath=$file;
		$oPath=~s/.+\///g;
		$oPath=$finalDir."/".$oPath;
		my $command = "mv $sPath $oPath";
		warn $command."\n";
		`$command`;
	}
	return 1;
}

