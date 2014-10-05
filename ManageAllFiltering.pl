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

die "usage : perl $0 <config file governing all filtering>\n\n" unless $#ARGV==0;

my $q = Thread::Queue->new();
my $config = Configuration->new($configFile);
my $threads = $config->get("OPTIONS","BWAThreads");
my @Groups = $config->getAll("FILTERING");

for(my $i=0;$i<=$#Groups;$i++){
	warn "enqueuing $i ($Groups[$i])\n";
	$q->enqueue($Groups[$i]);
}

for(my$i=0;$i<1;$i++){
	my $thr=threads->create(\&workerThread);
}
while(threads->list()>0){
	my @thr=threads->list();
	$thr[0]->join();
}


sub workerThread{
	while(my $work=$q->dequeue_nb()){
		my $grp=$work;
		my $DataDir = $config->get("DIRECTORIES","Data");
		my $OutDir  = $config->get("DIRECTORIES","Filtered");
		my @SubDirs = split(",",$config->get("FILTERING",$grp));
		my $workThreads = $config->get("OPTIONS","workThreads");
		my @CurrentSourcePaths;
		my @GarbageCollector;
		foreach my $sub (@SubDirs){
			my $dir=$DataDir."/".$sub;
			warn "Working with $dir...\n";
			opendir(DIR,$dir) || die "cannot open directory: $dir!\n$!\nexiting...\n";
			my @CurrentFiles=grep {m/fastq/} readdir DIR;
			closedir DIR;
			map {push @CurrentSourcePaths, $dir."/".$_} @CurrentFiles;
		}
		if($config->get("PIPELINE","Compressed")){
			warn "Treating files as compressed.\n";
			my @thisSet = @CurrentSourcePaths;
			@CurrentSourcePaths=();
			foreach my $file (@thisSet){
				die "file: $file is not a gzipped file.\n" unless $file=~m/\.gz$/;
				my $nPath=$file;
				$nPath=~s/\.gz//;
				my $command = "gunzip ".$file;
				warn $command."\n";
				`$command`;
				push @CurrentSourcePaths, $nPath;
				#push @GarbageCollector, $nPath;
			}	
		}
		if($config->get("PIPELINE","FivePrimeFilter")){
			warn "Trimming 5' end of sequences...\n";
			my @thisSet = @CurrentSourcePaths;
			@CurrentSourcePaths=();
			my $script=$config->get("PATHS","FivePrimeTrimmer");
			my $length=$config->get("OPTIONS","LengthOf5pTrim");
			foreach my $file (@thisSet){
				my $oPath=$file;
				$oPath=~s/fastq/5pTrim\.fastq/;
				$oPath=~s/fq/5pTrim\.fastq/;
				my $command="perl $script $file $oPath $length";
				warn $command."\n";
				`$command`;
				push @CurrentSourcePaths, $oPath;
				push @GarbageCollector, $oPath if $config->get("PIPELINE","Compressed");
			}
		}
		if($config->get("PIPELINE","ThreePrimeFilter")){
			warn "Trimming 3' end on quality...\n";
			my @thisSet = @CurrentSourcePaths;
			@CurrentSourcePaths=();
			my $script = $config->get("PATHS","fastq_quality_trimmer");
			my $MinL   = $config->get("OPTIONS","Min3pLength");
			my $MinQ   = $config->get("OPTIONS","Min3pQuality");
			my $opts   = "-t $MinQ -Q 33 -l $MinL";
			foreach my $file (@thisSet){
				my $oPath=$file;
				$oPath=~s/fastq/3pTrim\.fastq/;
				$oPath=~s/fq/3pTrim\.fastq/;
				my $command="$script $opts -i $file -o $oPath";
				warn $command."\n";
				`$command`;
				push @CurrentSourcePaths, $oPath;
				push @GarbageCollector, $oPath;
			}
		}
		if($config->get("PIPELINE","Paired")){
			warn "Parsing for Pairs...\n";
			my @thisSet = @CurrentSourcePaths;
			@CurrentSourcePaths=();
			my $script = $config->get("PATHS","PairsAndOrphans");
			my @R1=grep {m/\_R1\_/} @thisSet;
			my @R2=grep {m/\_R2\_/} @thisSet;
			my $T1=$DataDir."/$grp.TempRead1.fastq";
			my $T2=$DataDir."/$grp.TempRead2.fastq";
			my $command="cat ".join(" ",@R1)." > $T1";
			warn $command."\n";
			`$command`;
			$command   ="cat ".join(" ",@R2)." > $T2";
			warn $command."\n";
			`$command`;
			push @GarbageCollector, $T1;
			push @GarbageCollector, $T2;
			my $O=$DataDir."/$grp";
			my $OR1=$DataDir."/$grp".".R1.fastq";
			my $OR2=$DataDir."/$grp".".R2.fastq";
			my $ORO=$DataDir."/$grp".".orphan.fastq";
			$command="perl $script $T1 $T2 $O";
			warn $command."\n";
			`$command`;
			push @GarbageCollector, $OR1;
			push @GarbageCollector, $OR2;
			push @GarbageCollector, $ORO;
			push @CurrentSourcePaths, $OR1;
			push @CurrentSourcePaths, $OR2;
			push @CurrentSourcePaths, $ORO;
		}
		prepFinal($OutDir,@CurrentSourcePaths);
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

