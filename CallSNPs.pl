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
my $mchScriptDir = $FindBin::Bin."/Lib/SysCall_1_1/";
my $mchScript=$mchScriptDir."/SysCall.pl";

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
		my $bcftools	= $config->get("PATHS","bcftools");
		my $snpRate	= $config->get("OPTIONS","snpRate");
		my $minCov	= $config->get("OPTIONS","minCov");
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
			$command = "$samtools mpileup -F 0.00001 -g -C50 -d 10000000 -f $IndexPath $baseOutput.sorted.bam | $bcftools view -b -m 0.01 -p .99 - | $bcftools view - > $baseOutput.raw.vcf";
			warn $command."\n";
			`$command`;
			my %H = %{parseResults("$baseOutput.raw.vcf",$baseOutput.".filt.vcf",$snpRate,$minCov)};
			generateMeacham($IndexPath,$baseOutput.".mch.tab",\%H);
			$command = "perl $mchScript $baseOutput.mch.tab $baseOutput.sam $baseOutput $mchScriptDir";
			warn $command."\n";
			`$command`;
			push @GarbageCollector, $baseOutput.".mch.tab";
			push @GarbageCollector, $baseOutput.".reads_parsed";
			push @GarbageCollector, $baseOutput.".sys_errors.nn";
			push @GarbageCollector, $baseOutput.".heterozygous.nn";
			push @GarbageCollector, $baseOutput.".table_all";
			push @GarbageCollector, $baseOutput.".table_chosen";
			push @GarbageCollector, $baseOutput.".table_chosen.n";
			push @GarbageCollector, $baseOutput.".table_chosen.nn";
			push @GarbageCollector, $baseOutput.".filt.vcf";
			push @GarbageCollector, $baseOutput.".fasta.raw.vcf";
			my %bad=%{getSysErrors($baseOutput)};
			my $outFinal = $baseOutput.".final.vcf";
			my @outFinal;
			foreach my $key (keys %H){
				my $K=$H{$key}{"Chr"}."-".$H{$key}{"Pos"};
				if(defined($bad{$K})){
				}else{
	#				push @outFinal, parseToFinal($H{$key}{"L"});
					push @outFinal, $H{$key}{"L"};
				}
			}
			Tools->printToFile($outFinal,\@outFinal);
		}
		collectTheGarbage(@GarbageCollector);
	}
}

sub parseToFinal {
	my $line=shift;
	my @line=split(/\t/,$line);
	my $chr=$line[0];
	my $pos=$line[1];
	my $refBase=$line[3];
	next if $line[4] eq "X";
        my @posAlt =split(/\,/,$line[4]);
	my %I=%{parseInfo($line[7])};
			
}

sub getSysErrors {
	my $file=shift;
	my @file=@{Tools->LoadFile($file.".sys_errors")};
	my $head=shift @file;
	my %H;
	foreach my $line (@file){
		$line=~s/\s.+//;
		my ($chr,$coord)=split(/\:/,$line);
		my $k="$chr-$coord";
		$H{$k}=1;
	}
	return \%H;
}

sub generateMeacham {
	my $file=shift;
	my $out=shift;
	my %H=%{$_[0]};
	my @output;
	my %Fasta=%{Tools->LoadFasta($file)};
	foreach my $pos (keys %H){
		my $chr=$Fasta{$H{$pos}{"Chr"}};
		next if $pos <3;
		my $rpos=$pos-3;  ### -1 for coordinate correct, -2 for start of meacham.
		my $seq=substr($chr,$rpos,5);
		my @seq=split(//,$seq);
		push @output, $H{$pos}{"Chr"}." ".$pos." ".join(" ",@seq);
	}
	Tools->printToFile($out,\@output);
	return 1;
}

sub parseResults {
	my $file=shift;
	my $output=shift;
	my $rate=shift;
	my $minCov=shift;
	my @file=@{Tools->LoadFile($file)};
	my @out;
	my %R;
	foreach my $line (@file){
		my @line=split(/\t/,$line);
		my $chr=$line[0];
		my $pos=$line[1];
		my $refBase=$line[3];
		next if $line=~m/^\#/;
		next if $line[4] eq "X";
		next if $line =~m/INDEL/;
		my @posAlt =split(/\,/,$line[4]);
		my %I=%{parseInfo($line[7])};
		unless(defined $I{"DP"}){
			warn $line."\n";
			die "malformed output from SNP calling!\n";
		}
		unless(defined $I{"I16"}){
			warn $line."\n";
			die "malformed output from SNP calling!\n"
		}
		unless( defined $I{"QS"}){
			warn $line."\n";
			die "malformed output from SNP calling!\n";
		}
		my %r;
		$r{"Chr"}=$chr;
		$r{"Pos"}=$pos;
		$r{"DP"}=$I{"DP"};
		$r{"I16"}=$I{"I16"};
		$r{"QS"}=$I{"QS"};
		$r{"L"}=$line;
		if(checkInfo(\%r,$rate,$minCov)){
			$R{$pos}=\%r;
			push @out, $line;
		}
	}
#	Tools->printToFile($output,\@out);
	return \%R;
}

sub checkInfo {
	my %r=%{$_[0]};
	my $rate=$_[1];
	my $minCov=$_[2];
	my @I=split(/\,/,$r{"I16"});
	if($r{"DP"}<$minCov){
		return 0;
	}
	my $s=$I[2]+$I[3];
	if (($s/$r{"DP"})>$rate){
		return 1;
	}
	return 0;
}

sub parseInfo {
	my $info=shift;
	my @info=split(/\;/,$info);
	my %h;
	foreach my $stanza (@info){
		my ($key,$value)=split(/\=/,$stanza);
		$h{$key}=$value;
	}
	return \%h;
}

sub collectTheGarbage {
	my @files = @_;
	foreach my $file (@files){
		my $command="rm -rf $file";
		next unless -e $file;
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

