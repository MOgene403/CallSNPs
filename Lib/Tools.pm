#!/usr/bin/perl

package Tools;

use warnings;
use strict;
use FindBin;
use lib "$FindBin::Bin/../Lib";
use Statistics::RankCorrelation;
use threads;
use Benchmark ':hireswallclock';

our $Bench1=undef;
our $Bench2=undef;

sub BenchStart {
	my $self=shift;
	$Bench1=Benchmark->new();
	return 1;
}

sub BenchStop {
	my $self=shift;
	die "Cannot call BenchStop() without first calling BenchStart()!\n" unless defined $Bench1;
	$Bench2=Benchmark->new();
	my $td=timediff($Bench2,$Bench1);
	warn "Total work (real) time: ". timestr($td) ."\n";
	return 1;
}


sub cleanProbesetName {
	my $self=shift;
	my $probeset=shift;
	$probeset=~s/F_at//;
	$probeset=~s/R_st//;
	return $probeset;
}

sub log10 {
	my $self=shift;
	my $value=shift;
	if($value==0){
		return 0;
	}else{
		return log($value)/log(10);
	}
}

sub log2 {
	my $self=shift;
	my $value=shift;
	if($value==0){
		return 0;
	}else{
		return log($value)/log(2);
	}
}

sub averageNormalizeArray { 
	my $self=shift;
	my @A=@{$_[0]};
	my $mean=mean($self,@A);
	for(my$i=0;$i<=$#A;$i++){
		$A[$i]=$A[$i]/$mean;
	}
	return \@A;
}

sub averageCentralizeArray { 
	my $self=shift;
	my @A=@{$_[0]};
	my $mean=mean($self,@A);
	for(my$i=0;$i<=$#A;$i++){
		$A[$i]=$A[$i]-$mean;
	}
	return \@A;
}

sub binValues {
	my $self=shift;
	my $d=shift;
	my @values=@_;
	my %bins;
	foreach my $value (@values){
		my $this=sprintf("%.$d"."f",$value);
		if(defined($bins{$this})){
			$bins{$this}++;
		}else{
			$bins{$this}=1;
		}
	}
	return \%bins;
}

sub GetLongestProtein {
	my $self=shift;
	my @prots=@_;
	my $M=0;
	my $N=-1;
	for(my$i=0;$i<=$#prots;$i++){
		if(length($prots[$i])>$M){
			$M=length($prots[$i]);
			$N=$i;
		}
	}
	return $prots[$N];
}

sub meanLog2 {
	my $self=shift;
	my @values=@_;
	my $sum=0;
	my $N=scalar(@values);
	for(my$i=0;$i<=$#values;$i++){
		my $real=2**$values[$i];
		$sum+=$real;
	}
	my $Average=$sum/$N;
	return log($Average)/log(2);
}

sub max {
	my $self=shift;
	my $max=0;
	map {$max=$_ if $_>$max} @_;
	return $max;
}

sub Max {
	my $self=shift;
	my @maximized=sort {$b <=> $a} @_;
	return $maximized[0];
}

sub median {
	my $self=shift;
	my @Array=sort {$a <=> $b} @_;
	my $index=(int(scalar(@Array)/2));
	return $Array[$index];
}

sub abs_max {
	my $self=shift;
	my $max=0;
	map {$max=abs($_) if ((abs($_))>$max)} @_;
	return $max;
}

sub Min {
	my $self=shift; #### not drunk.
	my @minimized=sort {$a <=> $b} @_;
	return $minimized[0];
}

sub min {
	my $self=shift; #### were you drunk or something when you wrote this? This is silly.
	my $min=Min($self,@_);   ### workaround for stupidity.
#	my $min=100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
#	map {$min=$_ if $_<$min} @_;
	return $min;
}


sub sum {
	my $self=shift;
	my $sum=0;
	map {$sum+=$_} @_;
	return $sum;
}

sub mean {
	my $self=shift;
      my $sum=0;
      map {$sum+=$_} @_;
	return 1 if scalar(@_)==0;
      return $sum/scalar(@_);
}

sub Variance {
	#### computes the Sample Variance of a set of numbers (in that Var(X)=covariance(X,X))
	#### Or, Var(X)=E(X^2)-(E(X))^2;
	my $self=shift;
	my @array=@_;
	my $mean=mean($self,@array);
	my @sqd;
	for(my$i=0;$i<=$#array;$i++){
		$sqd[$i]=$array[$i]*$array[$i];
	}
	my $VarX=mean($self,@sqd)-($mean**2);
	return $VarX;
}

sub stdev {
	my $self=shift;
      my @array=@_;
      my $mean = mean($self,@array);
	my @diffs;
	foreach my $val (@array){
		my $this=($val-$mean)**2;
		push @diffs, $this;
	}
	my $ad=mean($self,@diffs);
	return sqrt($ad);
}

sub LoadDir {
	my $self=shift;
	my $dir=shift;
	opendir(DIR,$dir) || die "cannoe open $dir!\n$!\nexiting...\n";
	my @files=readdir(DIR);
	closedir DIR;
	return \@files;
}

sub printWig {
	my $self=shift;
	my $out=shift;
	my %wiggle=%{$_[0]};
	open(OUT,">",$out) || die "cannot open $out!\n$!\nexiting...\n";
	foreach my $chr (keys %wiggle){
		print OUT "fixedStep\tchrom=$chr\tstart=1\tstep=1\n";
		for(my$i=0;$i<=$#{$wiggle{$chr}};$i++){
			print OUT ${$wiggle{$chr}}[$i]."\n";
		}
	}
	close OUT;
	return 1;
}

sub LoadWig {
	my $self=shift;
	my $file=shift;
	my %wiggle;
	my $chr;
	my $index;
	my @chrome;
	open(WIG,"<",$file) || die "Cannot open $file!\n$!\nexiting...\n";
	until(eof(WIG)){
		my $line=<WIG>;
		chomp $line;
#		fixedStep       chrom=scaffold_23       start=1 step=1
		if($line=~m/fixedStep\s+chrom\=(\w+\d+)\s/i){
			$chr=$1;
			$wiggle{$chr}=[];
		}elsif($line=~m/fixedStep\schrom\=(\w+)\s/){
			$chr=$1;
			$wiggle{$chr}=[];
		}else{
		die "$line\n" unless defined $chr;
			push @{$wiggle{$chr}}, $line;	
		}
	}
	close WIG;
	return \%wiggle;
}

sub LoadSingleWigChrome {
	my $self=shift;
	my $file=shift;
	my $chrome=shift;
	my %wiggle;
	my $chr;
	my $index;
	my @chrome;
	open(WIG,"<",$file) || die "Cannot open $file!\n$!\nexiting...\n";
	until(eof(WIG)){
		my $line=<WIG>;
		chomp $line;
#		fixedStep       chrom=scaffold_23       start=1 step=1
		if($line=~m/fixedStep\tchrom\=(\w+=*\d+)\t/i){
			$chr=$1;
			next unless $chr eq $chrome;
			$wiggle{$chr}=[];
		}else{
			next unless $chr eq $chrome;
			push @{$wiggle{$chr}}, $line;		
		}
	}
	close WIG;
	return \%wiggle;

}

sub findLarger {
	my $self=shift;
	my $v1=shift;
	my $v2=shift;
	return $v1 if $v1>$v2;
	return $v2 if $v2>$v1;
	return $v1 if $v1==$v2;
	return undef;
}

sub findSmaller {
	my $self=shift;
	my $v1=shift;
	my $v2=shift;
	return $v1 if $v1<$v2;
	return $v2 if $v2<$v1;
	return $v1 if $v1==$v2;
	return undef;
}

sub launch {
	my $cmd=shift;
	`$cmd`;
	return 1;
}

sub getZ {
	my $self=shift;
	my $obs=shift;
	my @arr=@{$_[0]};
	my $stdev=stdev($self,@arr);
	my $mean=mean($self,@arr);
	my $Z=($obs-$mean)/$stdev;
	return $Z;

}

sub runCmd {
	my $self=shift;
	my $cmd=shift;
	my $thr=threads->create(\&launch, $cmd);
	$thr->join();
	return 1;
}

sub printToFile {
	my $self=shift;
	my $file=shift;
	my @array=@{$_[0]};
	open(OUT,">",$file) || die "cannot open $file!\n$!\nexiting...\n";
	print OUT join("\n",@array)."\n";
	close OUT;
	return 1;
}

sub LoadFile {
	my $self=shift;
	my $file=shift;
	my @file;
	open(FILE,"<",$file) || die "Cannot open $file!\n$!\nexiting...\n";
	while(<FILE>){
		chomp $_;
		push @file, $_;
	}
	close FILE;
	return \@file;
}

sub getOverlap {
	my $self=shift;
	my @a1=@{$_[0]};
	my @a2=@{$_[1]};
	my (%a1,%a2);
	map {$a1{$_}=1} @a1;
	map {$a2{$_}=1} @a2;
	my @over;
	foreach my $key (keys %a1){
		push @over, $key if(defined($a2{$key}));
	}
	return \@over;
}

sub getRandomIndicies {
	my $self=shift;
	my $max=shift;
	my $num=shift;
	my %chosen;
	my $selected=0;
	until($selected==$num){
		my $index=int(rand($max));
		next if defined($chosen{$index});
		$chosen{$index}=1;
		$selected++;
	}
	return \%chosen;
}

sub LoadFastq {
	my $self=shift;
	my $file=shift;
	my %hash;
	open(FILE,"<",$file) || die "cannot open $file!\n$!\nexiting...\n";
	until(eof(FILE)){
		my $h1=<FILE>;
		my $s1=<FILE>;
		my $h2=<FILE>;
		my $s2=<FILE>;
		chomp $h1;
		chomp $s1;
		chomp $h2;
		chomp $s2;
		if(defined($hash{$h1})){
			die "Found $h1 in $file twice. This is NOT allowed.\n";
		}else{
			$hash{$h1}={};
			$hash{$h1}{'s1'}=$s1;
			$hash{$h1}{'h2'}=$h2;
			$hash{$h1}{'s2'}=$s2;
		}
	}
	close FILE;
	return \%hash;
}

sub LoadFasta {
	my $self=shift;
	my $file=shift;
	my %hash;
	my $header;
	open(FILE,"<",$file) || die "cannot open $file!\n$!\nexiting...\n";
	while(<FILE>){
		chomp $_;
		if($_=~m/\>/){
			$header=$_;
			$header=~s/\>//;
			$header=~s/\s.+//;
			$header=~s/\:.+//;
#			$header=~s/\|.+//;
		}else{
			if($hash{$header}){
				$hash{$header}.=$_;
			}else{
				$hash{$header}=$_;
			}
		}
	}
	close FILE;
	return \%hash;
}

sub totalSumOfSquares {
	my $self=shift;
	my @array=@{$_[0]};
	my $mean=mean($self,@array);
	my $TSS=0;
	foreach my $value (@array){
		$TSS+=(($value-$mean)**2);
	}
	return $TSS;
}
sub residualSumOfSquares {
#my $RSS=Tools->residualSumOfSquares($nr1,$nr2);
# here, we calculate the RSS, assuming that the user has done the work for us
# of populating two ordered arrays - one for the predicted values, and one for the observed value
	my $self=shift;
	my @predicted=@{$_[0]};
	my @actual=@{$_[1]};
	my $RSS=0;
	for(my$i=0;$i<=$#predicted;$i++){
		$RSS+=(($actual[$i]-$predicted[$i])**2);
	}
	return $RSS;
}

sub SearchByProfile {
	my $self=shift;
	my @profile=@{$_[0]};
	my @data=@{$_[1]};
	my $cutoff=$_[2];
	my $type=$_[3];
	my @HitIDs;
	for(my $i=0;$i<=$#data;$i++){
		my @D=@{$data[$i]};
		my $id=shift @D;
		my $c;
		if($type eq "pearson"){
			$c=pearsonsR($self,\@profile,\@D);
		}elsif($type eq "spearman"){
			$c=spearmansRho($self,\@profile,\@D);
		}elsif($type eq "kendall"){
			$c=kendallsTau($self,\@profile,\@D);
		}elsif($type eq "csim"){
			$c=csim($self,\@profile,\@D);
		}else{
			die "Unrecognized correlation type! Acceptable: \"pearson\"\n\n";
		}
		push @HitIDs, $id if $c>$cutoff;
	}
	return \@HitIDs;
}

sub getHaystackPhase {
	my $self=shift;
	my $item=shift;
	my $phase=undef;
	if($item=~m/ct(\d\d)\-.+/){
		$phase=$1;
	}elsif($item=~m/ph\_(\d+)$/){
		$phase=$1;
	}
	return sprintf("%02s",$phase);
}

sub getAllvsAll {
	my $self=shift;
	my @array=@{$_[0]}; ### Array of array refs
	my $type=$_[1]; ### valid: "pearson",
	my %AvA;
	for(my$i=0;$i<=$#array;$i++){
		my ($ID,@data)=@{$array[$i]};
		$AvA{$ID}={};
		for(my$j=0;$j<=$#array;$j++){
			next if $j==$i;
			my ($nID,@nData)=@{$array[$j]};
			if($type eq "pearson"){
				$AvA{$ID}{$nID}=pearsonsR($self,\@data,\@nData);
			}elsif($type eq "spearman"){
				$AvA{$ID}{$nID}=spearmansRho($self,\@data,\@nData);
			}elsif($type eq "kendall"){
				$AvA{$ID}{$nID}=kendallsTau($self,\@data,\@nData);
			}elsif($type eq "csim"){
				$AvA{$ID}{$nID}=csim($self,\@data,\@nData);
			}else{
				die "Unrecognized correlation type! Acceptable: \"pearson\"\n\n";
			}
		}
	}
	return \%AvA;
}

sub _createEdgeList {
	my %AvsA=%{$_[0]};
	my $cutoff=$_[1];
	my %Edges;
	my $n=0;
	foreach my $n1 (keys %AvsA){
		foreach my $n2 (keys %{$AvsA{$n1}}){
			if($AvsA{$n1}{$n2}>=$cutoff){
				$n++;
				my @edges=sort {$a cmp $b} ($n1,$n2);
				my $edge=join("-",@edges);
				$Edges{$edge}=$AvsA{$n1}{$n2};
			}else{
			}
		}
	}
	#warn "No significant interactions stronger than cutoff of $cutoff.\n" if $n==0;
	return \%Edges;
}

sub _getMaxEdge {
	my %Edges=%{$_[0]};
	my $max=0;
	my $current=undef;
	foreach my $edge (keys %Edges){
		if($Edges{$edge}>$max){
			$max=$Edges{$edge};
			$current=$edge;
		}
	}
	return $current;
}

sub _cluster {
	my %AvsA=%{$_[0]};
	my @nodes=@{$_[1]};
	my $cutoff=$_[2];
	my @Nedges;
	my @Nnodes;
	#### this should return a list of edges, so that the edge master list can be edited appropriately
	foreach my $node (@nodes) {
		my %SubNodes=%{$AvsA{$node}};
		foreach my $subnode (keys %SubNodes){
			if($AvsA{$node}{$subnode}>=$cutoff){
				my @edge=sort {$a cmp $b} ($node,$subnode);
				my $edge=join("-",@edge);
				push @Nedges, $edge;
				push @Nnodes, $subnode;
			}
		}
	}
	return [\@Nedges,\@Nnodes];
}

sub parseGffCommentToHash {
	my $self=shift;
	my $comment=shift;
	my @components=split(/\;/,$comment);
	my %hash;
	foreach my $comp (@components){
		my ($key,$value,@theRest)=split(/\=/,$comp);
		die "Error parsing GFF comment:\n$comment\nHad more than one \'=\' in a single field!\n" if defined $theRest[0];
		$hash{$key}=$value;
	}
	return \%hash;
}

sub clusterAllvsAll {
	my $self=shift;
	my %AvsA=%{$_[0]};
	my $cutoff=$_[1];
	my @nodeList=keys %AvsA;
	my %edges=%{_createEdgeList(\%AvsA,$cutoff)};
	my @master;
	my %used;
	while((scalar(keys(%edges)))>0){
		my $nremaining=scalar(keys(%edges));
		warn "clustering... $nremaining edges remain\n";
		my $bestEdge=_getMaxEdge(\%edges);
		my @cluster=split(/\-/,$bestEdge);
		delete($edges{$bestEdge});
		while(1){   #### !!!!!
			my $startSize=scalar(@cluster);
			my $r=_cluster(\%AvsA,\@cluster,$cutoff);
			my @newEdges=@{$$r[0]};
			my @newNodes=@{$$r[1]};
			for(my$i=0;$i<=$#newEdges;$i++){
				my $nEdge=$newEdges[$i];
				my $nNode=$newNodes[$i];
				if((defined($edges{$nEdge}))){
					push @cluster, $nNode unless defined($used{$nNode});
					delete($edges{$nEdge});
					$used{$nNode}=1;
				}else{
				}
			}
			my $endSize=scalar(@cluster);
			last if $endSize==$startSize;
		}
		push @master, \@cluster;
	}
	return \@master;
}

sub normalize {
	my $self=shift;
	my @array=@{$_[0]};
	my $max=max($self,@array);
	for(my$i=0;$i<=$#array;$i++){
		$array[$i]=$array[$i]/$max;
	}
	return \@array;
}

sub csim {	
	my $self=shift;
	my ($array1ref,$array2ref) = @_;
	my $c=Statistics::RankCorrelation->new($array1ref,$array2ref,sorted=>1);
	my $CDet=$c->csim;
	return $CDet;
}

sub kendallsTau {
	my $self=shift;
	my ($array1ref,$array2ref) = @_;
	my $c=Statistics::RankCorrelation->new($array1ref,$array2ref,sorted=>1);
	my $CDet=$c->kendall;
	return $CDet;
}

sub spearmansRho {
	my $self=shift;
	my ($array1ref,$array2ref) = @_;
	my $c=Statistics::RankCorrelation->new($array1ref,$array2ref,sorted=>1);
	my $CDet=$c->spearman;
	return $CDet;
}
sub pearsonsR {
	my $self=shift;
	my ($array1ref,$array2ref) = @_;
#	my @arr1=@$array1ref;
#	my @arr2=@$array2ref;
	my @arr1=@{$_[0]};
	my @arr2=@{$_[1]};
	my $mean1=mean($self,@arr1);
	my $mean2=mean($self,@arr2);
	die "Cannot calculate pearson's R for arrays of unequal size!\n\n" unless ((scalar(@arr1))==(scalar(@arr2)));
	my ($N,$SQR1,$SQR2)=(0,0,0);
	for(my $i=0;$i<=$#arr1;$i++){
		my $a1=$arr1[$i];
		my $a2=$arr2[$i];
		my $v1=($a1-$mean1);
		my $v2=($a2-$mean2);
		$N+=($v1*$v2);
		$SQR1+=$v1**2;
		$SQR2+=$v2**2;
	}
	my $D=(sqrt($SQR1))*(sqrt($SQR2))+1e-6;	
	my $r=$N/$D;
	return $r;
}
sub clusterAllvsAllSimple {
	my $self=shift;
	my %AvsA=%{$_[0]};
	my $cutoff=$_[1];
	my @nodes=keys %AvsA;
	for(my$i=0;$i<=$#nodes;$i++){
		for(my$j=0;$j<=$#nodes;$j++){
			if($AvsA{$nodes[$i]}{$nodes[$j]}<$cutoff){
				delete($AvsA{$nodes[$i]}{$nodes[$j]});
			}
		}
	}
	### look at that. clustered.
	### now just need to print it
	my %usedEdges;
	foreach my $node (@nodes){
		my @edges;
		foreach my $N (keys %{$AvsA{$node}}){
			my @edge=sort {$a cmp $b} ($node,$N);
			my $edge=join("-",@edge);
		}
	}
}
sub FisherYatesShuffle {
	my $self=shift;
	my $array = shift;
      my $i;
	for ($i = @$array; --$i; ) {
		my $j = int rand ($i+1);
		next if $i == $j;
		@$array[$i,$j] = @$array[$j,$i];
	}
	return $array;
}

sub pearsonCorrelation {  ### legacy, for old scripts with old calls
	my $self=shift;
	my ($array1ref,$array2ref) = @_;
	my @arr1=@$array1ref;
	my @arr2=@$array2ref;
	my $mean1=mean($self,@arr1);
	my $mean2=mean($self,@arr2);
	die "Cannot calculate pearson's R for arrays of unequal size!\n\n" unless ((scalar(@arr1))==(scalar(@arr2)));
	my $numerator=0;
	for(my $i=0;$i<=$#arr1;$i++){
		$numerator+=(($arr1[$i]-$mean1)*($arr2[$i]-$mean2));
	}
	my $sqres1=0;
	my $sqres2=0;
	for(my$i=0;$i<=$#arr1;$i++){
		$sqres1+=($arr1[$i]-$mean1)**2;
		$sqres2+=($arr2[$i]-$mean2)**2;
	}
	my $denominator=(sqrt($sqres1))*(sqrt($sqres2));
	my $r=$numerator/$denominator;
	return $r;
}

sub _covariance {
	my @array1=@{$_[0]};
	my @array2=@{$_[1]};
	my $result=0;
	for(my$i=0;$i<=$#array1;$i++) {
		$result+=$array1[$i]*$array2[$i];
	}
	$result /= scalar(@array1);
	$result -= mean(@array1)*mean(@array2);
}

sub fisherYates {
	my $self=shift;
	my @array=@{$_[0]};
	for(my$i=$#array;$i>=1;$i--){
		my $j=int(rand($i+1));
		my $k=$array[$j];
		$array[$j]=$array[$i];
		$array[$i]=$k;
	}
	return \@array;
}

sub revcomp {
	my $self=shift;
	my $seq=shift;
	$seq=~tr/ATCGatcg/TAGCtagc/;
	$seq=reverse $seq;
	return $seq;
}

sub _getTranslationHash {
	my %translationHash =
         (gca => "A", gcg => "A", gct => "A", gcc => "A", gcn => "A",
          tgc => "C", tgt => "C",
          gat => "D", gac => "D",
          gaa => "E", gag => "E",
          ttt => "F", ttc => "F",
          gga => "G", ggg => "G", ggc => "G", ggt => "G", ggn => "G",
          cat => "H", cac => "H",
          ata => "I", att => "I", atc => "I",
          aaa => "K", aag => "K",
          cta => "L", ctg => "L", ctt => "L", ctc => "L", ctn => "L", tta => "L", ttg => "L",
          atg => "M",
          aat => "N", aac => "N",
          cca => "P", cct => "P", ccg => "P", ccc => "P", ccn => "P",
          caa => "Q", cag => "Q",
          cga => "R", cgg => "R", cgc => "R", cgt => "R", cgn => "R",
          aga => "R", agg => "R",
          tca => "S", tcg => "S", tcc => "S", tct => "S", tcn => "S",
          agc => "S", agt => "S",
          aca => "T", acg => "T", acc => "T", act => "T", acn => "T",
          gta => "V", gtg => "V", gtc => "V", gtt => "V", gtn => "V",
          tgg => "W",
          tat => "Y", tac => "Y",
          tag => "*", taa => "*", tga => "*");
	    return \%translationHash;
}

sub getOrf {
	my $self         =shift;
	my $directStrand =shift;
	my $sequenceTitle=shift;
	my $MINORF       =shift;
	my %translationHash=%{_getTranslationHash()};
	$directStrand=~ tr/GATCN/gatcn/; # converts to lowercase
	my $reverseComplement = reverse $directStrand;
      $reverseComplement      =~ tr/gatcn/ctagn/;
	my @arrayOfORFs = ();
	my @arrayOfTranslations = ();
      #Execute outer loop twice, once for direct strand and once for reverse complement
	for(my $i = 0; $i < 2; $i = $i + 1) { # start of outer loop
      	my $sequenceEntry = "";
		my $strand = "";
		if ($i == 0) {
			$sequenceEntry  = $directStrand;
			$strand = "+";
		}else{
			$sequenceEntry  = $reverseComplement;
			$strand = "-";
		}
		my @startsRF1 =();
		my @startsRF2 =();
		my @startsRF3 =();
		my @stopsRF1 = ();
		my @stopsRF2 = ();
		my @stopsRF3 = ();
		while ($sequenceEntry =~ m/atg/gi){
			my $matchPosition = pos($sequenceEntry) - 3;
			if(($matchPosition % 3) == 0) {
				push (@startsRF1, $matchPosition);
			}
			elsif((($matchPosition + 2) % 3) == 0) {
				push (@startsRF2, $matchPosition);
			}
			else{
				push (@startsRF3, $matchPosition);
			}
		}
		while ($sequenceEntry =~ m/tag|taa|tga/gi){
			my $matchPosition = pos($sequenceEntry);
			if(($matchPosition % 3) == 0) {
				push (@stopsRF1, $matchPosition);
                       }
			elsif((($matchPosition + 2) % 3) == 0) {
				push (@stopsRF2, $matchPosition);
			}
			else{
				push (@stopsRF3, $matchPosition);
			}
		}
		my $codonRange = "";
		my $startPosition = 0;
		my $stopPosition = 0;

		@startsRF1 = reverse(@startsRF1);
		@stopsRF1 = reverse(@stopsRF1);
		while (scalar(@startsRF1) > 0) {
			$codonRange     = "";
			$startPosition  = pop(@startsRF1);
			if($startPosition < $stopPosition) {
				next;
			}
			while(scalar(@stopsRF1) > 0) {
				$stopPosition = pop(@stopsRF1);
				if ($stopPosition > $startPosition) {
					last;
				}
			}
			if($stopPosition <= $startPosition) {
				$stopPosition = length($sequenceEntry) - (length($sequenceEntry) % 3);
				$codonRange = $strand . "1|" . $startPosition . "|" . $stopPosition;
				my $this=$stopPosition-$startPosition;
				push (@arrayOfORFs, $codonRange) if $this > $MINORF;
				last;
			}else{
				$codonRange = $strand . "1|" . $startPosition . "|" . $stopPosition;
				my $this=$stopPosition-$startPosition;
				push (@arrayOfORFs, $codonRange) if $this > $MINORF;
			}
		}

		$stopPosition = 0;
		@startsRF2 = reverse(@startsRF2);
		@stopsRF2 = reverse(@stopsRF2);
		while (scalar(@startsRF2) > 0) {
			$codonRange     = "";
			$startPosition  = pop(@startsRF2);
			if($startPosition < $stopPosition) {
				next;
			}
			while(scalar(@stopsRF2) > 0) {
				$stopPosition = pop(@stopsRF2);
				if ($stopPosition > $startPosition) {
					last;
				}
			}
			if($stopPosition <= $startPosition) {
				$stopPosition = length($sequenceEntry) - ((length($sequenceEntry) + 2) % 3);
				$codonRange = $strand . "2|" . $startPosition . "|" . $stopPosition;
				my $this=$stopPosition-$startPosition;
				push (@arrayOfORFs, $codonRange) if $this > $MINORF;
				last;
			}else{
				$codonRange = $strand . "2|" . $startPosition . "|" . $stopPosition;
				my $this=$stopPosition-$startPosition;
				push (@arrayOfORFs, $codonRange) if $this > $MINORF;
			}
		}
		$stopPosition = 0;
		@startsRF3 = reverse(@startsRF3);
		@stopsRF3 = reverse(@stopsRF3);
		while (scalar(@startsRF3) > 0) {
			$codonRange     = "";
			$startPosition  = pop(@startsRF3);
			if($startPosition < $stopPosition) {
				next;
			}
			while(scalar(@stopsRF3) > 0) {
				$stopPosition = pop(@stopsRF3);
				if ($stopPosition > $startPosition) {
					last;
				}
			}
			if($stopPosition <= $startPosition) {
				$stopPosition = length($sequenceEntry) - ((length($sequenceEntry) + 1) % 3);
				$codonRange = $strand . "3|" . $startPosition . "|" . $stopPosition;
				my $this=$stopPosition-$startPosition;
				push (@arrayOfORFs, $codonRange) if $this > $MINORF;
				last;
			}else{
				$codonRange = $strand . "3|" . $startPosition . "|" . $stopPosition;
				my $this=$stopPosition-$startPosition;
				push (@arrayOfORFs, $codonRange) if $this > $MINORF;
			}
		}
	}       #end of outer loop

	#Use a loop to examine each ORF description in @arrayOfORFs. Translate ORFs that are
	#longer than $MINORF, and add those translations to @arrayOfTranslations
	foreach(@arrayOfORFs)   {
		#Use the matching operator to copy parts of the ORF range into variables.
		#The strand will be stored in $1, the reading frame in $2, the start of
		#the ORF in $3, and the end of the ORF in $4.
		$_ =~ m/([\+\-])(\d)\|(\d+)\|(\d+)/;
		next if (($4 - $3) < $MINORF);
		my $ORFsequence = "";
		if ($1 eq "+") {
			$ORFsequence    = substr($directStrand, $3, $4 - $3);
		}
		else {
			$ORFsequence    = substr($reverseComplement, $3, $4 - $3);
		}
		#Now use a for loop to translate the ORF sequence codon by codon.
		#The amino acids are added as elements to the array @growingProtein.
		my @growingProtein = ();
			for (my $i = 0; $i <= (length($ORFsequence) - 2); $i = $i + 3) {
			my$codon = substr($ORFsequence, $i, 3);
			if(exists( $translationHash{$codon} )){
				push (@growingProtein, $translationHash{$codon});
			}else{
				push (@growingProtein, "X");
			}
		}
		my $joinedAminoAcids = join("",@growingProtein);
		push (@arrayOfTranslations, $joinedAminoAcids);
#		push @arrayOfTranslations, $ORFsequence;	
	}
	return (\@arrayOfORFs,\@arrayOfTranslations);
}

1;
