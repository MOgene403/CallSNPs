#!/usr/bin/perl
use warnings;
use strict;
use threads;
use FindBin;
use Thread::Queue;
use lib "$FindBin::Bin/Lib";
use Configuration;
use Tools;


my $usage = "perl $0 <Configuration>\n\n";
die $usage unless $#ARGV==0;

my $config = Configuration->new($ARGV[0]);
my @Groups = $config->getAll("GROUPS");
my $OutBase= $config->get("DIRECTORIES","Output");
foreach my $grp (@Groups){
	my @indicies = split(/\,/,$config->get("GROUPS",$grp));
	foreach my $index (@indicies){
		my $vcfFile = $OutBase."/".$grp."/$grp.vs.$index.Alignments.raw.vcf";
		warn "Can't find $vcfFile\n" unless -e $vcfFile;
		my $outFile = $OutBase."/".$grp."/$grp.byDist.csv";
		parseForDistribs($vcfFile,$outFile,$config->get("OPTIONS","snpRate"),$config->get("OPTIONS","minZ"));
	}
}

sub parseForDistribs {
	my $in_File=shift;
	my $outFile=shift;
	my $Scut = shift;
	my $Zcut = shift;
	my @Input = @{Tools->LoadFile($in_File)};
	my @insertfreqs;
	my @output;
	my @OL;
	my @OC;
	warn $in_File."\n";
	my $main=0;
	foreach my $line (@Input){
		next if $line=~m/^\#/;
		my %line=%{parseVCF($line)};
#Notch_3_HC      1514    .       A       G,T,X   0       .       DP=28329;I16=8212,19947,20,37,1023967,37747939,2006,73002,1404317,70080057,2528,112756,471290,10195328,960,20104;QS=0.998061,0.001880,0.000059,0.000000;VDB=3.929859e-01;RPB=-2.551181e-01      PL      0,255,255,255,255,255,255,255,255,255
#Notch_3_HC      1587    .       C       T,X     0       .       DP=14643;I16=1934,12587,3,15,540693,20391583,623,22221,698795,34197409,802,35920,230168,4534778,265,5131;QS=0.998836,0.001164,0.000000,0.000000;VDB=1.619116e-01;RPB=1.177923e+00       PL      0,255,255,255,255,255
#Notch_3_HC      206     .       T       C,X     0       .       DP=6211;I16=2797,2859,271,196,200591,7229181,16142,571640,282413,14105661,20858,933236,101978,2266550,8989,199349;QS=0.925546,0.074454,0.000000,0.000000;VDB=2.007198e-01;RPB=-1.554554e-01     PL      0,255,255,255,255,255
		my @I16=@{$line{I16}};
		my $tgr = $I16[0]+$I16[1];
		my $tgnr= $I16[2]+$I16[3];
		my @Q=@{$line{QS}};
		my @A=@{$line{alt}};
		my $jq=shift @Q;
		my $rq=shift @Q;
		my $rb=shift @A;
		my $pass=0;
		my $out = $line{chr}.",".$line{pos}.",".$rb.",".$line{DP}.",".$rq.",".$Q[0].",".$A[0];
		for(my$i=0;$i<=$#Q;$i++){
			if($Q[$i] >= $Scut){
#				push @output, $line;
		push @output, $out;
				$main++;
			}elsif($Q[$i]==0.0){
			}else{
				my $minRat = getQRatios($Scut,\@Q);
#				warn $minRat."\n".$line."\n";
				my @OD;
				push @OD, $line{pos};
				push @OD, $A[$i];
				push @OD, $Q[$i];
				push @insertfreqs, $Q[$i];
				push @OC, \%line;
				push @OL, \@OD;
				last;
			}
		}
	}
	my $mean = Tools->mean(@insertfreqs);
	my $stdev= Tools->stdev(@insertfreqs);
	print "mean: $mean\n";
	print "stdev: $stdev\n";
	my $hits=0;
	foreach my $ref (@OC){
		my %line = %{$ref};
		my $out = 
		my @Q = @{$line{QS}};
		my @A=@{$line{alt}};
		my $jq=shift @Q;
		for(my$i=0;$i<=$#Q;$i++){
			next if $Q[$i] >= $Scut;
			next if $Q[$i] == 0;
			my $minRat =  getQRatios($Scut,\@Q);
			my $N = $A[$i];
			my $Z = ($Q[$i] - $mean)/$stdev;
			if($Z>$Zcut){
				#push @output, $line{line}."\t".$Q[$i]."\t".$mean."\t".$stdev."\t".$Z;
				push @output, $line{chr}.",".$line{pos}.",".$A[0].",".$line{DP}.",".$Q[0].",".$Q[$i].",".$N.",".$Z.",".$minRat;
#				push @output, $N;
				$hits++;
			}
		}
	}
#	for(my$r=0;$r<=$#OL;$r++){
	# proc candidate non-zero variance locations
#		my @D = @{$OL[$r]};
#		my $pos = shift @D;
#		for(my$i=0;$i<=$#D-1;$i+=2){
#			my $N = $D[$i];
#			my $Q = $D[$i+1];
#			my $Z = ($Q-$mean)/$stdev;
#			if($Z>$Zcut){
#				push @output, $OC[$r]."\t".$Q."\t".$mean."\t".$stdev."\t".$Z;
#				$hits++;
#			}
#		}
#	}
	warn "$main locations exceeded base SNP call rate, $hits exceeded one-tail-Z>=2.33 (p<0.01)\n";
	Tools->printToFile($outFile,\@output);
}

sub getQRatios {
	my $cut=shift;
	my @Q=@{$_[0]};
	my @R;
	for(my$i=0;$i<=$#Q;$i++){
		for(my$j=0;$j<=$#Q;$j++){
			next if $Q[$j] == 0;
			next if $Q[$i] == 0;
			next if $Q[$i] > $cut;
			next if $Q[$j] > $cut;
			my $rat = $Q[$i]/$Q[$j];
			push @R, $rat;
		}
	}
	@R = sort {$b <=> $a} @R;
	return $R[0];
}

sub parseVCF {
	my $line=shift;
	my @line=split(/\t/,$line);
	my %line;
	$line{chr}=$line[0];
	$line{pos}=$line[1];
	$line{ref}=$line[3];
	$line{line}=$line;
	$line{alt}=[ split(/\,/,$line[4]) ];
	my %info=%{parseInfo($line[7])};
	$line{DP}=$info{"DP"};
	$line{I16}=$info{"I16"};
	$line{QS}=$info{"QS"};
	return \%line;
}

sub parseInfo {
	my $info=shift;
	my @info=split(/\;/,$info);
	my %info;
	foreach my $item (@info){
		my ($key,$val)=split(/\=/,$item);
		if($key eq "I16"){
			$info{$key}=[split(/\,/,$val)];
		}elsif($key eq "DP"){
			$info{$key}=$val;
		}elsif($key eq "QS"){
			$info{$key}=[split(/\,/,$val)];
		}else{
		}
	}
	return \%info;
}

#
#Notch_3_HC      1514    .       A       G,T,X   0       .       DP=28329;I16=8212,19947,20,37,1023967,37747939,2006,73002,1404317,70080057,2528,112756,471290,10195328,960,20104;QS=0.998061,0.001880,0.000059,0.000000;VDB=3.929859e-01;RPB=-2.551181e-01      PL      0,255,255,255,255,255,255,255,255,255
#Notch_3_HC      1587    .       C       T,X     0       .       DP=14643;I16=1934,12587,3,15,540693,20391583,623,22221,698795,34197409,802,35920,230168,4534778,265,5131;QS=0.998836,0.001164,0.000000,0.000000;VDB=1.619116e-01;RPB=1.177923e+00       PL      0,255,255,255,255,255
#Notch_3_HC      206     .       T       C,X     0       .       DP=6211;I16=2797,2859,271,196,200591,7229181,16142,571640,282413,14105661,20858,933236,101978,2266550,8989,199349;QS=0.925546,0.074454,0.000000,0.000000;VDB=2.007198e-01;RPB=-1.554554e-01     PL      0,255,255,255,255,255
#Notch_3_HC      823     .       C       G,T,A   0       .       DP=22419;I16=13616,8574,18,12,809817,30045033,1070,39182,1105853,55173431,1336,59566,383557,8178791,507,10373;QS=0.998679,0.001048,0.000222,0.000051;VDB=1.024477e-01;RPB=4.940619e-01  PL      0,255,255,255,255,255,255,255,255,255
#Notch_3_HC      223     .       T       G,C,X   0       .       DP=7297;I16=2775,3660,18,1,188615,6009625,338,6332,319350,15865272,945,47011,101923,2109925,251,4557;QS=0.998211,0.001588,0.000201,0.000000;VDB=1.498251e-01;RPB=4.807658e+00   PL      0,255,255,255,255,255,255,255,255,255
#
