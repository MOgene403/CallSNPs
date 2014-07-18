#!/usr/bin/perl
use warnings;
use strict;
use threads;
use FindBin;
use lib "$FindBin::Bin/Lib";
use Tools;

my $usage = "perl $0 <input VCF> <cutoff for including alternative bases (frequency)> <output CSV>\n\n";
die $usage unless $#ARGV==2;
my $file=$ARGV[0];
my $cut =$ARGV[1];
my $outF=$ARGV[2];

my @Input = @{Tools->LoadFile($file)};
my @output;

foreach my $line (@Input){
	my %line=parseVCF($line);
#Notch_3_HC      1514    .       A       G,T,X   0       .       DP=28329;I16=8212,19947,20,37,1023967,37747939,2006,73002,1404317,70080057,2528,112756,471290,10195328,960,20104;QS=0.998061,0.001880,0.000059,0.000000;VDB=3.929859e-01;RPB=-2.551181e-01      PL      0,255,255,255,255,255,255,255,255,255
#Notch_3_HC      1587    .       C       T,X     0       .       DP=14643;I16=1934,12587,3,15,540693,20391583,623,22221,698795,34197409,802,35920,230168,4534778,265,5131;QS=0.998836,0.001164,0.000000,0.000000;VDB=1.619116e-01;RPB=1.177923e+00       PL      0,255,255,255,255,255
#Notch_3_HC      206     .       T       C,X     0       .       DP=6211;I16=2797,2859,271,196,200591,7229181,16142,571640,282413,14105661,20858,933236,101978,2266550,8989,199349;QS=0.925546,0.074454,0.000000,0.000000;VDB=2.007198e-01;RPB=-1.554554e-01     PL      0,255,255,255,255,255
	my @I16=@{$line{I16}};
	my $tgr = $I16[0]+$I16[1];
	my $tgnr= $I16[2]+$I16[3];
	my $out = $line{chr}.",".$line{pos}.",".$line{ref}.",".$line{DP}.",".$tgr.",".$tgnr;
	my @Q=@{$line{QS}};
	my @A=@{$line{alt}};
	for(my$i=1;$i<=$#Q;$i++){
		if($Q[$i]>=$cut){
			$out.=",".$A[$i-1].",".$Q[$i];
		}else{
		}
	}
	push @output, $out;
}
Tools->printToFile($outF,\@output);

sub parseVCF {
	my $line=shift;
	my @line=split(/\t/,$line);
	my %line;
	$line{chr}=$line[0];
	$line{pos}=$line[1];
	$line{ref}=$line[3];
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
