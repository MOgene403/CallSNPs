#!/usr/bin/perl
use warnings;
use strict;
use threads;
use FindBin;
use lib "$FindBin::Bin/Lib";
use Thread::Queue;
use Tools;
use Configuration;

use threads;
use threads::shared;


my $usage="perl $0 <Config>\n\n";
die $usage unless $#ARGV==0;

my $q=Thread::Queue->new();
my $Config = Configuration->new($ARGV[0]);

my @lines = $Config->getAll("GROUPS");

my $bwa = $Config->get("PATHS","bwa");

my %Rs;
foreach my $g (@lines){
	my @R=split(/\,/,$Config->get("GROUPS",$g));
	map {$Rs{$_}=1} @R;
}

foreach my $r (keys %Rs){
	my $p = $Config->get("DIRECTORIES","References")."/".$r;
	my $command = $bwa." index $p";
	`$command`;
}

