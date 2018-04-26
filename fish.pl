#!/usr/bin/perl
#Informatic Biology departments of Beijing Genomics Institute (BGI) 
#Modified by Wanfei Liu
use strict;
use Getopt::Long;
my %opts;
my $program=`basename $0`;
chomp $program;
my $usage=<<USAGE; #******* Instruction of this program *********# 

Author:  Fan wei, <fanw\@genomics.org.cn>, new year 2006

Program: fishing in one file according to baits in another file

Usage: $program  <bait_file fish_file>
    -tbait<type>		type of bait file
	-tfish<type>		type of fish file
	-bait <colum>       colum as bait (1,2,3)
	-fish <colum>       colum as fish (1,2,3)
	-contrary			fish those not in the bait
	-details            output the middle results to screen
	-help               output help information

Comment: This is a frequently used program, which deals with the most frequetly used 
file format .fa, .psl, and .list(just like excel tables which have several colums and rows).
It gets the baits(key that is universal in both bait and fish file) from the bait file
at first, and then searches the baits in the fish file. If found, then retrive that item
from the fish file. If -contrary is specifed, retrive those items which are not exist
in the bait file.

USAGE

GetOptions(\%opts, "bait:s","fish:s","tbait:s","tfish:s","details!","contrary","help!");
die $usage if ( @ARGV==0 || defined($opts{"help"}));

#****************************************************************#
#--------------------Main-----Function-----Start-----------------#
#****************************************************************#
my $bait_file=shift;
my $fish_file=shift;
my %bait;


if ( $opts{tbait} eq 'psl'  || ( ! exists $opts{tbait} && $bait_file=~/.psl$/ ) ) {
	read_psl($bait_file,\%bait);
}elsif($opts{tbait} eq 'fa' || ( ! exists $opts{tbait} && $bait_file=~/.fa$/ ) ){
	read_fa($bait_file,\%bait);
}else{
	read_list($bait_file,\%bait);
}

print STDERR "read bait done\n" if (exists $opts{details});

if ($opts{tfish} eq 'psl' || ( ! exists $opts{tfish} && $fish_file=~/.psl$/ ) ) {
	out_psl($fish_file,\%bait);
}elsif($opts{tfish} eq 'fa' || ( ! exists $opts{tfish} && $fish_file=~/.fa$/ ) ){
	out_fa($fish_file,\%bait);
}else{
	out_list($fish_file,\%bait);
}

print STDERR "Fish out done\n" if (exists $opts{details});
#****************************************************************#
#------------------Children-----Functions-----Start--------------#
#****************************************************************#

sub read_list{
	my $file=shift;
	my $bait_hp=shift;
	my $bait_colum=1;
	my @bait_colum=(1);
	if (defined $opts{bait}) {
		$bait_colum=$opts{bait};
		@bait_colum=split /\,/,$bait_colum;
	}
	open(IN,$file)||die("fail to open $file\n");
	while (<IN>) {
		my @temp=split(/\s+/,$_);
		my $key="";
		foreach my $bait (@bait_colum) {
			$key.="$temp[$bait-1]\t";
		}
		$key=~s/\t$//;
		$$bait_hp{$key}=1;
	}
	close(IN);
}

sub read_psl{
	my $file=shift;
	my $bait_hp=shift;
	my $bait_colum=10;
	my @bait_colum=(10);
	if (defined $opts{bait}) {
		$bait_colum=$opts{bait};
		@bait_colum=split /\,/,$bait_colum;
	}
	open(IN,$file)||die("fail to open $file\n");
	while (<IN>) {
		my @temp=split(/\s+/,$_);
		my $key="";
		foreach my $bait (@bait_colum) {
			$key.="$temp[$bait-1]\t";
		}
		$key=~s/\t$//;
		$$bait_hp{$key}=1;
	}
	close(IN);
}

sub read_fa{
	my $file=shift;
	my $bait_hp=shift;
	my $bait_colum=1;
	my @bait_colum=(1);
	if (defined $opts{bait}) {
		$bait_colum=$opts{bait};
		@bait_colum=split /\,/,$bait_colum;
	}
	open(IN,$file)||die("fail to open $file\n");
	$/=">";<IN>;$/="\n";
	while (<IN>) {
		my $title=$_;
		$/=">";
		my $seq=<IN>;
		chomp $seq;
		$/="\n";
		my @temp=split(/\s+/,$title);
		my $key="";
		foreach my $bait (@bait_colum) {
			$key.="$temp[$bait-1]\t";
		}
		$key=~s/\t$//;
		$$bait_hp{$key}=1;
	}
	close(IN);
}



sub out_list{
	
	my $file=shift;
	my $fish_hp=shift;
	my $fish_colum=1;
	my @fish_colum=(1);
	if (defined $opts{fish}) {
		$fish_colum=$opts{fish};
		@fish_colum=split /\,/,$fish_colum;
	}
	my $output;
	open(IN,$file)||die("fail to open $file\n");
	while (<IN>) {
		my @temp=split(/\s+/,$_);
		my $key="";
		foreach my $fish (@fish_colum) {
			$key.="$temp[$fish-1]\t";
		}
		$key=~s/\t$//;
		$output .= $_ if ( exists $$fish_hp{$key} && !exists $opts{contrary});
		$output .= $_ if ( !exists $$fish_hp{$key} && exists $opts{contrary});
	}
	close(IN);
	print $output;
}

sub out_psl{
	my $file=shift;
	my $fish_hp=shift;
	my $fish_colum=10;
	my @fish_colum=(10);
	if (defined $opts{fish}) {
		$fish_colum=$opts{fish};
		@fish_colum=split /\,/,$fish_colum;
	}
	my $output;
	open(IN,$file)||die("fail to open $file\n");
	while (<IN>) {
		my @temp=split(/\s+/,$_);
		my $key="";
		foreach my $fish (@fish_colum) {
			$key.="$temp[$fish-1]\t";
		}
		$key=~s/\t$//;
		$output .= $_ if ( exists $$fish_hp{$key} && !exists $opts{contrary});
		$output .= $_ if ( !exists $$fish_hp{$key} && exists $opts{contrary});
	}
	close(IN);
	print $output;
}

sub out_fa{
	my $file=shift;
	my $fish_hp=shift;
	my $fish_colum=1;
	my @fish_colum=(1);
	if (defined $opts{fish}) {
		$fish_colum=$opts{fish};
		@fish_colum=split /\,/,$fish_colum;
	}
	my $output;
	open(IN,$file)||die("fail to open $file\n");
	$/=">";<IN>;$/="\n";
	while (<IN>) {
		my $title=$_;
		$/=">";
		my $seq=<IN>;
		chomp $seq;
		$/="\n";
		my @temp=split(/\s+/,$title);
		my $key="";
		foreach my $fish (@fish_colum) {
			$key.="$temp[$fish-1]\t";
		}
		$key=~s/\t$//;
		$output .= ">".$title.$seq if ( exists $$fish_hp{$key} && !exists $opts{contrary});
		$output .= ">".$title.$seq if ( !exists $$fish_hp{$key} && exists $opts{contrary});
	}
	close(IN);
	print $output;
}
