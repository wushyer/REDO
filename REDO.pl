#!/usr/bin/perl
#Copyright (c) BIG
#Author: Wu Shuangyang <wushy@big.ac.cn>& Liu Wanfei <liuwf@big.ac.cn>
#Date:4/12/2017
#Description:This program can get RNA editing sites based on vcf file.
#use strict;
#modified pvalue default value and movied depth cutoff before err and depth statistics 
use warnings;
my $version="1.0 version";
use Getopt::Long;
use Text::NSP::Measures::2D::Fisher::left;
use Text::NSP::Measures::2D::Fisher::twotailed;

my %opts;
GetOptions(\%opts,"g:s","v:s@","t:s","o:s","d:s","c:s","p:s","w:s","s:s","a:s","l:s","f:s","dv:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{g} or !defined $opts{v} or !defined $opts{t} or !defined $opts{o}) {
	die "************************************************************************
	Usage: REDO.pl -g genome_sequence -v variant_file -t tbl_file -o prefix_of_outfile -d reads_depth -c minimum_coverage_of_alt_allele -p alt_proportion -w window_size -s minimum_splice_anchor -a minimum_alt_distance -l llr_value -f fisher_pvalue -dv depth_variant.
	   -g: the absolute path of genome sequence.
	   -v: the absolute path of variant file/files (for example, \"-v test1.vcf -v test2.vcf -v test3.vcf ......\").
	   -t: the absolute path of tbl_file (feature table file, http://www.ncbi.nlm.nih.gov/projects/Sequin/table.html).
	   -o: the prefix of output file.
	   -d: minimum reads depth (4).
	   -c: minimum coverage of altanative allele (3).
	   -p: minimum alt proportion (0.1).
	   -w: minimum window size for calculating error rate and average depth (10).
	   -s: minimum splice_anchor size (0).
	   -a: minimum alt distance (0).
	   -l: minimum likelihood ratio (LLR) (10).
	   -f: maximum pvalue of fisher exact test (0.01).
	  -dv: the depth variant (0.2).
	   -i: identify indel (y or n, default n).
************************************************************************\n";
}
my $genomefile=$opts{g};
my @vcffile=@{$opts{v}};
my $tblfile=$opts{t};
my $prefix=$opts{o};

my $depth_cutoff=4;
if (exists $opts{d}) {
	$depth_cutoff=$opts{d};
}
my $alt_cutoff=3;
if (exists $opts{c}) {
	$alt_cutoff=$opts{c};
}
my $alt_proportion=0.1;
if (exists $opts{p}) {
	$alt_proportion=$opts{p};
}
my $window=10;
if (exists $opts{w}) {
	$window=$opts{w};
}
my $splice_anchor=2;
if (exists $opts{s}) {
	$splice_anchor=$opts{s};
}
my $alt_distance=3;
if (exists $opts{a}) {
	$alt_distance=$opts{a};
}
my $dep_variant=0.2;
if (exists $opts{dv}) {
	$dep_variant=$opts{dv};
}
my $fisher_pvalue=0.01;
if (exists $opts{f}) {
	$fisher_pvalue=$opts{f};
}
my $llr_value=10;
if (exists $opts{l}) {
	$llr_value=$opts{l};
}

my $indel="n";
if (exists $opts{i}) {
	$indel=$opts{i};
}

my ($id,%chr,@prefix,%vcf,%hydro,);

open (IN,"<$genomefile")||die("fail to open $genomefile.\n");
while(<IN>){
	next if(/^\#/); #ignore header
	chomp;
	if (/>(\S+)/) {
		$id=$1;
		next;
	}
	$chr{$id}.=$_;
}
close IN;

for (my $i=0;$i<@vcffile;$i++) {
	my $file=$vcffile[$i];
	$file=~/([^\/]+)$/;
	$file=$1;
	while ($file=~/^(.+)\.[^\.]+$/) {
		$file=~s/^(.+)\.[^\.]+$/$1/;
	}
	push (@prefix,$file);
	my $distance=undef;
	my @depth=();
	my $depth=0;
	my @err=();
	my $err=0;
	open (IN,"<$vcffile[$i]")||die("fail to open $vcffile[$i].\n");
	while(<IN>){
		next if(/^\#/); #ignore header
		chomp;
		#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SRR1063404_cp.sorted.bam
		#cp      1       .       A       .       38.9945 .       DP=3;MQ0F=0;AF1=0;AC1=0;DP4=3,0,0,0;MQ=33;FQ=-35.98
		my @list=split /\t/,$_;
		my %attr=();
		$attr{"CHROM"}=$list[0];
		$attr{"POS"}=$list[1];
		$attr{"REF"}=$list[3];
		$attr{"ALT"}=$list[4];
		my @ALT=split /\,/,$list[4];
		for (my $j=0;$j<@ALT;$j++) {
			$attr{"ALT$j"}=$ALT[$j];
		}
		$attr{"QUAL"}=$list[5];
		my @attributes = split /;/, $list[7];
		foreach my $attr ( @attributes) {
			if ($attr =~ /^(\S+)\=(\S+)$/) {
				my $c_type  = $1;
				my $c_value = $2;
				$attr{$c_type} = $c_value;
				my @c_value=split /\,/,$c_value;
				for (my $j=0;$j<@c_value;$j++) {
					$attr{"$c_type$j"}=$c_value[$j];
				}
			}elsif ($attr =~ /^(\S+)$/) {
				my $c_type  = $1;
				$attr{$c_type} = "";
			}
		}
		my @attributes1 = split /:/, $list[8];
		my @attributes2 = split /:/, $list[9];
		for (my $j=0;$j<@attributes1;$j++) {
			my $c_type  = $attributes1[$j];
			my $c_value = $attributes2[$j];
			$attr{$c_type} = $c_value;
			my @c_value=split /\,/,$c_value;
			for (my $k=0;$k<@c_value;$k++) {
				$attr{"$c_type$k"}=$c_value[$k];
			}
		}
		my $ref=0;
		my $alt=0;
		if (defined $attr{"DP4"}) {
			$ref=$attr{"DP40"}+$attr{"DP41"};
			$alt=$attr{"DP42"}+$attr{"DP43"};
		}elsif (defined $attr{"AD"}) {
			$ref=$attr{"AD0"};
			$alt=$attr{"AD1"};
		}
		#filter: 1 low quality 2 depth 
		next if ($alt+$ref<$depth_cutoff);
		if ($attr{"ALT"} eq ".") {
			if (@depth<$window) {
				push (@depth,$ref+$alt);
				my $total=0;
				for (my $j=0;$j<@depth;$j++) {
					$total+=$depth[$j];
				}
				$depth=$total/@depth;
			}else {
				shift @depth;
				push (@depth,$ref+$alt);
				my $total=0;
				for (my $j=0;$j<@depth;$j++) {
					$total+=$depth[$j];
				}
				$depth=$total/@depth;
			}
			if (@err<$window) {
				push (@err,$alt/($ref+$alt));
				my $total=0;
				for (my $j=0;$j<@err;$j++) {
					$total+=$err[$j];
				}
				$err=$total/@err;
			}else {
				shift @err;
				push (@err,$alt/($ref+$alt));
				my $total=0;
				for (my $j=0;$j<@err;$j++) {
					$total+=$err[$j];
				}
				$err=$total/@err;
			}
		#filter multiple alt
		}elsif ($attr{"ALT"} ne "." and @ALT==1) {
			#filter indel 
			if ($indel eq "n") {
				next if (length($attr{"REF"})>1 or length($attr{"ALT"})>1);
			}
			#filter: 3 alt coverage 4 alt proportion
			next if ($alt<$alt_cutoff);
			next if ($alt/($ref+$alt)<$alt_proportion);
			#filter depth variant 
			next if ($ref+$alt<$depth*$dep_variant);
			#filter possible mismatch site comparing with adjacent sites
			next if ($alt/($ref+$alt)-$err<$alt_proportion);
			if (!defined $distance) {
				$distance=$attr{"POS"};
			}else {
				#filter the RNA editing sites with short distance from each other
				if ($attr{"POS"}-$distance<$alt_distance) {
					delete $vcf{$i}{$attr{"CHROM"}}{$distance} if (exists $vcf{$i}{$attr{"CHROM"}}{$distance});
					delete $vcf{$i}{$attr{"CHROM"}}{$attr{"POS"}} if (exists $vcf{$i}{$attr{"CHROM"}}{$attr{"POS"}});
				}else {
					#a likelihood ratio (LLR) test
					my $perr=0;
					if ($err<0.0001) {
						$perr=((1-0.0001)**$ref)*((0.0001)**$alt);
					}else {
						$perr=((1-$err)**$ref)*(($err)**$alt);
					}
					my $pedit=(($ref/($ref+$alt))**$ref)*((1-$alt/($ref+$alt))**$alt);
					my $llr=0;
					if ($perr==0 or $pedit==0) {
						$llr=1000;
					}else {
						$llr=&log10($pedit/$perr);
					}
					#fisher exact test
					my $pfisher=0;
					$pfisher=&fisher_left(0,$ref+$alt,(sprintf "%.0f",$alt-($ref+$alt)*$err),$ref);
					next if ($llr<$llr_value or $pfisher>$fisher_pvalue);
					$vcf{$i}{$attr{"CHROM"}}{$attr{"POS"}}=$attr{"REF"}."\t".$attr{"ALT"}."\t".($alt/($ref+$alt))."\t$err\t".($ref+$alt)."\t$depth\t$ref\t$alt\t".($attr{"POS"}-$distance)."\t$llr\t$pfisher";
				}
				$distance=$attr{"POS"};
			}
		}
	}
	close IN;
}

my $acc=0;
my $name="";
my $str="";
my $type="";
my $product="";
my $strand="";
my $exon_num=0;
my $exon_start="";
my $exon_end="";
open (OUT,">$prefix.out")||die("fail to open $prefix.out.\n");
print OUT "#Chr\tStrand\tName\tGenome_pos\tGene_pos\tAA_pos\tSample_num\tPhase\t".(join "\t",@prefix)."\tExon_num\tExon_start\tExon_end\tFunction\n";
for (my $k=0;$k<@prefix;$k++) {
	my $filehandle="OUT$k";
	open ($filehandle,">$prefix\_$prefix[$k].out")||die("fail to open $prefix\_$prefix[$k].out.\n");
	print $filehandle "#Chr\tStrand\tName\tGenome_pos\tGene_pos\tAA_pos\tPhase\tRef->Alt\tRefCodon->AltCodon\tRefAA->AltAA\tAltRatio\tErrorRatio\tDepth\tAverageDepth\tRefRead\tAltRead\tAdjacentAlt\tLLR\tPvalue(fisher)\tRNAeditGC\tGeneGC\tExon_num\tExon_start\tExon_end\tFunction\n";
}
open (IN,"<$tblfile")||die("fail to open $tblfile.\n");
while (<IN>) {
	chomp;
	if (/^>\S+\s+(\S+)$/) {
		$id=$1;
	}elsif (/^\<?(\d+)\t(\d+)\tgene$/) {
		if ($str ne "") {
			my @start=split /,/,$exon_start;
			my @end=split /,/,$exon_end;
			my $gene_pos=0;
			my $gene_gc=&gc($str);
			for (my $i=0;$i<@start;$i++) {
				if ($start[$i]<$end[$i]) {
					for (my $j=$start[$i];$j<=$end[$i];$j++) {
						$strand="+";
						$gene_pos++;
						#filter according to the splice anchor size
						if (@start>1) {
							if ($i==0) {
								next if ($end[$i]-$j+1<$splice_anchor);
							}elsif ($i<@start-1) {
								next if ($j-$start[$i]+1<$splice_anchor or $end[$i]-$j+1<$splice_anchor);
							}else {
								next if ($j-$start[$i]+1<$splice_anchor);
							}
						}
						my $print_str=&RNA_alt($j,$gene_pos,$strand,$str,$gene_gc,$type);
						if ($print_str ne "") {
							print OUT "$print_str";
						}
					}
				}else {
					for (my $j=$start[$i];$j>=$end[$i];$j--) {
						$strand="-";
						$gene_pos++;
						#filter according to the splice anchor size
						if (@start>1) {
							if ($i==0) {
								next if ($j-$end[$i]+1<$splice_anchor);
							}elsif ($i<@start-1) {
								next if ($start[$i]-$j+1<$splice_anchor or $j-$end[$i]+1<$splice_anchor);
							}else {
								next if ($start[$i]-$j+1<$splice_anchor);
							}
						}
						my $print_str=&RNA_alt($j,$gene_pos,$strand,$str,$gene_gc,$type);
						if ($print_str ne "") {
							print OUT "$print_str";
						}
					}
				}
			}
		}
		$type="";$name="";$strand="";$exon_num=0;$exon_start="";$exon_end="";$product="";$str="";
		$acc++;
	}elsif (/^\t\t\t(gene|locus_tag)\t(.*)$/) {
		if ($name eq "") {
			$name=$2;
		}
	}elsif (/^\<?(\d+)\t(\d+)\t(CDS|tRNA|rRNA)$/) {
		$type=$3;
		if ($type eq "CDS" or $type eq "tRNA" or $type eq "rRNA") {
			my $pos1=$1;
			my $pos2=$2;
			if ($pos1<$pos2) {
				my $substr=substr($chr{$id},$pos1-1,$pos2-$pos1+1);
				$str.=$substr;
			}else {
				my $substr=substr($chr{$id},$pos2-1,$pos1-$pos2+1);
				$substr=reverse $substr;
				$substr=~tr/AaTtGgCc/TtAaCcGg/;
				$str.=$substr;
			}
			$exon_num++;
			$exon_start.=$pos1.",";
			$exon_end.=$pos2.",";
		}
	}elsif (/^\<?(\d+)\t(\d+)$/) {
		if ($type eq "CDS" or $type eq "tRNA" or $type eq "rRNA") {
			my $pos1=$1;
			my $pos2=$2;
			if ($pos1<$pos2) {
				my $substr=substr($chr{$id},$pos1-1,$pos2-$pos1+1);
				$str.=$substr;
			}else {
				my $substr=substr($chr{$id},$pos2-1,$pos1-$pos2+1);
				$substr=reverse $substr;
				$substr=~tr/AaTtGgCc/TtAaCcGg/;
				$str.=$substr;
			}
			$exon_num++;
			$exon_start.=$pos1.",";
			$exon_end.=$pos2.",";
		}
	}elsif (/^\t\t\tproduct\t(.*)$/) {
		$product=$1;
	}
}
if ($str ne "") {
	my @start=split /,/,$exon_start;
	my @end=split /,/,$exon_end;
	my $gene_pos=0;
	my $gene_gc=&gc($str);
	for (my $i=0;$i<@start;$i++) {
		if ($start[$i]<$end[$i]) {
			for (my $j=$start[$i];$j<=$end[$i];$j++) {
				$strand="+";
				$gene_pos++;
				#filter according to the splice anchor size
				if (@start>1) {
					if ($i==0) {
						next if ($end[$i]-$j+1<$splice_anchor);
					}elsif ($i<@start-1) {
						next if ($j-$start[$i]+1<$splice_anchor or $end[$i]-$j+1<$splice_anchor);
					}else {
						next if ($j-$start[$i]+1<$splice_anchor);
					}
				}
				my $print_str=&RNA_alt($j,$gene_pos,$strand,$str,$gene_gc,$type);
				if ($print_str ne "") {
					print OUT "$print_str";
				}
			}
		}else {
			for (my $j=$start[$i];$j>=$end[$i];$j--) {
				$strand="-";
				$gene_pos++;
				#filter according to the splice anchor size
				if (@start>1) {
					if ($i==0) {
						next if ($j-$end[$i]+1<$splice_anchor);
					}elsif ($i<@start-1) {
						next if ($start[$i]-$j+1<$splice_anchor or $j-$end[$i]+1<$splice_anchor);
					}else {
						next if ($start[$i]-$j+1<$splice_anchor);
					}
				}
				my $print_str=&RNA_alt($j,$gene_pos,$strand,$str,$gene_gc,$type);
				if ($print_str ne "") {
					print OUT "$print_str";
				}
			}
		}
	}
}
close IN;
close OUT;

for (my $k=0;$k<@prefix;$k++) {
	my $filehandle="OUT$k";
	close $filehandle;
}

my @change=("A->C","A->G","A->T","C->G","C->T","C->A","G->T","G->A","G->C","T->A","T->C","T->G");
open (OUT,">$prefix.change_matrix")||die("fail to open $prefix.change_matrix.\n");
print OUT "#Prefix";
for (my $i=0;$i<@change;$i++) {
	print OUT "\t$change[$i](percent)";
}
print OUT "\n";
for (my $k=0;$k<@vcffile;$k++) {
	my $total=0;
	foreach my $alter (keys %{$change{$k}}) {
		$total+=$change{$k}{$alter};
	}
	print OUT "$prefix[$k]";
	foreach my $alter (@change) {
		my $percent=0;
		if ($total>0) {
			my $change=0;
			$change=$change{$k}{$alter} if (exists $change{$k}{$alter});
			$percent=(sprintf "%.4f",$change/$total)*100;
		}
		my $subtotal=0;
		if (exists $change{$k}{$alter}) {
			$subtotal=$change{$k}{$alter};
		}
		print OUT "\t$subtotal($percent)";
	}
	print OUT "\n";
}
close OUT;

#print hydrophobic and hydrophilic statistic
open (OUT,">$prefix.stat")||die("fail to open $prefix.stat.\n");
print OUT "#Prefix\tTotal\tOther\tPhase(1,2,3)\tSilent\tNonSilent";
for (my $i=0;$i<@change;$i++) {
	print OUT "\t$change[$i](silent)\t$change[$i](non_silent)\t$change[$i](1:hydrophobic2hydrophobic,2:hydrophilic2hydrophilic,3:hydrophobic2hydrophilic,4:hydrophilic2hydrophobic)";
}
print OUT "\n";
for (my $k=0;$k<@vcffile;$k++) {
	my $total=0;
	$total=$change{$k}{""} if (exists $change{$k}{""});
	my $other=0;
	$other=$change{$k}{""} if (exists $change{$k}{""});
	my $silent=0;
	my $non_silent=0;
	my $phase1=0;
	my $phase2=0;
	my $phase3=0;
	foreach my $alter (@change) {
		$phase1+=$hydro{$k}{$alter}{"1"}  if (exists $hydro{$k}{$alter}{"1"});
		$phase2+=$hydro{$k}{$alter}{"2"}  if (exists $hydro{$k}{$alter}{"2"});
		$phase3+=$hydro{$k}{$alter}{"3"}  if (exists $hydro{$k}{$alter}{"3"});
		$total+=$change{$k}{$alter} if (exists $change{$k}{$alter});
		$silent+=$hydro{$k}{$alter}{"silent"}  if (exists $hydro{$k}{$alter}{"silent"});
		$non_silent+=$hydro{$k}{$alter}{"non_silent"}  if (exists $hydro{$k}{$alter}{"non_silent"});
	}
	print OUT "$prefix[$k]\t$total\t$other\t$phase1,$phase2,$phase3\t$silent\t$non_silent";
	foreach my $alter (@change) {
		my $si=0;
		$si=$hydro{$k}{$alter}{"silent"} if (exists $hydro{$k}{$alter}{"silent"});
		my $no_si=0;
		$no_si=$hydro{$k}{$alter}{"non_silent"} if (exists $hydro{$k}{$alter}{"non_silent"});
		my $phopho=0;
		$phopho=$hydro{$k}{$alter}{"phobicphobic"} if (exists $hydro{$k}{$alter}{"phobicphobic"});
		my $phiphi=0;
		$phiphi=$hydro{$k}{$alter}{"philicphilic"} if (exists $hydro{$k}{$alter}{"philicphilic"});
		my $phophi=0;
		$phophi=$hydro{$k}{$alter}{"phobicphilic"} if (exists $hydro{$k}{$alter}{"phobicphilic"});
		my $phipho=0;
		$phipho=$hydro{$k}{$alter}{"philicphobic"} if (exists $hydro{$k}{$alter}{"philicphobic"});
		print OUT "\t$si\t$no_si\t$phopho,$phiphi,$phophi,$phipho";
	}
	print OUT "\n";
}
close OUT;

for (my $i=0;$i<@prefix;$i++) {
	my $num=(keys %{$change{$i}})-1;
	open (R,">$prefix\_$prefix[$i].R")||die("fail to open $prefix\_$prefix[$i].R.\n");
	print R "jpeg(\"$prefix\_$prefix[$i].jpeg\",width=12000,height=12000,res=600);\n";
	print R "nopar<-par(mfrow=c(4,4),cex.axis=1.5);\n";
	print R "A<-read.table(\"$prefix\_$prefix[$i].out\",sep=\"\\t\");\n";
	#Chr	Accession	Name	Genome_pos	Gene_pos	AA_pos	Phase	Ref->Alt	RefCodon->AltCodon	RefAA->AltAA	AltRatio	ErrorRatio	Depth	AverageDepth	RefRead	AltRead	AdjacentAlt	Exon_num	Exon_start	Exon_end	Function
	print R "plot(density(A\$V11),xlim=c(-0.2,1.2),type=\"l\",lwd=2,cex.axis=1.5,cex.lab =1.5,cex.main =1.8,col=\"red\",main=\"$prefix[$i]\",xlab=\"Alt allele proportion\",ylab=\"Density\");\n";
	print R "lines(density(A\$V12),col=\"green\");\n";
	print R "legend(\"topright\",cex=1.5,lwd=2,legend=c(\"AltAllele\",\"Error\"),col=c(\"red\",\"green\"),lty=c(1,1),bty=\"n\");\n";
	print R "plot(density(log10(A\$V13)),type=\"l\",lwd=2,cex.axis=1.5,cex.lab =1.5,cex.main =1.8,col=\"red\",main=\"$prefix[$i]\",xlab=\"Read depth(log10)\",ylab=\"Density\");\n";
	print R "lines(density(log10(A\$V14)),col=\"green\");\n";
	print R "legend(\"topright\",cex=1.5,lwd=2,legend=c(\"ReadDepth\",\"AverageDepth\"),col=c(\"red\",\"green\"),lty=c(1,1),bty=\"n\");\n";
	print R "plot(density(A\$V16/(A\$V15+A\$V16)),xlim=c(-0.2,1.2),type=\"l\",lwd=2,cex.axis=1.5,cex.lab =1.5,cex.main =1.8,col=\"red\",main=\"$prefix[$i]\",xlab=\"Read proportion\",ylab=\"Density\");\n";
	print R "lines(density(A\$V15/(A\$V15+A\$V16)),col=\"green\");\n";
	print R "legend(\"topright\",cex=1.5,lwd=2,legend=c(\"Alt\",\"Ref\"),col=c(\"red\",\"green\"),lty=c(1,1),bty=\"n\");\n";
	print R "plot(density(log10(A\$V17)),type=\"l\",lwd=2,cex.axis=1.5,cex.lab =1.5,cex.main =1.8,col=\"red\",main=\"$prefix[$i]\",xlab=\"Distance to adjacent allele(log10)\",ylab=\"Density\");\n";
	print R "plot(density(A\$V18),type=\"l\",lwd=2,cex.axis=1.5,cex.lab =1.5,cex.main =1.8,col=\"red\",main=\"$prefix[$i]\",xlab=\"Log10 likelihood ratio\",ylab=\"Density\");\n";
	print R "plot(density(log10(A\$V19)),type=\"l\",lwd=2,cex.axis=1.5,cex.lab =1.5,cex.main =1.8,col=\"red\",main=\"$prefix[$i]\",xlab=\"Pvalue(Fisher exect test,log10)\",ylab=\"Density\");\n";
	print R "plot(density(A\$V21),xlim=c(0,100),type=\"l\",lwd=2,cex.axis=1.5,cex.lab =1.5,cex.main =1.8,col=\"green\",main=\"$prefix[$i]\",xlab=\"GC content\",ylab=\"Density\");\n";
	print R "lines(density(A\$V20),col=\"red\");\n";
	print R "legend(\"topright\",cex=1.5,lwd=2,legend=c(\"RNAedit\",\"Gene\"),col=c(\"red\",\"green\"),lty=c(1,1),bty=\"n\");\n";
	print R "boxplot(A\$V11 ~ A\$V8,las=2,col=rainbow($num),lwd=2,cex.axis=1.5,cex.lab =1.5,cex.main =1.8,main=\"$prefix[$i]\",ylab=\"Alt allele proportion\");\n";
	print R "boxplot(A\$V12 ~ A\$V8,las=2,col=rainbow($num),lwd=2,cex.axis=1.5,cex.lab =1.5,cex.main =1.8,main=\"$prefix[$i]\",ylab=\"Error proportion\");\n";
	print R "boxplot(log10(A\$V13) ~ A\$V8,las=2,col=rainbow($num),lwd=2,cex.axis=1.5,cex.lab =1.5,cex.main =1.8,main=\"$prefix[$i]\",ylab=\"Total read depth(log10)\");\n";
	print R "boxplot(log10(A\$V16) ~ A\$V8,las=2,col=rainbow($num),lwd=2,cex.axis=1.5,cex.lab =1.5,cex.main =1.8,main=\"$prefix[$i]\",ylab=\"Alt read depth(log10)\");\n";
	print R "boxplot(log2(A\$V13/A\$V14) ~ A\$V8,las=2,col=rainbow($num),lwd=2,cex.axis=1.5,cex.lab =1.5,cex.main =1.8,main=\"$prefix[$i]\",ylab=\"Read depth/Average depth(log2)\");\n";
	print R "boxplot(log10(A\$V17) ~ A\$V8,las=2,col=rainbow($num),lwd=2,cex.axis=1.5,cex.lab =1.5,cex.main =1.8,main=\"$prefix[$i]\",ylab=\"Distance to adjacent allele(log10)\");\n";
	print R "boxplot(A\$V18 ~ A\$V8,las=2,col=rainbow($num),lwd=2,cex.axis=1.5,cex.lab =1.5,cex.main =1.8,main=\"$prefix[$i]\",ylab=\"Log10 likelihood ratio\");\n";
	print R "boxplot(log10(A\$V19) ~ A\$V8,las=2,col=rainbow($num),lwd=2,cex.axis=1.5,cex.lab =1.5,cex.main =1.8,main=\"$prefix[$i]\",ylab=\"Pvalue(Fisher exect test,log10)\");\n";
	print R "boxplot(A\$V20 ~ A\$V8,las=2,col=rainbow($num),lwd=2,cex.axis=1.5,cex.lab =1.5,cex.main =1.8,main=\"$prefix[$i]\",ylab=\"GC content\");\n";
	print R "dev.off();\n";
	close R;
	system "R CMD BATCH $prefix\_$prefix[$i].R";
}

if (@prefix>1) {
	open (DPR,">$prefix\_DPR.out")||die("fail to open $prefix\_DPR.out.\n");
	print DPR "Chr\tPosition\tProportion1\tAlt1\tRef1\tProportion2\tAlt2\tRef2\tFold(log2)\tPvalue\n";
	my %cluster=();
	for (my $i=0;$i<@prefix-1;$i++) {
		for (my $j=1;$j<@prefix;$j++) {
			print DPR "\#$prefix[$i] vs $prefix[$j]\n";
			my %sample=();
			my %sample1=();
			my %sample2=();
			open (IN1,"$prefix\_$prefix[$i].out")||die("fail to open $prefix\_$prefix[$i].out.\n");
			while (<IN1>) {
				chomp;
				#Chr	Accession	Name	Genome_pos	Gene_pos	AA_pos	Phase	Ref->Alt	RefCodon->AltCodon	RefAA->AltAA	AltRatio	ErrorRatio	Depth	AverageDepth	RefRead	AltRead	AdjacentAlt	LLR	Pvalue(fisher)	RNAeditGC	GeneGC	Exon_num	Exon_start	Exon_end	Function
				#cp	4	rpl23	2143	71	24	2	C->T	TCT->TTT	S->F	0.920792079207921	0	202	197.2	16	186	18	1000	4.31666359754032e-97	25.00	37.59	1	2213,	1932,	ribosomal protein L23
				next if (/^\#/);
				my @list=split /\t/,$_;
				$sample1{$list[0]}{$list[3]}="$list[10]\t$list[14]\t$list[15]";
				$sample{$list[0]}{$list[3]}++;
			}
			close IN1;
			open (IN2,"$prefix\_$prefix[$j].out")||die("fail to open $prefix\_$prefix[$j].out.\n");
			while (<IN2>) {
				chomp;
				#Chr	Accession	Name	Genome_pos	Gene_pos	AA_pos	Phase	Ref->Alt	RefCodon->AltCodon	RefAA->AltAA	AltRatio	ErrorRatio	Depth	AverageDepth	RefRead	AltRead	AdjacentAlt	LLR	Pvalue(fisher)	RNAeditGC	GeneGC	Exon_num	Exon_start	Exon_end	Function
				#cp	4	rpl23	2143	71	24	2	C->T	TCT->TTT	S->F	0.920792079207921	0	202	197.2	16	186	18	1000	4.31666359754032e-97	25.00	37.59	1	2213,	1932,	ribosomal protein L23
				next if (/^\#/);
				my @list=split /\t/,$_;
				$sample2{$list[0]}{$list[3]}="$list[10]\t$list[14]\t$list[15]";
				$sample{$list[0]}{$list[3]}++;
			}
			close IN2;
			#fisher exact test
			foreach my $chr (keys %sample) {
				foreach my $pos (keys %{$sample{$chr}}) {
					my $proportion1=0;my $proportion2=0;my $ref1=0;my $ref2=0;my $alt1=0;my $alt2=0;
					($proportion1,$ref1,$alt1)=split /\t/,$sample1{$chr}{$pos} if (exists $sample1{$chr}{$pos});
					($proportion2,$ref2,$alt2)=split /\t/,$sample2{$chr}{$pos} if (exists $sample2{$chr}{$pos});
					$cluster{$chr}{$pos}{$i}=$proportion1;
					$cluster{$chr}{$pos}{$j}=$proportion2;
					my $pfisher=&fisher_twotailed($alt1,$ref1,$alt2,$ref2);
					my $fold="NA";
					if ($proportion1>0 and $proportion2>0) {
						$fold=sprintf ("%.2f",&log2($proportion2/$proportion1));
					}
					print DPR "$chr\t$pos\t$proportion1\t$alt1\t$ref1\t$proportion2\t$alt2\t$ref2\t$fold\t$pfisher\n";
				}
			}
			print DPR "\n\n";
		}
	}
	close DPR;

	open (OUT,">$prefix\_cluster.out")||die("fail to open $prefix\_cluster.out.\n");
	for (my $k=0;$k<@prefix;$k++) {
		print OUT "\t$prefix[$k]";
	}
	print OUT "\n";
	foreach my $chr (keys %cluster) {
		foreach my $pos (keys %{$cluster{$chr}}) {
			print OUT "$chr\-$pos";
			for (my $i=0;$i<@prefix;$i++) {
				print OUT "\t$cluster{$chr}{$pos}{$i}";
			}
			print OUT "\n";
		}
	}
	close OUT;

	open (R,">$prefix\_cluster.R")||die("fail to open $prefix\_cluster.R.\n");
	print R "jpeg(\"$prefix\_cluster.jpeg\",width=6000,height=9000,res=600);\n";
	print R "A<-read.table(\"$prefix\_cluster.out\",sep=\"\\t\",header=TRUE,row.names=1);\n";
	print R "require(graphics);\n";
	print R "require(grDevices);\n";
	print R "x  <- as.matrix(A);\n";
	print R "rc <- rainbow(nrow(x), start = 0, end = .3);\n";
	print R "cc <- rainbow(ncol(x), start = 0, end = .3);\n";
	print R "hv <- heatmap(x, col = cm.colors(256), scale = \"column\",RowSideColors = rc, ColSideColors = cc, margins = c(5,10),xlab = \"Samples\", ylab =  \"RNA editing\",main = \"heatmap(RNA editing)\");\n";
	print R "dev.off();\n";
	close R;
	system "R CMD BATCH $prefix\_cluster.R";
}

#&RNA_alt($genome_pos,$gene_pos,$strand,$str);
sub RNA_alt {
	my $genome_pos=shift;
	my $gene_pos=shift;
	my $strand=shift;
	my $str=shift;
	my $gene_gc=shift;
	my $genetype=shift;

	my $num=0;
	my $substr="";
	my $printstr="";
	if ($genetype eq "CDS") {
		my $pos_aa=int(($gene_pos-1)/3)+1;
		my $phase;
		if ($gene_pos%3==1) {
			$phase=1;
		}elsif ($gene_pos%3==2) {
			$phase=2;
		}elsif ($gene_pos%3==0) {
			$phase=3;
		}
		for (my $k=0;$k<@vcffile;$k++) {
			my $filehandle="OUT$k";
			if (exists $vcf{$k}{$id}{$genome_pos}) {
				my ($ref_base,$alt_base,$ratio,$err,$depth,$dep,$ref,$alt,$dis,$llr,$pfisher)=split /\t/,$vcf{$k}{$id}{$genome_pos};
				if ($strand eq "-") {
					$ref_base=~tr/AaTtGgCc/TtAaCcGg/;
					$alt_base=~tr/AaTtGgCc/TtAaCcGg/;
				}
				if (length($ref_base)==length($alt_base) and length($ref_base)==1) {
					my $ref_codon=substr($str,$gene_pos-$phase,3);
					my @alt_codon=split //,$ref_codon;
					my $alt_codon=undef;
					for (my $l=0;$l<@alt_codon;$l++) {
						if ($l==$phase-1) {
							$alt_codon.=$alt_base;
						}else {
							$alt_codon.=$alt_codon[$l];
						}
					}
					my $gcstr="";
					if ($gene_pos-$window>=0) {
						$gcstr=substr($str,$gene_pos-$window,$window*2);
					}else {
						$gcstr=substr($str,0,$window*2);
					}
					my $gc=&gc($gcstr);
					$hydro{$k}{"$ref_base\->$alt_base"}{"$phase"}++;
					my $ref_aa=&codon2amino($ref_codon);
					my $alt_aa=&codon2amino($alt_codon);
					my $ref_hydro=&hydrophobic($ref_aa);
					my $alt_hydro=&hydrophobic($alt_aa);
					my $type;
					if ($ref_base eq "C" and $alt_base eq "T") {
						$type="+";
					}else {
						$type="-";
					}
					if ($ref_aa eq $alt_aa) {
						if ($type eq "+") {
							if ($ratio>=0.5) {
							}else {
								my $flag=&phase_phobic("-",$phase,$ref_hydro,$alt_hydro,$type);
								next if ($flag eq "-");
							}
						}else {
							my $flag=&phase_phobic("-",$phase,$ref_hydro,$alt_hydro,$type);
							next if ($flag eq "-");
						}
						$hydro{$k}{"$ref_base\->$alt_base"}{"silent"}++;
					}elsif ($ref_aa ne $alt_aa) {
						if ($type eq "+") {
							if ($ratio>=0.5) {
							}else {
								my $flag=&phase_phobic("+",$phase,$ref_hydro,$alt_hydro,$type);
								next if ($flag eq "-");
							}
						}else {
							my $flag=&phase_phobic("+",$phase,$ref_hydro,$alt_hydro,$type);
							next if ($flag eq "-");
						}
						$hydro{$k}{"$ref_base\->$alt_base"}{"non_silent"}++;
						if ($ref_hydro eq $alt_hydro) {
							if ($ref_hydro eq "+") {
								$hydro{$k}{"$ref_base\->$alt_base"}{"phobicphobic"}++;
							}else {
								$hydro{$k}{"$ref_base\->$alt_base"}{"philicphilic"}++;
							}
						}else {
							if ($ref_hydro eq "+") {
								$hydro{$k}{"$ref_base\->$alt_base"}{"phobicphilic"}++;
							}else {
								$hydro{$k}{"$ref_base\->$alt_base"}{"philicphobic"}++;
							}
						}
					}
					$substr.="$ref_base->$alt_base:$ref_codon->$alt_codon:$ref_aa->$alt_aa:$ratio:$err:$depth:$dep:$ref:$alt:$dis:$llr:$pfisher:$gc:$gene_gc\t";
					$num++;
					$change{$k}{"$ref_base\->$alt_base"}++;
					print $filehandle "$id\t$strand\t$name\t$genome_pos\t$gene_pos\t$pos_aa\t$phase\t$ref_base\->$alt_base\t$ref_codon\->$alt_codon\t$ref_aa->$alt_aa\t$ratio\t$err\t$depth\t$dep\t$ref\t$alt\t$dis\t$llr\t$pfisher\t$gc\t$gene_gc\t$exon_num\t$exon_start\t$exon_end\t$product\n";
				}elsif ($indel ne "n") {
					my $gcstr="";
					if ($gene_pos-$window>=0) {
						$gcstr=substr($str,$gene_pos-$window,$window*2);
					}else {
						$gcstr=substr($str,0,$window*2);
					}
					my $gc=&gc($gcstr);
					if ($ratio>=0.5) {
						$substr.="$ref_base->$alt_base:?->?:?->?:$ratio:$err:$depth:$dep:$ref:$alt:$dis:$llr:$pfisher:$gc:$gene_gc\t";
						$num++;
						$change{$k}{""}++;
						print $filehandle "$id\t$strand\t$name\t$genome_pos\t$gene_pos\t$pos_aa\t$phase\t$ref_base\->$alt_base\t?\->?\t?->?\t$ratio\t$err\t$depth\t$dep\t$ref\t$alt\t$dis\t$llr\t$pfisher\t$gc\t$gene_gc\t$exon_num\t$exon_start\t$exon_end\t$product\n";
					}
				}
			}else {
				$substr.="\t";
			}
		}
		if ($num>0) {
			$printstr="$id\t$strand\t$name\t$genome_pos\t$gene_pos\t$pos_aa\t$num\t$phase\t$substr$exon_num\t$exon_start\t$exon_end\t$product\n";
		}
	}else {
		for (my $k=0;$k<@vcffile;$k++) {
			my $filehandle="OUT$k";
			if (exists $vcf{$k}{$id}{$genome_pos}) {
				my ($ref_base,$alt_base,$ratio,$err,$depth,$dep,$ref,$alt,$dis,$llr,$pfisher)=split /\t/,$vcf{$k}{$id}{$genome_pos};
				if ($strand eq "-") {
					$ref_base=~tr/AaTtGgCc/TtAaCcGg/;
					$alt_base=~tr/AaTtGgCc/TtAaCcGg/;
				}
				if (length($ref_base)==length($alt_base) and length($ref_base)==1) {
					my $gcstr="";
					if ($gene_pos-$window>=0) {
						$gcstr=substr($str,$gene_pos-$window,$window*2);
					}else {
						$gcstr=substr($str,0,$window*2);
					}
					my $gc=&gc($gcstr);
					if ($ratio>=0.5) {
						if ($ref_base eq "C" and $alt_base eq "T") {
							$substr.="$ref_base->$alt_base:?->?:?->?:$ratio:$err:$depth:$dep:$ref:$alt:$dis:$llr:$pfisher:$gc:$gene_gc\t";
							$num++;
							$change{$k}{""}++;
							print $filehandle "$id\t$strand\t$name\t$genome_pos\t$gene_pos\t?\t?\t$ref_base\->$alt_base\t?\->?\t?->?\t$ratio\t$err\t$depth\t$dep\t$ref\t$alt\t$dis\t$llr\t$pfisher\t$gc\t$gene_gc\t$exon_num\t$exon_start\t$exon_end\t$type\n";
						}
					}
				}elsif ($indel ne "n") {
					my $gcstr="";
					if ($gene_pos-$window>=0) {
						$gcstr=substr($str,$gene_pos-$window,$window*2);
					}else {
						$gcstr=substr($str,0,$window*2);
					}
					my $gc=&gc($gcstr);
					if ($ratio>=0.5) {
						$substr.="$ref_base->$alt_base:?->?:?->?:$ratio:$err:$depth:$dep:$ref:$alt:$dis:$llr:$pfisher:$gc:$gene_gc\t";
						$num++;
						$change{$k}{""}++;
						print $filehandle "$id\t$strand\t$name\t$genome_pos\t$gene_pos\t?\t?\t$ref_base\->$alt_base\t?\->?\t?->?\t$ratio\t$err\t$depth\t$dep\t$ref\t$alt\t$dis\t$llr\t$pfisher\t$gc\t$gene_gc\t$exon_num\t$exon_start\t$exon_end\t$type\n";
					}
				}
			}else {
				$substr.="\t";
			}
		}
		if ($num>0) {
			$printstr="$id\t$strand\t$name\t$genome_pos\t$gene_pos\t?\t$num\t?\t$substr$exon_num\t$exon_start\t$exon_end\t$type\n";
		}
	}
	return $printstr;
}

sub chiX2 {
	my $n11=shift;
	my $n12=shift;
	my $n21=shift;
	my $n22=shift;

	my $np1=$n11+$n21;
	my $np2=$n12+$n22;
	my $n1p=$n11+$n12;
	my $n2p=$n21+$n22;
	my $npp=$np1+$np2;
	my $num1=0;
	my $num2=0;
	my $num3=0;
	my $num4=0;
	$num1=($n11-$n1p/$npp*$np1)**2/($n1p/$npp*$np1) if ($npp*$np1>0);
	$num2=($n12-$n1p/$npp*$np2)**2/($n1p/$npp*$np2) if ($npp*$np2>0);
	$num3=($n21-$n2p/$npp*$np1)**2/($n2p/$npp*$np1) if ($npp*$np1>0);
	$num4=($n22-$n2p/$npp*$np2)**2/($n2p/$npp*$np2) if ($npp*$np2>0);
	my $chisquare=($num1+$num2+$num3+$num4); 
	return (sprintf "%.2f",$chisquare);
}

sub pvalue {
	my $n = shift;
	$n=sprintf "%.1f",$n;

	my $p=0.001;
	if (exists $pvalue{$n}) {
		$p=$pvalue{$n};
	}
	return $p;
}

my %pvalue = (
	"0.0" => "1.000",
	"0.1" => "0.752",
	"0.2" => "0.655",
	"0.3" => "0.584",
	"0.4" => "0.527",
	"0.5" => "0.480",
	"0.6" => "0.439",
	"0.7" => "0.403",
	"0.8" => "0.371",
	"0.9" => "0.343",
	"1.0" => "0.317",
	"1.1" => "0.294",
	"1.2" => "0.273",
	"1.3" => "0.254",
	"1.4" => "0.237",
	"1.5" => "0.221",
	"1.6" => "0.206",
	"1.7" => "0.192",
	"1.8" => "0.180",
	"1.9" => "0.168",
	"2.0" => "0.157",
	"2.1" => "0.147",
	"2.2" => "0.138",
	"2.3" => "0.129",
	"2.4" => "0.121",
	"2.5" => "0.114",
	"2.6" => "0.107",
	"2.7" => "0.100",
	"2.8" => "0.094",
	"2.9" => "0.089",
	"3.0" => "0.083",
	"3.1" => "0.078",
	"3.2" => "0.074",
	"3.3" => "0.069",
	"3.4" => "0.065",
	"3.5" => "0.061",
	"3.6" => "0.058",
	"3.7" => "0.054",
	"3.8" => "0.051",
	"3.9" => "0.048",
	"4.0" => "0.046",
	"4.1" => "0.043",
	"4.2" => "0.040",
	"4.3" => "0.038",
	"4.4" => "0.036",
	"4.5" => "0.034",
	"4.6" => "0.032",
	"4.7" => "0.030",
	"4.8" => "0.028",
	"4.9" => "0.027",
	"5.0" => "0.025",
	"5.1" => "0.024",
	"5.2" => "0.023",
	"5.3" => "0.021",
	"5.4" => "0.020",
	"5.5" => "0.019",
	"5.6" => "0.018",
	"5.7" => "0.017",
	"5.8" => "0.016",
	"5.9" => "0.015",
	"6.0" => "0.014",
	"6.1" => "0.014",
	"6.2" => "0.013",
	"6.3" => "0.012",
	"6.4" => "0.011",
	"6.5" => "0.011",
	"6.6" => "0.010",
	"6.7" => "0.010",
	"6.8" => "0.009",
	"6.9" => "0.009",
	"7.0" => "0.008",
	"7.1" => "0.008",
	"7.2" => "0.007",
	"7.3" => "0.007",
	"7.4" => "0.007",
	"7.5" => "0.006",
	"7.6" => "0.006",
	"7.7" => "0.006",
	"7.8" => "0.005",
	"7.9" => "0.005",
	"8.0" => "0.005",
	"8.1" => "0.004",
	"8.2" => "0.004",
	"8.3" => "0.004",
	"8.4" => "0.004",
	"8.5" => "0.004",
	"8.6" => "0.003",
	"8.7" => "0.003",
	"8.8" => "0.003",
	"8.9" => "0.003",
	"9.0" => "0.003",
	"9.1" => "0.003",
	"9.2" => "0.002",
	"9.3" => "0.002",
	"9.4" => "0.002",
	"9.5" => "0.002",
	"9.6" => "0.002",
	"9.7" => "0.002",
	"9.8" => "0.002",
	"9.9" => "0.002",
	"10.0" => "0.002",
	"10.1" => "0.001",
	"10.2" => "0.001",
	"10.3" => "0.001",
	"10.4" => "0.001",
	"10.5" => "0.001",
	"10.6" => "0.001",
	"10.7" => "0.001",
	"10.8" => "0.001",
	"10.9" => "0.001",
	"11.0" => "0.001",
	"11.1" => "0.001",
	"11.2" => "0.001",
	"11.3" => "0.001",
	"11.4" => "0.001",
	"11.5" => "0.001",
	"11.6" => "0.001",
	"11.7" => "0.001",
	"11.8" => "0.001",
	"11.9" => "0.001",
	"12.0" => "0.001",
	"12.1" => "0.001",
	"12.2" => "0.000",	
);

sub codon2amino {
	my $codon=$_[0];
	$codon = uc($codon);
	my(%genetic_code)=(
	'TCA' => 'S', # Serine ££Ë¿°±Ëá 
	'TCC' => 'S', # Serine ££Ë¿°±Ëá 
	'TCG' => 'S', # Serine ££Ë¿°±Ëá 
	'TCT' => 'S', # Serine ££Ë¿°±Ëá 
	'TTC' => 'F', # Phenylalanine ££±½±û°±Ëá 
	'TTT' => 'F', # Phenylalanine ££±½±û°±Ëá 
	'TTA' => 'L', # Leucine ££ÁÁ°±Ëá 
	'TTG' => 'L', # Leucine ££ÁÁ°±Ëá 
	'TAC' => 'Y', # Tyrosine ££ÀÒ°±Ëá 
	'TAT' => 'Y', # Tyrosine ££ÀÒ°±Ëá 
	'TAA' => '*', # Stop ££Í£Ö¹ 
	'TAG' => '*', # Stop ££Í£Ö¹ 
	'TGC' => 'C', # Cysteine ££°ëë×°±Ëá 
	'TGT' => 'C', # Cysteine ££°ëë×°±Ëá 
	'TGA' => '*', # Stop ££Í£Ö¹ 
	'TGG' => 'W', # Tryptophan ££É«°±Ëá 
	'CTA' => 'L', # Leucine ££ÁÁ°±Ëá 
	'CTC' => 'L', # Leucine ££ÁÁ°±Ëá 
	'CTG' => 'L', # Leucine ££ÁÁ°±Ëá 
	'CTT' => 'L', # Leucine ££ÁÁ°±Ëá 
	'CCA' => 'P', # Proline ££¸¬°±Ëá 
	'CAT' => 'H', # Histidine ££×é°±Ëá 
	'CAA' => 'Q', # Glutamine ££¹È°±õ£°· 
	'CAG' => 'Q', # Glutamine ££¹È°±õ£°· 
	'CGA' => 'R', # Arginine ££¾«°±Ëá 
	'CGC' => 'R', # Arginine ££¾«°±Ëá 
	'CGG' => 'R', # Arginine ££¾«°±Ëá 
	'CGT' => 'R', # Arginine ££¾«°±Ëá 
	'ATA' => 'I', # Isoleucine ££ÒìÁÁ°±Ëá 
	'ATC' => 'I', # Isoleucine ££ÒìÁÁ°±Ëá 
	'ATT' => 'I', # Isoleucine ££ÒìÁÁ°±Ëá 
	'ATG' => 'M', # Methionine ££µ°°±Ëá 
	'ACA' => 'T', # Threonine ££ËÕ°±Ëá 
	'ACC' => 'T', # Threonine ££ËÕ°±Ëá 
	'ACG' => 'T', # Threonine ££ËÕ°±Ëá 
	'ACT' => 'T', # Threonine ££ËÕ°±Ëá 
	'AAC' => 'N', # Asparagine ££Ìì¶¬õ£°· 
	'AAT' => 'N', # Asparagine ££Ìì¶¬õ£°· 
	'AAA' => 'K', # Lysine ££Àµ°±Ëá 
	'AAG' => 'K', # Lysine ££Àµ°±Ëá 
	'AGC' => 'S', # Serine ££Ë¿°±Ëá 
	'AGT' => 'S', # Serine ££Ë¿°±Ëá 
	'AGA' => 'R', # Arginine ££¾«°±Ëá 
	'AGG' => 'R', # Arginine ££¾«°±Ëá 
	'CCC' => 'P', # Proline ££¸¬°±Ëá 
	'CCG' => 'P', # Proline ££¸¬°±Ëá 
	'CCT' => 'P', # Proline ££¸¬°±Ëá 
	'CAC' => 'H', # Histidine ££×é°±Ëá 
	'GTA' => 'V', # Valine ££çÓ°±Ëá 
	'GTC' => 'V', # Valine ££çÓ°±Ëá 
	'GTG' => 'V', # Valine ££çÓ°±Ëá 
	'GTT' => 'V', # Valine ££çÓ°±Ëá 
	'GCA' => 'A', # Alanine ££±û°±Ëá 
	'GCC' => 'A', # Alanine ££±û°±Ëá 
	'GCG' => 'A', # Alanine ££±û°±Ëá 
	'GCT' => 'A', # Alanine ££±û°±Ëá 
	'GAC' => 'D', # Aspartic Acid ££Ìì¶¬°±Ëá 
	'GAT' => 'D', # Aspartic Acid ££Ìì¶¬°±Ëá 
	'GAA' => 'E', # Glutamic Acid ££¹È°±Ëá 
	'GAG' => 'E', # Glutamic Acid ££¹È°±Ëá 
	'GGA' => 'G', # Glycine ££¸Ê°±Ëá 
	'GGC' => 'G', # Glycine ££¸Ê°±Ëá 
	'GGG' => 'G', # Glycine ££¸Ê°±Ëá 
	'GGT' => 'G' # Glycine ££¸Ê°±Ëá 
	);
	if (defined $genetic_code{$codon}) {
		return $genetic_code{$codon};
	}else {
		return "";
	}
}

sub fisher_left {
	my $n11=shift;
	my $n12=shift;
	my $n21=shift;
	my $n22=shift;
	my $np1=$n11+$n21;
	my $np2=$n12+$n22;
	my $n1p=$n11+$n12;
	my $n2p=$n21+$n22;
	my $npp=$np1+$np2;
	my $left_value = Text::NSP::Measures::2D::Fisher::left::calculateStatistic( n11=>$n11,n1p=>$n1p,np1=>$np1,npp=>$npp);
	return $left_value;
}

sub fisher_twotailed {
	my $n11=shift;
	my $n12=shift;
	my $n21=shift;
	my $n22=shift;
	my $np1=$n11+$n21;
	my $np2=$n12+$n22;
	my $n1p=$n11+$n12;
	my $n2p=$n21+$n22;
	my $npp=$np1+$np2;
	my $left_value = Text::NSP::Measures::2D::Fisher::twotailed::calculateStatistic( n11=>$n11,n1p=>$n1p,np1=>$np1,npp=>$npp);
	return $left_value;
}

sub log10 {
	my $n = shift;
	return log($n)/log(10);
}

sub log2 {
	my $n = shift;
	return log($n)/log(2);
}

sub hydrophobic {
	my $aa=shift;
	my(%hydropathy)=(
	'G' => 'G',
	'A' => 'A',
	'V' => 'V',
	'L' => 'L',
	'I' => 'I',
	'F' => 'F',
	'W' => 'W',
	'M' => 'M',
	'P' => 'P',
	);
	if (exists $hydropathy{$aa}) {
		return "+";
	}else {
		return "-";
	}
}

sub gc {
	 my $seq = $_[0];
	 my @seqarray = split('',$seq);
	 my $count = 0;
	 foreach my $base (@seqarray) {
	   $count++ if $base =~ /[G|C]/;
	 }
	 my $len = $#seqarray+1;
	 my $num=$count/$len*100;
	 my $dec=sprintf "%.2f",$num;
	 return $dec;
}

#Hydrophobic G A V L I F W M P
#Hydrophilic D H N E K Q R S T C Y

#my $flag=&phase_phobic($change,$phase,$ref_hydro,$alt_hydro);
sub phase_phobic {
	my $change = $_[0];
	my $phase = $_[1];
	my $refhydro = $_[2];
	my $althydro = $_[3];
	my $flag = $_[4];
	my $type;

	if ($change eq "-") {
		if ($phase==2) {
			$type="+";
		}else {
			if ($althydro eq "+") {
				$type="+";
			}else {
				if ($flag eq "+" and $phase==1) {
					$type="+";
				}else {
					$type="-";
				}
			}
		}
	 }else {
		 if ($phase==1) {
			 if ($refhydro eq "-" and $althydro eq "+") {
				 $type="+";
			 }elsif ($refhydro eq "+" and $althydro eq "+") {
				 $type="+";
			 }else {
				 if ($flag eq "+") {
					  $type="+";
				 }else {
					$type="-";
				 }
			 }
		 }elsif ($phase==2) {
			 $type="+";
		 }else {
			 if ($refhydro eq "-" and $althydro eq "+") {
				 $type="+";
			 }elsif ($refhydro eq "+" and $althydro eq "+") {
				 if ($flag eq "+") {
					  $type="+";
				 }else {
					 $type="-";
				 }
			 }elsif ($refhydro eq "-" and $althydro eq "-") {
				  if ($flag eq "+") {
					  $type="+";
				 }else {
					 $type="-";
				 }
			 }else {
				 $type="-";
			 }
		 }
	 }
	 return $type;
}
	 
