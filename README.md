# REDO
RNA Editing Detection in Organelle based on variant calling results
NAME
====

REDO - RNA Editing Detection in Organelle 

VERSION
=======

Version 1.0

COPYRIGHT
=========

Copyright 2017 - Shuangyang Wu & Wanfei Liu & Qiang Lin - Beijing Institute of Genomics, Chinese Academy of Sciences.

This tool is a free package; you can redistribute it and/or modify it freely.

INTRODUCTION
============

REDO is a comprehensive application tool for identifying RNA editing events in organelles based on variant call results (VCF files). It is a suite of Perl scripts and can work easily and directly in any operating system installed Perl and R Environment. The stringent rule depended filters and statistical filters are used in REDO for reducing false positive rate. It can provide detailed annotations, statistics and figures for RNA editing sites. REDO also can detect RNA editing events in multiple samples simultaneously and identify differential proportion of RNA editing events in different samples. Moreover, the genome variation can be easily removed by a subprogram fish.pl in our package.


REDO only uses the variant calling format (VCF) files (records for all sites), the genome sequence file (for annotation) and the gene annotation file (for annotation, feature table file, http://www.ncbi.nlm.nih.gov/projects/Sequin/table.html) as input, then the raw variants are filtered by rule-depended filters and statistics filters for reducing the false positive rate as following steps. 1. Quality control filter: the low quality sites are filtered according to the reads quality (DP=0 in GATK (reads with MQ=255 or with bad mates are filtered) and DP4=0 in samtools (number of high-quality ref-forward, ref-reverse, alt-forward and alt-reverse bases)). 2. Depth filter: three measures are used in this filter, which are total reads depth (<4), alt reads depth (<3), and the total reads depth variation between alt site and ad-jacent sites without variant in specific window (total reads depth<0.2*average reads depth of adjacent sites). 3. Alt proportion filter: two measures are used in this filter, including alt proportion (<0.1) and the reduced alt proportion after removing possible sequencing error proportion obtained from the adjacent sites without variant in specific window (alt proportion – error proportion<0.1). 4. Multiple alt filter: only the variant with one alt allele is kept for RNA editing detection. 5. Distance filter: the variant sites with short distance from each other (<3) are filtered due to the possible position obstruction for RNA editing. 6. Splice junction filter: variants within short splice anchor (<2) are removed. 7. Indel filter: the indel variants are removed in default (indels can be kept when assigning "-i yes" option). 8. Likelihood ratio test filter: a likelihood ratio (LLR) test is used for RNA editing sites (Chepelev, 2012; Sun, et al., 2016). LLR test is a probabilistic test incorporating error probability of bases (error probability is obtained using adjacent sites without variant in specific window) for detecting RNA editing sites (LLR<10). 9. Fisher exact test filter: we assessed the significance for a given RNA editing site (alt reads, ref reads) by comparing its expected levels (0, alt reads + ref reads) using the Fisher’s exact test and the p-value of fisher exact test is used for filter (p-value>0.01). 10. The complicated model filter: based on the statistics results for ex-periment validated RNA editing sites and the attributes of codon table, a complicated filter model was built according to five characteristics of RNA editing sites, which are RNA editing types (C->T, A->G, T->C, etc), alt proportion, amino acids change, codon phase and hydropho-bic/hydrophilic change. 

PREREQUISITES
=============

Before running, the following software or program are required:
-	Perl Environment (perl v5.18.4 or later (it can be lower), tested for Windows and Linux)
-	Perl module: Getopt::Long; Text::NSP
-	R Environment (tested on R 3.1.2 or later (it can be lower) for Windows and Linux)
-	R module: graphics; grDevices

REDO have been tested on Windows and Linux.

INSTALLATION
============

This package doesn't need install. You can run it using absolute path after decompressed it. If you want to use it without absolute path, you should add the absolute path of REDO directory in PATH of your bash profile.

COMMAND LINE
============
	
Run REDO.pl without any parameter, it will print below usage on the screen.

*************
*1.0 version*
*************
************************************************************************
        Usage: REDO.pl -g genome_sequence -v variant_file -t tbl_file -o prefix_of_outfile -d reads_depth -c minimum_coverage_of_alt_allele -p alt_proportion -w window_size -s minimum_splice_anchor -a minimum_alt_distance -l llr_value -f fisher_pvalue -dv depth_variant -i identify_indel.
           -g: the absolute path of genome sequence.
           -v: the absolute path of variant file/files (for example, "-v test1.vcf -v test2.vcf -v test3.vcf ......").
           -t: the absolute path of tbl file (feature table file, http://www.ncbi.nlm.nih.gov/projects/Sequin/table.html).
           -o: the prefix of output file.
           -d: minimum reads depth (4).
           -c: minimum coverage of alternative (alt) allele (3).
           -p: minimum alt proportion (0.1).
           -w: minimum window size for calculating error rate and average depth (10).
           -s: minimum splice anchor size (2).
           -a: minimum alt distance (3).
           -l: minimum likelihood ratio (LLR) (10).
           -f: maximum p-value of fisher exact test (0.01).
          -dv: the depth variant (0.2).
           -i: identify indel (y or n, default n).
************************************************************************

Note: all files should be uncompressed files. For VCF file, it should include all positions with reads supports (samtools/bcftools default VCF output file or GATK output file with “-ERC BP_RESOLUTION”). We can provide one VCF file or multiple VCF files simultaneously. tbl file is the standard feature table file in NCBI (http://www.ncbi.nlm.nih.gov/projects/Sequin/table.html), we can download it from each genome sequence record by click “Send -> Complete Record -> File -> Format -> Feature Table -> Create File” easily in nucleotide database of NCBI. If you want run REDO multiple times with different parameters, please use “-o” option to assign different prefix name for output files. For calculating the error rate and the average reads depth, a default 10bp size window is used. To filter false positive sites, we use 4 as the default minimum reads depth and 3 as the minimum alt allele reads coverage, because most of false positive variants identified by GATK and samtools come from the low depth and alt coverage records according to our previous study. RNA editing may need second structure or motif to recognize related enzymes and we didn’t find the RNA editing sites with short distance from each other in several experiment data except for Arabidopsis (<3bp). Therefore, we use 3 as the minimum alt distance and 2 as the minimum splice anchor size. According to the test using experiment data, the likelihood ratio (LLR) should be equal or larger than 10 (default value). Moreover, according to the statistics results for experiment validated RNA editing sites and the attributes of codon table, we applied a complicated filter model. The filter model based on five characteristics of RNA editing sites, which are RNA editing types (C->T, A->G, T->C, etc), alt proportion, amino acids change, codon phase and hydrophobic/hydrophilic change. Using this filter model, we lost 0%~1.34% sensitivity, but increased 0.96%~10.45% precision rate and 0%~0.13% specificity percent. Even though, our sensitivity is comparable with REDItools. Furthermore, the chromosome ID in genome sequence file should be identical with vcf and tbl file. At present, we only identify C->T RNA editing sites with at least 0.50 editing proportion in tRNA and rRNA region.

INPUTFILES
==========

REDO need three input files:
-	VCF file (Variant Call Format, the result of GATK and samtools). REDO can accept one or multiple VCF files by option "-v".
-	Reference genome in Fasta format. Chromosome/region names must be equal to those in VCF file and tble file.
-	tbl file (feature table file, http://www.ncbi.nlm.nih.gov/projects/Sequin/table.html). It can be downloaded from NCBI directly.

Note: the chromosome ID in genome sequence file should be identical with vcf and tbl file.

OUTPUTFILES
===========

REDO mainly produces output files in simple tab-delimited table files.

1.	*_sample.out
The RNA editing result for each sample.
Main fields are:
-	Chr: chromosome name;
-	Strand: strand information ("+" or "-");
-	Name: gene name;
-	Genome_pos: nucleotide position in chromosome (1-based);
-	Gene_pos: nucleotide position in gene (1-based);
-	AA_pos: amino acid position in gene (1-based);
-	Phase: codon phase of editing site (1-based);
-	Ref->Alt: reference allele to editing allele;
-	RefCodon->AltCodon: reference codon to editing codon;
-	RefAA->AltAA: reference amino acid to editing amino acid;
-	AltRatio: editing reads proportion;
-	ErrorRatio: sequencing error rate obtained from the adjacent sites without variant in specific window (default window size is 10);
-	Depth: total reads number;
-	AverageDepth: average total reads number of adjacent sites without variant in specific window (default window size is 10);
-	RefRead: reference allele reads number;
-	AltRead: editing allele reads number;
-	AdjacentAlt: distance to adjacent candidate editing site;
-	LLR: a likelihood ratio (LLR) is a probabilistic test incorporating error probability of bases (error probability is obtained using adjacent sites without variant in specific window) for detecting RNA editing sites;
-	 Pvalue(fisher): p-value is obtained by assessing the significance for a given RNA editing site (alt reads, ref reads) by comparing its expected levels (0, alt reads + ref reads) using the Fisher’s exact test;
-	RNAeditGC: GC content of RNA editing site in specific window (default 20bp);
-	GeneGC: gene GC content;
-	Exon_num: exon number of gene;
-	Exon_start: exon start/starts of gene (multiple exon starts are separated by comma) (1-based);	
-	Exon_end: exon end/ends of gene (multiple exon ends are separated by comma) (1-based);
-	Function: gene function.

2.	*.out
The integrated RNA editing result for all samples.
Main fields are:
-	Chr: chromosome name;
-	Strand: strand information;
-	Name: gene name;
-	Genome_pos: nucleotide position in chromosome (1-based);
-	Gene_pos: nucleotide position in gene (1-based);
-	AA_pos: amino acid position in gene (1-based);
-	Sample_num: total sample number of supporting RNA editing site;
-	Phase: codon phase of editing site (1-based);
-	Sample_name1: the detailed editing information in the sample, including Ref->Alt, RefCodon->AltCodon, RefAA->AltAA, AltRatio, ErrorRatio, Depth, AverageDepth, RefRead, AltRead, AdjacentAlt, LLR, Pvalue(fisher), RNAeditGC, and GeneGC separated by colon;
-	Sample_name2;
-	......
-	Exon_num: exon number of gene;
-	Exon_start: exon start/starts of gene (multiple exon starts are separated by comma) (1-based);	
-	Exon_end: exon end/ends of gene (multiple exon ends are separated by comma) (1-based);
-	Function: gene function.

3.	*.change_matrix
The statistics result file for RNA editing in each sample, including total number and percent for each RNA editing type.
Main fields are:
-	Prefix: prefix of sample name;
-	Column 2 to column 13 are 12 possible RNA editing types, the total RNA editing number and percent for each type are shown.

4.	*.stat
The statistics result file for RNA editing in each sample, including total editing number, indel editing number, editing number in different codon phases, silent number, non-silent number, and the detailed statistics for each RNA editing type.  
Main fields are:
-	Prefix: prefix of sample name;
-	Total: total RNA editing number;
-	Other: non-substitution (indel) RNA editing number;
-	Phase(1,2,3): total RNA editing number in phase 1, 2 and 3 of codons for substitution RNA editing;
-	Silent: silent RNA editing number for substitution RNA editing;
-	NonSilent: non-silent RNA editing number for substitution RNA editing;
-	A->C(silent): silent RNA editing number for A->C RNA editing;
-	A->C(non_silent): non-silent RNA editing number for A->C RNA editing;
-	A->C(1:hydrophobic2hydrophobic,2:hydrophilic2hydrophilic,3:hydrophobic2hydrophilic,4:hydrophilic2hydrophobic): the number of hydrophobic/hydrophilic changes due to A->C RNA editing in four types separated by comma; 
-	Column 10 to column 42 are RNA editing number in silent, non-silent and 4 types of hydrophobic/hydrophilic changes for other 11 RNA editing types.

5.	*_sample.R and *_sample.jpeg
*_sample.R is the R script for drawing the attributes of RNA editing in each sample while *_sample.jpeg is the figure produced by *_sample.R, which shown 16 subfigures to illustrate the RNA editing attributes.

6.	*_DPR.out
The pairwise sample comparison for differential proportion of RNA editing events.
Main fields are:
-	#sample name1 vs sample name2: the sample name of comparison;
-	Chr: chromosome name;
-	Position: nucleotide position in chromosome (1-based);
-	Proportion1: editing allele proportion in sample one;
-	Alt1: editing allele reads in sample one; 
-	Ref1: reference allele reads in sample one;
-	Proportion2: editing allele proportion in sample two;
-	Alt2: editing allele reads in sample two; 
-	Ref2: reference allele reads in sample two;
-	Fold(log2): log2 fold change comparing sample two to one;
-	Pvalue: p-value is obtained by assessing the significance for a RNA editing site in sample one (alt1 reads, ref1 reads) by comparing to sample two (alt2 reads, ref2 reads) using the Fisher’s exact test.

7.	*_cluster.out, *_cluster.R and *_cluster.jpeg
*_cluster.out is a matrix for editing allele proportion of each RNA editing site in all samples. The columns are samples and the rows are RNA editing sites. *_cluster.R is the R script for doing the cluster analysis and *_cluster.jpeg is the cluster heatmap for RNA editing sites.

TEST & EXAMPLE
==============

REDO has been tested on coconut cp genome (KX028884) and related RNA-seq data (SRR1063404, SRR1063407, SRR1137438, SRR1173229, SRR1265939, SRR1273070, SRR1273180 and SRR606452, for comparing with REDItools) and Arabidopsis mt (NC_001284) and related RNA-seq data (SRR2079784, for comparing with REDItools). The cp RNA-seq data were mapped by GSNAP while mt RNA-seq data was mapped by BWA. All variants were identified by samtools. The small organelle genome, small VCF file and low memory request can make REDO run very fast in almost any computer (less than 2 minutes for 8 samples in coconut cp genome and less than 1 minute for 1 sample in Arabidopsis mt).

When you get the REDO package, decompress the directory, then you can test the program by following commands (just selected one command for each step).

1. Data pre-processing
1.1	Filter
java -jar trimmomatic-0.33.jar PE -threads 8 -trimlog SRR1063404.trimlog SRR1063404_1.fastq SRR1063404_2.fastq -baseout SRR1063404Filtered.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36&
1.2	Mapping
GSNAP:
gsnap -d coco_cp -a off -N 1 -A sam -t 8 --force-xs-dir --pairmax-rna=1000000 --split-output= SRR1063404_cp SRR1063404_1.paired.fq SRR1063404_2.paired.fq
BWA:
bwa mem -t 16 -M ath SRR2079784_1.paired.fastq.gz SRR2079784_2.paired.fastq.gz 1>SRR2079784.sam 2> SRR2079784.maperr &
1.3	Variant calling
Sam to bam: samtools view -b SRR2079784.sam -o SRR2079784.bam &
Sam sort: samtools sort -o SRR2079784.sorted.bam -O bam -T SRR2079784.temp SRR2079784.bam &
Sam mpileup: samtools mpileup -v -f arabidopsis_mt.fa SRR2079784.sorted.bam --output SRR2079784.sorted.vcf&
bcftools: nohup bcftools call -o SRR2079784.sorted.call.vcf -O v -c SRR2079784.sorted.vcf &

Note: we only used *.concordant_mult and *.concordant_uniq files as the final GSNAP mapping result while BWA default result was used as mapping result.

Table 1 The information of all test RNA-seq data
Species	Organelle	SRA accession No.	Length	Original fragments	High quality fragments	Percent	Mapping fragment	Percent
Coconut	cp	SRR1063404	202	36,009,632	32,555,041	90.41%	1,499,762	4.61%
		SRR1063407	202	35,467,948	32,141,745	90.62%	101,308	0.32%
		SRR1137438	152	50,839,994	42,267,444	83.14%	68,693	0.16%
		SRR1173229	202	119,333,177	113,394,045	95.02%	249,845	0.22%
		SRR1265939	202	51,540,183	48,892,847	94.86%	34,510	0.07%
		SRR1273070	337	40,564,276	37,752,443	93.07%	10,508	0.03%
		SRR1273180	252	60,030,680	54,291,251	90.44%	1,016,376	1.87%
		SRR606452	180	27,465,703	27,063,513	98.54%	548,364	2.03%
Arabidopsis	mt	SRR2079784	132	21,839,182	20,471,554	93.74%	1,796,028	8.77%

2. RNA editing identification:
2.1	REDItools
python REDItoolDenovo.py -i SRR2079784.sorted.bam -f arabidopsis_mt.fa -o SRR2079784 -c 4 >SRR2079784.nohup &
2.2	REDO
Coconut cp:
perl REDO.pl -v SRR1063404_cp.vcf -v SRR1063407_cp.vcf -v SRR1137438_cp.vcf -v SRR1173229_cp.vcf -v SRR1265939_cp.vcf -v SRR1273070_cp.vcf -v SRR1273180_cp.vcf -v SRR606452_cp.vcf -t cp.tbl -g cp.fsa -o coco_cp &
Arabidopsis mt:
perl REDO.pl -g mt.fa -v SRR2079784_mt.vcf -t mt.tbl -o ath_mt -a 0 -s 0&

Note: because Arabidopsis has a lot of RNA editing sites with short distance (adjacent RNA editing sites) or located in the junction points of transcripts in experiment data, we changed the default parameter -a and -s to 0 to recover these RNA editing sites.

Table 2 The information of RNA editing identification in test RNA-seq data
Species	Organelle	SRA accession No.	Samtools Variant	REDItools editing sites	Percent	REDO editing sites	Percent
Coconut	cp	SRR1063404	210	126	60.00%	96	45.71%
		SRR1063407	387	-	-	66	17.05%
		SRR1137438	197	-	-	51	25.89%
		SRR1173229	590	-	-	129	21.86%
		SRR1265939	477	-	-	49	10.27%
		SRR1273070	1426	375	26.30%	175	12.27%
		SRR1273180	229	-	-	88	38.43%
		SRR606452	247	-	-	68	27.53%
Arabidopsis	mt	SRR2079784	1219	553	45.37%	488	40.03%

3. RNA editing result:
3.1	Figures
	Figure 1 - The attributes of RNA editing sites in Arabidopsis mt.
	Figure 2 - The attributes of RNA editing sites in coconut SRR1063404.
	Figure 3 - The attributes of RNA editing sites in coconut SRR1273070.
	Figure 4 - The cluster of RNA editing sites in coconut cp based on 8 RNA-seq data according to the RNA editing proportions.
	Figure 5 – The attributes comparison between true positive and false positive RNA editing sites in mt (not produced by this program).

3.2	Tables
	Table 3 - ath_mt_SRR2079784_mt.out (partial).
	Table 4 - coco_cp.out (partial).
	Table 5 - coco_cp.change_matrix.
	Table 6 - coco_cp.stat.
	Table 7 - coco_cp_DPR.out (partial).
	Table 8 - The RNA editing sites identification between REDO and REDItools in three analyzed RNA-seq datasets (without applying the complicated filter model).
	Table 9 - The RNA editing sites identification between REDO and REDItools in three analyzed RNA-seq datasets (with applying the complicated filter model).

3.2.1	Table 3 - ath_mt_SRR2079784_mt.out (partial).
#Chr	Strand	Name	Genome_pos	Gene_pos	AA_pos	Phase	Ref->Alt	RefCodon->AltCodon	RefAA->AltAA	AltRatio	ErrorRatio	Depth	AverageDepth	RefRead	AltRead	AdjacentAlt	LLR	Pvalue(fisher)	RNAeditGC	GeneGC	Exon_num	Exon_start	Exon_end	Function
ath	-	nad2	81279	19	7	1	C->T	CGG->TGG	R->W	0.18914845516202	0.000623221129939135	1327	1271.3	1076	251	9	1000	1.15601898625219e-81	55.00	39.93	5	81297,80132,333105,330306,328078,	81113,79740,332945,329735,327890,	NADH dehydrogenase subunit 2
ath	-	nad2	81270	28	10	1	C->T	CCA->TCA	P->S	0.552412645590682	0.000497481358414786	1202	1048.9	538	664	31	1000	1.19926472602608e-256	55.00	39.93	5	81297,80132,333105,330306,328078,	81113,79740,332945,329735,327890,	NADH dehydrogenase subunit 2
ath	-	nad2	81239	59	20	2	C->T	TCC->TTC	S->F	0.973550356052899	0.00103016776810268	983	930.3	26	957	30	1000	0	55.00	39.93	5	81297,80132,333105,330306,328078,	81113,79740,332945,329735,327890,	NADH dehydrogenase subunit 2
ath	-	nad2	81209	89	30	2	C->T	TCC->TTC	S->F	0.930120481927711	0.000477897934247639	415	417.6	29	386	1	1000	1.33197151054271e-203	45.00	39.93	5	81297,80132,333105,330306,328078,	81113,79740,332945,329735,327890,	NADH dehydrogenase subunit 2
......

3.2.2	Table 4 - coco_cp.out (partial).
#Chr	Strand	Name	Genome_pos	Gene_pos	AA_pos	Sample_num	Phase	SRR1063404_cp	SRR1063407_cp	SRR1137438_cp	SRR1173229_cp	SRR1265939_cp	SRR1273070_cp	SRR1273180_cp	SRR606452_cp	Exon_num	Exon_start	Exon_end	Function
cp	-	rpl23	2143	71	24	2	2	C->T:TCT->TTT:S->F:0.920792079207921:0:202:197.2:16:186:18:1000:4.31666359754032e-97:25.00:37.59						C->T:TCT->TTT:S->F:0.968599033816425:0.000216450216450216:414:451.9:13:401:18:1000:4.2221509942939e-224:25.00:37.59		1	2213,	1932,	ribosomal protein L23
cp	-	trnM-CAU	2399	49	?	1	?						C->T:?->?:?->?:0.75:0:12:12:3:9:27:28.7754103989239:0.000168259523489029:40.00:45.95			1	2447,	2374,	tRNA
cp	+	ycf2	2886	371	124	1	2						G->A:AGA->AAA:R->K:0.6:0.025:20:19.9:8:12:416:11.3538827969065:2.25475753840605e-05:45.00:37.66			1	2516,	9400,	hypothetical chloroplast RF21
cp	+	ycf2	2997	482	161	2	2				C->T:CCG->CTG:P->L:0.333333333333333:0:18:15:12:6:5:20.8308785164354:0.00953079178885632:50.00:37.66		C->T:CCG->CTG:P->L:0.75:0:12:13:3:9:5:28.7754103989239:0.000168259523489029:50.00:37.66			1	2516,	9400,	hypothetical chloroplast RF21
......

3.2.3	Table 5 - coco_cp.change_matrix (number and percent).
#Prefix	A->C(percent)	A->G(percent)	A->T(percent)	C->G(percent)	C->T(percent)	C->A(percent)	G->T(percent)	G->A(percent)	G->C(percent)	T->A(percent)	T->C(percent)	T->G(percent)
SRR1063404_cp	0(0)	0(0)	0(0)	0(0)	96(100)	0(0)	0(0)	0(0)	0(0)	0(0)	0(0)	0(0)
SRR1063407_cp	0(0)	0(0)	0(0)	1(1.52)	59(89.39)	0(0)	0(0)	5(7.58)	0(0)	0(0)	1(1.52)	0(0)
SRR1137438_cp	0(0)	0(0)	0(0)	0(0)	50(98.04)	0(0)	0(0)	0(0)	0(0)	0(0)	1(1.96)	0(0)
SRR1173229_cp	2(1.55)	5(3.88)	0(0)	2(1.55)	103(79.84)	0(0)	0(0)	10(7.75)	0(0)	1(0.78)	3(2.33)	2(1.55)
SRR1265939_cp	0(0)	0(0)	0(0)	0(0)	44(84.62)	2(3.85)	0(0)	2(3.85)	0(0)	0(0)	1(1.92)	0(0)
SRR1273070_cp	5(2.79)	9(5.03)	1(0.56)	1(0.56)	85(47.49)	3(1.68)	3(1.68)	43(24.02)	3(1.68)	2(1.12)	16(8.94)	4(2.23)
SRR1273180_cp	0(0)	0(0)	0(0)	0(0)	87(97.75)	0(0)	0(0)	1(1.12)	0(0)	0(0)	0(0)	0(0)
SRR606452_cp	0(0)	0(0)	1(1.47)	0(0)	65(95.59)	0(0)	0(0)	1(1.47)	0(0)	0(0)	0(0)	1(1.47)

3.2.4	Table 6 - coco_cp.stat #(1:hydrophobic2hydrophobic,2:hydrophilic2hydrophilic,3:hydrophobic2hydrophilic,4:hydrophilic2hydrophobic).
#Prefix	Total	Other	Phase(1,2,3)	Silent	NonSilent	A->C(silent)	A->C(non_silent)	A->C(1,2,3,4)	A->G(silent)	A->G(non_silent)	A->G(1,2,3,4)	A->T(silent)	A->T(non_silent)	A->T(1,2,3,4)	C->G(silent)	C->G(non_silent)	C->G(1,2,3,4)	C->T(silent)	C->T(non_silent)	C->T(1,2,3,4)	C->A(silent)	C->A(non_silent)	C->A(1,2,3,4)	G->T(silent)	G->T(non_silent)	G->T(1,2,3,4)	G->A(silent)	G->A(non_silent)	G->A(1,2,3,4)	G->C(silent)	G->C(non_silent)	G->C(1,2,3,4)	T->A(silent)	T->A(non_silent)	T->A(1,2,3,4)	T->C(silent)	T->C(non_silent)	T->C(1,2,3,4)	T->G(silent)	T->G(non_silent)	T->G(1,2,3,4)
SRR1063404_cp	96	0	15,73,9	9	87	0	0	0,0,0,0	0	0	0,0,0,0	0	0	0,0,0,0	0	0	0,0,0,0	9	87	15,10,1,61	0	0	0,0,0,0	0	0	0,0,0,0	0	0	0,0,0,0	0	0	0,0,0,0	0	0	0,0,0,0	0	0	0,0,0,0	0	0	0,0,0,0
SRR1063407_cp	66	0	14,47,7	7	59	0	0	0,0,0,0	0	0	0,0,0,0	0	0	0,0,0,0	0	1	1,0,0,0	6	53	9,7,1,36	0	0	0,0,0,0	0	0	0,0,0,0	0	5	3,1,1,0	0	0	0,0,0,0	0	0	0,0,0,0	1	0	0,0,0,0	0	0	0,0,0,0
SRR1137438_cp	51	0	10,37,5	4	47	0	0	0,0,0,0	0	0	0,0,0,0	0	0	0,0,0,0	0	0	0,0,0,0	3	47	9,6,1,31	0	0	0,0,0,0	0	0	0,0,0,0	0	0	0,0,0,0	0	0	0,0,0,0	0	0	0,0,0,0	1	0	0,0,0,0	0	0	0,0,0,0
SRR1173229_cp	129	1	36,82,32	20	108	1	1	1,0,0,0	2	3	2,0,0,1	0	0	0,0,0,0	1	1	1,0,0,0	14	89	18,11,3,57	0	0	0,0,0,0	0	0	0,0,0,0	0	10	4,1,5,0	0	0	0,0,0,0	0	1	1,0,0,0	2	1	1,0,0,0	0	2	0,0,2,0
SRR1265939_cp	52	3	12,37,8	4	45	0	0	0,0,0,0	0	0	0,0,0,0	0	0	0,0,0,0	0	0	0,0,0,0	3	41	9,6,0,26	0	2	0,1,1,0	0	0	0,0,0,0	0	2	1,0,1,0	0	0	0,0,0,0	0	0	0,0,0,0	1	0	0,0,0,0	0	0	0,0,0,0
SRR1273070_cp	179	4	64,87,72	49	126	3	2	0,1,0,1	2	7	1,2,0,4	1	0	0,0,0,0	0	1	1,0,0,0	23	62	15,13,4,30	0	3	0,2,1,0	2	1	1,0,0,0	10	33	9,14,10,0	0	3	2,0,0,1	1	1	1,0,0,0	7	9	2,0,6,1	0	4	2,0,2,0
SRR1273180_cp	89	1	15,69,4	5	83	0	0	0,0,0,0	0	0	0,0,0,0	0	0	0,0,0,0	0	0	0,0,0,0	5	82	14,10,1,57	0	0	0,0,0,0	0	0	0,0,0,0	0	1	0,0,1,0	0	0	0,0,0,0	0	0	0,0,0,0	0	0	0,0,0,0	0	0	0,0,0,0
SRR606452_cp	68	0	15,48,8	8	60	0	0	0,0,0,0	0	0	0,0,0,0	1	0	0,0,0,0	0	0	0,0,0,0	6	59	10,6,2,41	0	0	0,0,0,0	0	0	0,0,0,0	0	1	1,0,0,0	0	0	0,0,0,0	0	0	0,0,0,0	0	0	0,0,0,0	1	0	0,0,0,0

3.2.5	Table 7 - coco_cp_DPR.out (partial)
Chr	Position	Proportion1	Alt1	Ref1	Proportion2	Alt2	Ref2	Fold(log2)	Pvalue
#SRR1063404_cp vs SRR1063407_cp
cp	136348	0.5	6	6	0	0	0	NA	1
cp	143536	0.936986301369863	342	23	1	24	0	0.09	0.382540126150565
cp	76241	0.571428571428571	40	30	1	7	0	0.81	0.0386037932652878
cp	32743	0.466531440162272	230	263	0.8125	13	3	0.80	0.00902812345271057
......

3.2.6	Table 8 - The RNA editing sites identification between REDO and REDItools in three analyzed RNA-seq datasets (without applying complicated filter model).
Samples 	Software	Type	True	Total	True+	True-	False+	False-	Precision	Sensitivity	Specificity
Coconut cp1
(SRR1063404)	REDO	CTa	75	83	67	33349	16	8	80.72%	89.33%	99.95%
	REDItools			92	69	33340	23	6	75.00%	92.00%	99.93%
	Common			80	65	33352	15	10	81.25%	86.67%	99.96%
	REDO	Allb		84	67	33348	17	8	79.76%	89.33%	99.95%
	REDItools			98	69	33334	29	6	70.41%	92.00%	99.91%
	Common			84	67	33347	17	8	79.76%	89.33%	99.95%
Coconut cp2
(SRR1273070)	REDO	CT		54	23	33378	31	52	42.59%	30.67%	99.91%
	REDItools			57	22	33375	35	53	38.60%	29.33%	99.90%
	Common			49	20	33383	29	55	40.82%	26.67%	99.90%
	REDO	All		112	23	33320	89	52	20.54%	30.67%	99.73%
	REDItools			121	22	33311	99	53	18.18%	29.33%	99.70%
	Common			98	20	33334	78	55	20.41%	26.67%	99.74%
Arabidopsis mt
(SRR2079784)	REDO	CT	428c	413	358	29041	55	70	86.68%	83.64%	99.81%
	REDItools			430	362	29024	68	66	84.19%	84.58%	99.77%
	Common			409	356	29045	53	72	87.04%	83.18%	99.81%
	REDO	All		430	358	29024	72	70	83.26%	83.64%	99.75%
	REDItools			469	362	28985	86	66	80.80%	84.58%	99.70%
	Common			426	356	29028	70	72	83.57%	83.18%	99.76%

Note: REDItools was run with “-c 4” option (total reads coverage>=4) and the result was filtered by p-value<=0.01. RED was run using “-a 0, -s 0” option (minimum alt distance==0 and minimum junction anchor size==0) for mt data and default for cp data without applying the complicated filter model. a: CT means C->T type. b: All means all possible RNA editing types. c: There are 455 C->T RNA editing sites in Arabidopsis by cDNA clones sequencing while 428 C->T RNA editing sites are located in CDS region. We only used these 428 C->T RNA editing sites for mt RNA editing sites test.

3.2.7	Table 9 - The RNA editing sites identification between REDO and REDItools in three analyzed RNA-seq datasets (with applying complicated filter model).
Samples 	Software	Type	True	Total	True+	True-	False+	False-	Precision	Sensitivity	Specificity
Coconut cp1
(SRR1063404)	REDO	CTa	75	83	67	16	8	80.72%	89.33%	33349	99.95%
	REDItools			92	69	23	6	75.00%	92.00%	33340	99.93%
	Common			80	65	15	10	81.25%	86.67%	33352	99.96%
	REDO	Allb		83	67	16	8	80.72%	89.33%	33349	99.95%
	REDItools			98	69	29	6	70.41%	92.00%	33334	99.91%
	Common			80	65	15	10	81.25%	86.67%	33352	99.96%
Coconut cp2
(SRR1273070)	REDO	CT		44	22	22	53	50.00%	29.33%	33388	99.93%
	REDItools			57	22	35	53	38.60%	29.33%	33375	99.90%
	Common			39	19	20	56	48.72%	25.33%	29415	99.93%
	REDO	All		71	22	49	53	30.99%	29.33%	29383	99.83%
	REDItools			121	22	99	53	18.18%	29.33%	29333	99.66%
	Common			65	19	46	56	29.23%	25.33%	29389	99.84%
Arabidopsis mt
(SRR2079784)	REDO	CT	428c	406	356	50	72	87.68%	83.18%	29048	99.83%
	REDItools			430	362	68	66	84.19%	84.58%	29024	99.77%
	Common			402	354	48	74	88.06%	82.71%	29052	99.84%
	REDO	All		419	356	63	72	84.96%	83.18%	29035	99.78%
	REDItools			469	362	86	66	80.80%	84.58%	28985	99.70%
	Common			415	354	61	74	85.30%	82.71%	29039	99.79%

Note: REDItools was run with “-c 4” option (total reads coverage>=4) and the result was filtered by p-value<=0.01. RED was run using “-a 0, -s 0” option (minimum alt distance==0 and minimum junction anchor size==0) for mt data and default for cp data with applying the complicated filter model. a: CT means C->T type. b: All means all possible RNA editing types. c: There are 455 C->T RNA editing sites in Arabidopsis by cDNA clones sequencing while 428 C->T RNA editing sites are located in CDS region. We only used these 428 C->T RNA editing sites for mt RNA editing sites test.

SUBPROGRAM
==========

We also provide a subprogram named fish.pl for processing the REDO result tables easily and quickly. It can be used to find common and unique sites in two samples or remove genome variant in REDO result file by comparing with genome variant file. 

It can get unique bait/baits (one or multiple items in column) in bait file at first (such as chromosome, position, strand, and others), then these bait/baits are used to search the corresponding column in fish file in two measures, existing in bait file (de-fault) or not exiting in bait file ("-contrary" option).

CONTACT
=======

If you have any question or suggestion, please contact us.

Shuangyang Wu &Wanfei Liu & Qiang Lin
Email: ,<wushy@big.ac.cn> & <liuwf@big.ac.cn> & <linq@big.ac.cn>

If you find REDO is useful in your research. Please cite 
REDO: RNA Editing Detection in Plant Organelles Based on Variant Calling Results.J Comput Biol. 2018 Apr 11. doi: 10.1089/cmb.2017.0214
PMID: 29641228 

