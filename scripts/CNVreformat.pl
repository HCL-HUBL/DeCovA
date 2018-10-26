#!/usr/bin/perl

=pod

CNV_output_tuning:

	1): from : smpl.summary.txt:

	ex:
	#Region	CNV-type	N_clean_intervals	clean_ratio_to_avg	N_dirty_intervals	dirty_ratio_to_avg	overlapping_Samples	qual
	chr1:1361474-1361790	DUP	1	1.448	.	.	1	low

	-> #Region	Gene(s)	gene-related-infos	CNV-type	CNV-qual	smplX_N_intervals	smplX_ratio_to_avg	overlapping_Samples


	2): from CNV_fam files:

	ex:
	#Region	CNV-type	boundaries_smpl1(A)	boundaries_smpl2(H)	boundaries_smpl3(H)	in_other_groups(tot=10)
	chr1:1361474-1361790	dup	1361474-1361790	.	.	0

	-> #Region	Gene(s)	gene-related-infos	type	qual	smplX_N_intervals	smplX_ratio_to_avg	overlapping_Samples	in_other_groups

	for N_intervals : sum of intervals
	for ratio: mean


	then, from : smpl.summary.txt:
	ex:
	#Region	CNV	N_clean_intervals	clean_ratio_to_avg	N_dirty_intervals	dirty_ratio_to_avg	overlapping_Samples	qual
	chr1:1361474-1361790	DUP	1	1.448	.	.	1	low

	-> get clean-val, and dirty-val if exists
	-> overlapping_Samples : take max val


	from smpl.highQual.txt:

	#Chrom	Start	End	Length	Info	interval_Order	CNV_type	ratio_to_avg	ratio_to_std	occurences	min	max	med	avg	std
	chr1	16370965	16371114	150 bp	CLCNKB	2359	DUP	1.595	8.903	4	206.0	379.6	241.1	238.1	15.9

	-> Infos?
	or annotation with other tool: with snpEff?, with DeCovA?: requires a bed format

=cut

use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use File::Find;



##get parameters
Getopt::Long::Configure ("bundling");
my %opt = ();
GetOptions (\%opt,
	"file|f=s@",
	"fileLs|F=s",
	"typeC|C=i",
	"posC|P=i",
	"dir|d=s@",
	"suffix|s=s@",
	"int|i=s",
	"help|h")
or die "!!! err in command line arguments\n$!\n";

unless(%opt) { usage(); }
if (exists $opt{help}) { usage(); }
sub usage {
die "usage: perl script.pl [options]
or: perl script.pl [options] CNVfile1 CNVfile2 ...
    options:
    -f / --file [file]: input CNV files (comma separated, or set several times),
                        either in bed format (chr<tab>start(0-based)<tab>end) or region format (chr:start(1-based)-end)
                        if region format, opt --posC is mandatory
    -F / --fileLs [file]: file with list of input CNV files, one per line
    -C / --typeC [int]: index(1-based) of CNV type field (def: 2nd col)
    -P / --posC [int]: index(1-based) of position field, if CNV file in region format
    -d / --dir [dir]: directory(ies) where to look for input files (comma separated, or set several times)
    -s / --suffix [str]: used to select input files in dir(s) (comma separated, or set several times)(def: .summary.txt)
    -i / --int [file]: file with infos and raw fold-changes per intervals
    -h / --help: this help text\n";
}

##input files
my @dir = ();
if (exists $opt{dir}) {
	foreach (@{ $opt{dir} }) {
		push(@dir, split(/,/, $_));
	}
	foreach (@dir) {
		unless (-d $_) { die "!!! err: $_ directory not found\n"; }
		unless ($_ =~ /\/$/) { $_ .= "/"; } 
	}
}

my @suff = ();
if (exists $opt{suffix}) {
	foreach (@{ $opt{suffix} }) {
		push(@suff, split(/,/, $_));
	}
}
unless (@suff) { push(@suff , ".summary.txt"); }
if (@dir) {
	print STDERR "## looking for files ending with:\n";
	foreach (@suff) { print STDERR "\t$_\n"; }
}

my $fh;
my (%File,@Files);
foreach (@ARGV) {
	push(@Files, $_);
}
if (exists $opt{fileLs}) {
	open($fh, "<", $opt{fileLs}) or die $!;
	while (my $line = <$fh>) {
		unless ($line =~ /^\s*$/) {
			$line =~ s/^\s+//; $line =~ s/\s+$//;
			push (@Files, $line);
		}
	}
	close $fh;
}
if (exists $opt{file}) {
	foreach (@{ $opt{file} }) {
		push(@Files, split(/,/, $_));
		}
	}
foreach (@dir) {
	find(\&subDir, $_);
}
sub subDir {
	foreach my $s (@suff) {
		if (-f $_ && $_=~/$s$/) { 
			push(@Files, $File::Find::name);
			$File{$File::Find::name}{"suff"} = $s;
		}
	}
}

if (@Files) {
	print STDERR "## file(s) to be analysed:\n";
	foreach (@Files) {
		if (exists $File{$_}{"suff"}) {
			($File{$_}{"name"},$File{$_}{"dir"}) = fileparse($_, $File{$_}{"suff"});
			print STDERR "\t".$File{$_}{"name"}.$File{$_}{"suff"}."\n";
		} else {
			($File{$_}{"name"},$File{$_}{"dir"}) = fileparse($_, $suff[0]);
			print STDERR "\t".$File{$_}{"name"}.$suff[0]."\n";
		}
		unless (-f $_) { die "!!! err: $_ not found\n"; }
	}
} else {
	die "!!! err: no input file found/provided\n";
}

my $CNVfmt = "bed";
my (%idx);
$idx{"pos"} = 0;
if ($opt{posC}) {
	$CNVfmt = "reg";
	$idx{"pos"} = $opt{posC} - 1;
}
$idx{"type"} = 1;
if ($opt{typeC}) {
	$idx{"type"}  = $opt{typeC} - 1;
}

my $rawIntervalFile;
if (exists $opt{"int"}) {
	if (-f $opt{"int"}) { $rawIntervalFile = $opt{"int"};
	} else { die "!!! err: ".$opt{"int"}." not found\n"; }
}







my (%rawIntervals,%sortedRawIntervals);
if ($rawIntervalFile) {
	my (@chrOrder, @chrHash,$prevChr,$prevPos);
	print STDERR "## reading $rawIntervalFile\n";
	open($fh, "<", $rawIntervalFile) or die $!;
	my $line = <$fh>;
	#Chrom   Start   End     Infos   Numero_Region   15A2728 15A5068 17A2472 17A2473 17A2474 17A3622 17A5566 17A5854 17A5859 17A5860 17A5941 17A5942 17A6155 17A6209 17A6485 17A6486 17A6487 17A6512 17A6513 18A684  18A685  18A724  18A725  18A726  Statut_Region
	#chr1    17345   17452   interG:OR4F5:NM_001005484--51638        0       0.834(-0.942)   1.302(1.714)    0.920(-0.457)   0.953(-0.267)   0.988(-0.067)   0.979(-0.118)   1.009(0.054)    1.179(1.014)    1.655(3.722)    0.895(-0.594)   0.699(-1.710)   1.096(0.545)    1.067(0.380)    1.135(0.766)    1.070(0.396)    1.207(1.174)    1.007(0.039)    0.853(-0.834)   0.592(-2.319)   1.297(1.685)    0.941(-0.335)   0.866(-0.760)   0.508(-2.791)   1.112(0.636)    OK
	my %rawIntervalIdx;
	my $i = 0;
	foreach (split(/\t/,$line)) {
		$rawIntervalIdx{$_} = $i; #print STDERR "$_ ; $i\n";
		$i++;
	}
	my %foundSmpl;
	foreach my $smpl (@Samples) {
		if (exists $rawIntervalIdx{$smpl}) {
		$foundSmpl{$smpl} = 1
		}
	}
	my ($chr,$start,$end,@tab);
	while ($line = <$fh>) {
		if ($line !~ /^#/ && $line =~ /^(?:[Cc][Hh][Rr])?(.+?)\t(\d+)\t(\d+)\t/) {
			$chr = $1;
			$start = $2+1;	# Bed -> 1-based
			$end = $3;
			if ($prevChr) {
				if ($chr ne $prevChr) {
					if (exists $chrHash{$chr}) {
						print STDERR "$rawIntervalFile not sorted\n";
					} else {
						push(@chrOrder, $chr);
						$chrHash{$chr} = 1;
					}
				}
			} else { $chr = $prevChr; }
			$rawIntervals{$chr}{$start}{"end"} = $end;
			chomp $line;
			@tab = split(/\t/,$line);
			if ( $tab[$rawIntervalIdx{"Infos"}] =~ /^(.+?):/) {
				$rawIntervals{$chr}{$start}{"gene"} = $1;
			}
			foreach my $smpl (keys%foundSmpl) {	 #print STDERR "$smpl\n"; unless (exists $rawIntervalIdx{$smpl}) {print STDERR "not found\n";exit;}
				if ($tab[$rawIntervalIdx{$smpl}] =~ /^(\d+\.?\d*)/) {
					$rawIntervals{$chr}{$start}{"diff"}{$smpl} = $1;
				}
			}
		}
	}
	close $fh;
	foreach my $chr (keys%rawIntervals) {		#print STDERR "raw : $chr\n";
		foreach my $start (sort{$a<=>$b}keys%{ $rawIntervals{$chr} }) {
			push(@{ $sortedRawIntervals{$chr} }, $start);
		}
	}
}
