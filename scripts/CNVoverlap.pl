#!/usr/bin/perl

## print CNVs from 1 file if overlapped by CNVs from other files

my $VERSION = "1.4";

use warnings;
use strict;
use Getopt::Long;


##get parameters
Getopt::Long::Configure ("bundling");
my%opt = ();
GetOptions (\%opt,
	"ref|r=s",
	"test|t=s@",
	"minRef|m=f",
	"minTest|M=f",
	"choice|c=s",
	"print|p=s",
	"help|h")
or die "!!! err in command line arguments\n$!\n";

unless(%opt) { usage(); }
if (exists $opt{help}) { usage(); }
sub usage {
die "usage: perl script.pl [options]
	options:
	-r / --ref [file]: input CNV file used to report matched intervals; format sample<tab>CNV_type<tab>pos(chr:start-end)
	-t / --test [file]: input CNV file(s) used to look for matching intervals (comma separated, or set several times)(same format as above)
	-m / --minRef [0-1]: CNVs from ref file shared if overlapped by x fraction of their length (def: shared if overlapped by 1 base)
	-M / --minTest [0-1]: CNVs from test file(s) shared if overlapped by x fraction of their length (def: shared if overlapped by 1 base)
	-c / --choice [A/O]: if both minRef and minTest selected, CNVs shared if both are overlapped enough (And), or if only one is overlapped (Or)
	-p / --print [A/O]: print All lines from ref file (def), or only those Overlapped by tested ones
	-h: help\n";
}

##arg
my($ref);
if (exists $opt{"ref"}) {
	if (-f $opt{"ref"}) { $ref = $opt{"ref"}; }
	else { die "ref file $opt{ref} not found\n"; }
	}
else { die "ref file required\n"; }

my(%File,@Files);
if (exists $opt{"test"}) {
	foreach (@{ $opt{"test"} }) {
		if (-f $_) { push(@Files, split(/,/, $_)); }
		else { die "file $_ not found\n"; }
		}
	}
unless (@Files) { die "no tested CNV bed file provided\n"; }

my$minRatioRef = 0;
if (exists $opt{minRef}) { $minRatioRef = $opt{minRef}; }

my$minRatioTest = 0;
if (exists $opt{minTest}) { $minRatioTest = $opt{minTest}; }

my$choice = "and";
if (exists $opt{choice}) {
	if ($opt{choice} =~ /^a/i) { $choice = "and"; }
	elsif ($opt{choice} =~ /^o/i) { $choice = "or"; }
	else { die "opt \"--choice\" not recognized (\"A\" or \"O\")\n"; }
	}

my$printAll = 1;
if ($opt{"print"}) {
	if ($opt{"print"} =~ /^a/i) { $printAll = 1; }
	elsif ($opt{"print"} =~ /^o/i) { $printAll = ""; }
	else { die "error with opt -p (enter \"A\" or \"O\"\n"; }
	}


my(%idx);
$idx{"smpl"} = 0;
$idx{"type"} = 1;
$idx{"pos"} = 2;


####
my(%RefCNV);
print STDERR"reading ref file: $ref\n";
open(my$fh, "<", "$ref") or die "cannot open $ref";
while (my$line = <$fh>) {
	if ($line !~ /^#/ && $line !~ /^\s*$/) {
		chomp $line;
		my@tab = split(/\t/,$line);
		my@pos = $tab[$idx{"pos"}] =~ m/^(.+):(\d+)-(\d+)$/;	#print STDERR "$tab[2] : $pos[0] , $pos[1] , $pos[2]\n";
		##		$RefCNV{sample}{chrom}{start} 
		$RefCNV{$tab[$idx{"smpl"}]}{$pos[0]}{$pos[1]}{"end"} = $pos[2];
		$RefCNV{$tab[$idx{"smpl"}]}{$pos[0]}{$pos[1]}{"type"} = $tab[$idx{"type"}];
		}
	}
close $fh;

####
my%trueRefOverL;			#push(@{ ${$trueRefOverL}{"$chr:$refStart-$refEnd"} }, "$testFile($testStart-$testEnd)")
foreach my$file (@Files) {
	my$fileName = $file;
	$fileName =~ s/^.+\///;
	$fileName =~ s/\.(txt|bed)$//;
	my%testCNV;
	print STDERR "reading tested file: $file\n";
	open($fh, "<", "$file") or die "cannot open $file\n";
	while (my$line = <$fh>) {
		if ($line !~ /^#/ && $line !~ /^\s*$/) {
			chomp $line;
			my@tab = split(/\t/,$line);
			my$smpl = $tab[$idx{"smpl"}];
			$smpl =~ s/-(.+?)$//;	#print STDERR "$tab[0] : $1\n";
			my@pos = $tab[$idx{"pos"}] =~ m/^(.+):(\d+)-(\d+)$/;
			$testCNV{$smpl}{$pos[0]}{$pos[1]}{"end"} = $pos[2];
			$testCNV{$smpl}{$pos[0]}{$pos[1]}{"type"} = $tab[$idx{"type"}];
			}
		}
	close $fh;
	foreach my$sample (keys%RefCNV) {
		foreach my$chr (keys%{ $RefCNV{$sample} }) {
			if (exists $RefCNV{$sample}{$chr} && exists $testCNV{$sample}{$chr}) {
				my@refStarts = sort{$a<=>$b}keys%{ $RefCNV{$sample}{$chr} };
				my$c = 0;	#idx of @refStarts
				foreach my$testStart (sort{$a<=>$b}keys%{ $testCNV{$sample}{$chr} }) {
					my$testEnd = $testCNV{$sample}{$chr}{$testStart}{"end"};
					if ($testStart > $RefCNV{$sample}{$chr}{$refStarts[-1]}{"end"}) { last; }
					while ( ($c < $#refStarts) && ($testStart > $RefCNV{$sample}{$chr}{$refStarts[$c]}{"end"}) ) {
						$c++; 
						}
					my$c2 = $c;
					while ( ($c2 < scalar@refStarts) && ($testEnd >= $refStarts[$c2]) ) {
						if ($minRatioRef || $minRatioTest) {
							testMinR($minRatioRef,$minRatioTest,$choice,$chr,$refStarts[$c2],$RefCNV{$sample}{$chr}{$refStarts[$c2]}{"end"},$RefCNV{$sample}{$chr}{$refStarts[$c2]}{"type"},$fileName,$testStart,$testEnd,$testCNV{$sample}{$chr}{$testStart}{"type"},\%{ $trueRefOverL{$sample} }); 
							}
						else {
							push(@{ $trueRefOverL{$sample}{"$chr:".$refStarts[$c2]."-".$RefCNV{$sample}{$chr}{$refStarts[$c2]}{"end"}}{$fileName} }, "$testStart-$testEnd");
							}
						$c2++;
						}
=pod
					if ($testEnd >= $refStarts[$c2]) {
						if ($minRatioRef || $minRatioTest) {
							testMinR($minRatioRef,$minRatioTest,$choice,$chr,$refStarts[$c2],$RefCNV{$sample}{$chr}{$refStarts[$c2]}{"end"},$RefCNV{$sample}{$chr}{$refStarts[$c2]}{"type"},$fileName,$testStart,$testEnd,$testCNV{$sample}{$chr}{$testStart}{"type"},\%{ $trueRefOverL{$sample} });
							}
						else {
							push(@{ $trueRefOverL{$sample}{"$chr:".$refStarts[$c2]."-".$RefCNV{$sample}{$chr}{$refStarts[$c2]}{"end"}}{$fileName} }, "$testStart-$testEnd");
							}
						}
=cut
					}
				}
			}
		}
	}

####
open($fh, "<", "$ref") or die "cannot open $ref\n";
while (my$line = <$fh>) {
	chomp $line;
	if ($line =~ /^#/) {
		print "$line\tOccurences\n";
		}
	elsif($line =~ /^\s*$/) { print "$line\n"; }
	else {
		my@tab = split(/\t/,$line);
		my$smpl = $tab[$idx{"smpl"}];
		my$pos = $tab[$idx{"pos"}];
		if ($printAll) {
			print "$line";
			if (exists $trueRefOverL{$smpl}{$pos}) {
				print "\t".scalar(keys%{ $trueRefOverL{$smpl}{$pos} });
				foreach my$testFile (sort(keys%{ $trueRefOverL{$smpl}{$pos} })) {
					print "\t$testFile";
					foreach (@{ $trueRefOverL{$smpl}{$pos}{$testFile} }) { print " ($_)"; }
					}
				}
			else { print "\t0"; }
			print "\n";
			}
		else {
			if (exists $trueRefOverL{$smpl}{$pos}) {
				print "$line\t".scalar(keys%{ $trueRefOverL{$smpl}{$pos} });
				#foreach (@{ $trueRefOverL{$smpl}{$pos} }) { print "\t$_"; }
				foreach my$testFile (sort(keys%{ $trueRefOverL{$smpl}{$pos} })) {
					print "\t$testFile";
					foreach (@{ $trueRefOverL{$smpl}{$pos}{$testFile} }) { print " ($_)"; }
					}
				print "\n";
				}
			}
		}
	}


exit;


####
sub testMinR {
my($minRatioRef,$minRatioTest,$choice,$chr,$refStart,$refEnd,$refCNV,$fileName,$testStart,$testEnd,$testCNV,$trueRefOverL) = @_;
my$RefOverL = 1;
if ($minRatioRef && $minRatioTest) { $RefOverL = 2; }
my($minStart, $minEnd);
if ($testStart > $refStart) { $minStart = $testStart; }
else { $minStart = $refStart; }
if ($testEnd < $refEnd) { $minEnd = $testEnd; }
else { $minEnd = $refEnd; }
if ($minRatioRef && (($minEnd-$minStart)/($refEnd-$refStart)) < $minRatioRef) { $RefOverL--; }
if ($minRatioTest && (($minEnd-$minStart)/($testEnd-$testStart)) < $minRatioTest) { $RefOverL--; }
if ($minRatioRef && $minRatioTest) {
	if ($choice eq "and") {
		if ($RefOverL < 2) { $RefOverL = 0; }
		}
	}
if ($RefOverL && $refCNV =~ /^$testCNV/i) {
	push(@{ ${$trueRefOverL}{"$chr:$refStart-$refEnd"}{$fileName} }, "$testStart-$testEnd");
	}
}



