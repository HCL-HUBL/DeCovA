#!/usr/bin/perl


## print intersections of CNVs from all files provided or found


my $VERSION = "1.5e";

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
	"bed|b=s",
	"outName|o=s",
	"name|n=s@",
	"nameLs|N=s",
	"fam=s",
	"minRatio|m=f",
	"choice|c=s",
	"union|u",
	"transitivity|t",
	"help|h")
or die "!!! err in command line arguments\n$!\n";

unless(%opt) { usage(); }
if (exists $opt{help}) { usage(); }
sub usage {
die "usage: perl script.pl [options]
or: perl script.pl [options] CNVfile1 CNVfile2 ...
	options:
	-f / --file [file]: input CNV files (comma separated, or set several times)
						in bed format (chr<tab>start(0-based)<tab>end) or region format (chr:start(1-based)-end)
	-F / --fileLs [file]: file with list of input CNV files, one per line
	-C / --typeC [int]: index(1-based) of CNV type field (def: 2nd col)
	-P / --posC [int]: index(1-based) of position field, if not a bed file
	-d / --dir [dir]: directory(ies) where to look for input files (comma separated, or set several times)
	-s / --suffix [str]: used to select input files in dir(s) (comma separated, or set several times)(def: .summary.txt)
	-b / --bed [str]: intervals file to restrict analysis
	-o / --outName [str]: prefix for output files (def: none in current working dir)
	-n / --name : [str]: list of sample names to select (comma separated or set several times) (def: all found)
				names can be followed by a status: name:A/H for Affected/Healthy
	-N / --nameLs : [str]: file with list of sample names to select, one per line (def: all found)
				names can be followed by a status: name:A/H for Affected/Healthy
	--fam : [str]: file with list of families, one per line, then sample names and status
				famName<tab>sample1:A/H<tab>sample2:A/H...
	-m / --minRatio : [0-1]: CNVs from 2 samples shared if overlapped by x fraction of their length (def: shared if overlapped by 1 base)
	-c / --choice : [A/O]: CNVs shared if both samples are overlapped enough (And), or if only one is significantly overlapped (Or)
	-u / --union : overlap == union of intervals (def: intersection of intervals)
	-t / --transitivity : if A overlaps B and B overlaps C, then A overlaps C (def: no)
	-h: help\n";
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
	push(@Files, "$_");
	}
if (exists $opt{fileLs}) {
	open($fh, "<", $opt{fileLs}) or die $!;
	while (my $line = <$fh>) {
		unless ($line =~ /^\s*$/) {
			$line =~ s/^\s+//; $line =~ s/\s+$//;
			push (@Files, "$line");
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
	print STDERR "file(s) to analysed:\n";
	foreach (@Files) {
		if (exists $File{$_}{"suff"}) {
			($File{$_}{"name"},$File{$_}{"dir"}) = fileparse($_, $File{$_}{"suff"});
			print STDERR "\t".$File{$_}{"name"}.$File{$_}{"suff"}."\n";
			}
		else {
			($File{$_}{"name"},$File{$_}{"dir"}) = fileparse($_, $suff[0]);
			print STDERR "\t".$File{$_}{"name"}.$suff[0]."\n";
			}
		unless (-f $_) { die "!!! err: $_ not found\n"; }
		}
	}
else { die "!!! err: no input file found/provided\n"; }

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

my ($bedFile,$Bed1based);
if (exists $opt{bed}) {
	if (-f $opt{bed}) { $bedFile = $opt{bed}; }
	else { die "!!! err: ".$opt{bed}." not found\n"; }
	}

##outName
my $outName = "CNVoverlap";
if (exists $opt{outName}) {
	if (-d $opt{outName}) {	#existing dir, no file name
		$outName =~ s/\/$//;
		$outName = $opt{outName}."/".$outName;
		}
	else {
		$outName = $opt{outName};
		my ($outFile,$outPath) = fileparse($outName);
		unless (-d $outPath) {
			die "!!! err: $outPath dir not found\n";
			}
		$outName =~ s/_$//;
		}
	}

my (@selectedSmpl,%selectedSmpl,%Fams);
if (exists $opt{nameLs}) {
	open($fh, "<", $opt{nameLs}) or die $!;
	while (my $line = <$fh>) {
		unless ($line =~ /^\s*$/ || $line =~ /^#/) {
			$line =~ s/^\s+//; $line =~ s/\s+$//;
			push (@selectedSmpl, "$line");
			}
		}
	close $fh;
	}
if (exists $opt{name}) {
	foreach (@{ $opt{name} }) {
		push(@selectedSmpl, split(/,/, $_));
		}
	}
foreach (@selectedSmpl) {
	if ($_ =~ /^(.+)\s*:\s*([ah])$/i) {
		my $smpl = $1; my $status = $2;
		$Fams{"fam1"}{"H_all"}{$smpl} = 1;
		push(@{ $Fams{"fam1"}{"A_all"} }, $smpl);
		if ($status =~ /a/i) { $Fams{"fam1"}{"A"}{$smpl} = 1; }
		else { $Fams{"fam1"}{"H"}{$smpl} = 1; }
		}
	else { $selectedSmpl{$_} = 1; }
	}
if (%Fams && !%selectedSmpl) { @selectedSmpl = (); }

if (exists $opt{fam}) {
	open($fh, "<", $opt{fam}) or die $!;
	while (my $line = <$fh>) {
		unless ($line =~ /^\s*$/ || $line =~ /^#/) {
			my @Fam = split(/\t/, $line);
			foreach (@Fam) { $_ =~ s/^\s+//; $_ =~ s/\s+$//; }
			my $fam = shift@Fam;
			foreach (@Fam) {
				if ($_ =~ /^(.+)\s*:\s*([ah])$/i) {
					my $smpl = $1; my $status = $2;
					$Fams{$fam}{"H_all"}{$smpl} = 1;
					push(@{ $Fams{$fam}{"A_all"} }, $smpl);
					if ($status =~ /a/i) { $Fams{$fam}{"A"}{$smpl} = 1; }
					else { $Fams{$fam}{"H"}{$smpl} = 1; }
					}
				else { die "!!! err: smpl format not recognized within ".$opt{fam}." file\n"; }
				}
			}
		}
	close $fh;
	}

my $minRatio = 0;
if (exists $opt{minRatio}) { $minRatio = $opt{minRatio}; }

my $choice = "and";
if (exists $opt{choice}) {
	if ($opt{choice} =~ /^a/i) { $choice = "and"; }
	elsif ($opt{choice} =~ /^o/i) { $choice = "or"; }
	else { die "!!! opt \"--choice\" not recognized (\"A\" or \"O\")\n"; }
	}


my $longest = 0;
if (exists $opt{union}) {
	$longest = 1;
	}
my $transitivity = 0;	#Transitive relation, a binary relation in which if A is related to B and B is related to C, then A is related to C
if (exists $opt{transitivity}) {
	$transitivity = 1;
	$longest = 1;
	}






## reads all CNV files, filling %allCNV{sample}{CNVtype}{chr}{start} = end

my (@Samples,%allCNV,%chromName);
foreach my $file (@Files) {
	my $smplName = $File{$file}{"name"};
	$smplName =~ s/^CNV_//;
	if (!@selectedSmpl || exists $selectedSmpl{$smplName}) {
		push(@Samples,$smplName);
		if (exists $allCNV{$smplName}) { die "!!! err: sample $smplName found several times\n"; }
		print STDERR "reading CNV file: $file\n";
		open($fh, "<", "$file") or die $!;
		my$someLinesOK = 0;
		while (my $line = <$fh>) {
			if ($line !~ /^#/ && $line !~ /^\s*$/) {
				chomp $line;
				my @tab = split(/\t/,$line);
				if ($CNVfmt eq "bed" && $line =~ m/^(\w+)\t(\d+)\t(\d+)/) {
					$someLinesOK++;
					my $chrName = $1;
					my $start = $2 + 1;
					my $end = $3;
					my $chr = $chrName;
					$chr =~ s/^chr//i;
					$chromName{$chr} = $chrName;
					$allCNV{$smplName}{lc($tab[$idx{"type"}])}{$chr}{$start}{"end"} = $end;
					}
				elsif ($CNVfmt eq "reg" && $tab[$idx{"pos"}] =~ m/^(\w+) *: *(\d+) *- *(\d+)/) {
					$someLinesOK++;
					my $chrName = $1;
					my $start = $2;
					my $end = $3;
					my $chr = $chrName;
					$chr =~ s/^chr//i;
					$chromName{$chr} = $chrName;
					$allCNV{$smplName}{lc($tab[$idx{"type"}])}{$chr}{$start}{"end"} = $end;
					}
				}
			}
		close $fh;
		unless ($someLinesOK) { print STDERR "no correctly formated line found in $file (format = $CNVfmt)\n"; }
		}
	}
if (@selectedSmpl) {
	foreach (@selectedSmpl) {
		unless (exists $allCNV{$_}) { die "!!! err: sample $_ not found\n"; }
		}
	}
if (%Fams) {
	foreach my $fam (keys%Fams) {
		foreach (@{ $Fams{$fam}{"A_all"} }) {
			unless (exists $allCNV{$_}) { die "!!! err: sample $_ not found\n"; }
			}
		}
	}

## order starts
my (%allStarts);
foreach my $smpl (keys%allCNV) {
	foreach my $type (keys%{ $allCNV{$smpl} }) {
		foreach my $chr (keys%{ $allCNV{$smpl}{$type} }) {
			@{ $allStarts{$smpl}{$type}{$chr} } = sort{$a<=>$b}keys%{ $allCNV{$smpl}{$type}{$chr} };
			}
		}
	}

##fill Overlap1 : exactly overlapping intervals
my %donePairs = ();
my $overlap1 = {};	# %{ ${$overlap1}{$type}{$chr}{$overStart}{$overEnd}{$smpl1}}{$S1start."-".$S1end} = 1
foreach my $smpl1 (keys%allCNV) {
	foreach my $smpl2 (keys%allCNV) {
		if ($smpl2 ne $smpl1 && !exists $donePairs{$smpl1}{$smpl2} && !exists $donePairs{$smpl2}{$smpl1}) {
			print STDERR "overlapping $smpl1 and $smpl2\n";
			foreach my $type (keys%{ $allCNV{$smpl1} }) {
				foreach my $chr (keys%{ $allCNV{$smpl1}{$type} }) {
					if (exists $allCNV{$smpl2}{$type}{$chr}) {
						my $S1start = \@{ $allStarts{$smpl1}{$type}{$chr} };
						my $c = 0;	#idx of @{$S1start}
						foreach my $S2start (@{ $allStarts{$smpl2}{$type}{$chr} }) {
							my $S2end = $allCNV{$smpl2}{$type}{$chr}{$S2start}{"end"};
							if ($S2start > $allCNV{$smpl1}{$type}{$chr}{${$S1start}[-1]}{"end"}) { last; }
							while ( ($c < (scalar@{$S1start}-1)) && ($S2start > $allCNV{$smpl1}{$type}{$chr}{${$S1start}[$c]}{"end"}) ) {
								$c++;
								}
							my $c2 = $c;
							my ($shortestStart,$shortestEnd);
							while ( ($c2 < scalar@{$S1start}) && ($S2end >= ${$S1start}[$c2]) ) {
								my $significantOverlap = 1;
								if ($minRatio) {
									$significantOverlap = testMinR($minRatio,$choice,${$S1start}[$c2],$allCNV{$smpl1}{$type}{$chr}{${$S1start}[$c2]}{"end"},$S2start,$S2end);
									}
								if ($significantOverlap) {
									fillOverlap1($overlap1,$type,$chr,$smpl1,$smpl2,${$S1start}[$c2],$allCNV{$smpl1}{$type}{$chr}{${$S1start}[$c2]}{"end"},$S2start,$S2end,$longest);
									}
								$c2++;
								}
							}
						}
					}
				}
			$donePairs{$smpl1}{$smpl2} = 1; $donePairs{$smpl2}{$smpl1} = 1;
			}
		}
	}

=pod
open($fh, ">", $outName."_S1.txt") or die $!;
foreach my $type (sort(keys%{$overlap1})) {
	print $fh "##\n\n$type :\n\nPos";
	foreach (@Samples) { print $fh "\t$_"; }
	print $fh "\n";
	foreach my $chr (sort(keys%{ ${$overlap1}{$type} })) {
		foreach my $start (sort{$a<=>$b}keys%{ ${$overlap1}{$type}{$chr} }) {
			foreach my $end (sort{$a<=>$b}keys%{ ${$overlap1}{$type}{$chr}{$start} }) {
				print $fh $chromName{$chr}.":$start-$end";
				foreach my $smpl (@Samples) {
					print $fh "\t";
					if (exists ${$overlap1}{$type}{$chr}{$start}{$end}{$smpl}) {
						my $txt = join(";", keys%{ ${$overlap1}{$type}{$chr}{$start}{$end}{$smpl} });
						#foreach (keys%{ ${$overlap1}{$type}{$chr}{$start}{$end}{$smpl} }) { $txt .= "$_;"; }
						#chop $txt;
						print $fh $txt;
						}
					else { print $fh "."; }
					}
				print $fh "\n";
				}
			}
		}
	}
close $fh;
=cut


## merging overlapping overlaps
my (%overlap2);

if ($transitivity) {

	foreach my $type (sort(keys%{$overlap1})) {
		foreach my $chr (sort(keys%{ ${$overlap1}{$type} })) {
			my @Starts = sort{$a<=>$b}keys%{ ${$overlap1}{$type}{$chr} };
			my $c = 0;	#idx of @Starts
			my $currStart = $Starts[$c];
			my @ends_1 = sort{$a<=>$b}keys%{ ${$overlap1}{$type}{$chr}{$Starts[$c]} };
			my $currEnd = $ends_1[-1];
			$overlap2{$type}{$chr}{$currStart}{"printStart"} = $currStart;
			$overlap2{$type}{$chr}{$currStart}{"printEnd"} = $currEnd;
			%{ $overlap2{$type}{$chr}{$currStart}{"overlaps"} } = ();
			for my $e (0..$#ends_1) {
				%{ $overlap2{$type}{$chr}{$currStart}{"overlaps"} } = ( %{ $overlap2{$type}{$chr}{$currStart}{"overlaps"} } , %{ ${$overlap1}{$type}{$chr}{$Starts[$c]}{$ends_1[$e]} } );
				}
			while ($c < $#Starts) {
				$c++;
				if ($Starts[$c] < $currEnd) {
					foreach my $end2 (keys%{ ${$overlap1}{$type}{$chr}{$Starts[$c]} }) {
						%{ $overlap2{$type}{$chr}{$currStart}{"overlaps"} } = ( %{ $overlap2{$type}{$chr}{$currStart}{"overlaps"} } , %{ ${$overlap1}{$type}{$chr}{$Starts[$c]}{$end2} } );
						}
					}
				else {
					$currStart = $Starts[$c];
					@ends_1 = sort{$a<=>$b}keys%{ ${$overlap1}{$type}{$chr}{$Starts[$c]} };
					$currEnd = $ends_1[-1];
					$overlap2{$type}{$chr}{$currStart}{"printStart"} = $currStart;
					$overlap2{$type}{$chr}{$currStart}{"printEnd"} = $currEnd;
					%{ $overlap2{$type}{$chr}{$currStart}{"overlaps"} } = ();
					for my $e (0..$#ends_1) {
						%{ $overlap2{$type}{$chr}{$currStart}{"overlaps"} } = ( %{ $overlap2{$type}{$chr}{$currStart}{"overlaps"} } , %{ ${$overlap1}{$type}{$chr}{$Starts[$c]}{$ends_1[$e]} } );
						}
					}
				}
			}
		}
	}

else {

	foreach my $type (sort(keys%{$overlap1})) {#				print STDERR "$type\n";
		foreach my $chr (sort(keys%{ ${$overlap1}{$type} })) {#			print STDERR "\t$chr\n";
			my @Starts = sort{$a<=>$b}keys%{ ${$overlap1}{$type}{$chr} };
			my $c = 1;	#idx of @Starts
			my (@ends_1,$e);
			while ($c < scalar@Starts) {					#print STDERR "\t\t$c\n";print STDERR "\t\t\t".$ends_1[0]."\n";print STDERR "\t\t\t".$Starts[$c]."\n";
				@ends_1 = sort{$a<=>$b}keys%{ ${$overlap1}{$type}{$chr}{$Starts[($c-1)]} };
				$e = 0;	#idx of @ends_1
				##reaching next overlap, filling %overlap2 with overlap1 if no overlap, but merging hash of different ends for same start
				while ( ($e < scalar@ends_1) && ($Starts[$c] > $ends_1[$e]) ) {
					if (!exists $overlap2{$type}{$chr}{$Starts[($c-1)]}) {
						$overlap2{$type}{$chr}{$Starts[($c-1)]}{"shortEnd"} = $ends_1[0];
						if ($longest) {
							$overlap2{$type}{$chr}{$Starts[($c-1)]}{"printEnd"} = $ends_1[-1];
							}
						else {
							$overlap2{$type}{$chr}{$Starts[($c-1)]}{"printEnd"} = $ends_1[0];
							}
						$overlap2{$type}{$chr}{$Starts[($c-1)]}{"printStart"} = $Starts[($c-1)];
						%{ $overlap2{$type}{$chr}{$Starts[($c-1)]}{"overlaps"} } = ();
						for my $i (0..$#ends_1) {
							%{ $overlap2{$type}{$chr}{$Starts[($c-1)]}{"overlaps"} } = ( %{ $overlap2{$type}{$chr}{$Starts[($c-1)]}{"overlaps"} } , %{ ${$overlap1}{$type}{$chr}{$Starts[($c-1)]}{$ends_1[$i]} } );
							}
						}
					$e++;
					}
				#merging overlapping intervals
				while ($e < scalar@ends_1) {
					my $c2 = $c;
					while ( ($c2 < scalar@Starts) && ($Starts[$c2] <= $ends_1[$e]) ) {
						if (!exists $overlap2{$type}{$chr}{$Starts[$c2]}) {
							## finding shortest/longest overlap end ; overlap start = $Starts[$c2]
							my @ends_2 = sort{$a<=>$b}keys%{ ${$overlap1}{$type}{$chr}{$Starts[($c2)]} };
							my $overEnd;
							if ($ends_2[0] < $ends_1[$e]) { $overEnd = $ends_2[0]; }
							else { $overEnd = $ends_1[$e]; }
							foreach my $c3 ($c..($c2-1)) {
								foreach my $end3 (sort{$a<=>$b}keys%{ ${$overlap1}{$type}{$chr}{$Starts[($c3)]} }) {
									if ($Starts[$c2] <= $end3) {
										if ($end3 < $overEnd) { $overEnd = $end3; }
										last;
										}
									}
								}
							$overlap2{$type}{$chr}{$Starts[$c2]}{"shortEnd"} = $overEnd;
							if ($longest) {
								my $longEnd;
								if ($ends_2[-1] > $ends_1[-1]) { $longEnd = $ends_2[-1]; }
								else { $longEnd = $ends_1[-1]; }
								foreach my $c3 ($c..($c2-1)) {
									my @ends_3 = sort{$a<=>$b}keys%{ ${$overlap1}{$type}{$chr}{$Starts[($c3)]} };
									if ( ($Starts[$c2] <= $ends_3[-1]) && ($ends_3[-1] > $longEnd) ) {
										$longEnd = $ends_3[-1];
										}
									}
								$overlap2{$type}{$chr}{$Starts[$c2]}{"printEnd"} = $longEnd;
								## longest overlap start
								$overlap2{$type}{$chr}{$Starts[$c2]}{"printStart"} = $Starts[($c-1)];
								}
							else {
								$overlap2{$type}{$chr}{$Starts[$c2]}{"printEnd"} = $overEnd;
								## shortest overlap start
								$overlap2{$type}{$chr}{$Starts[$c2]}{"printStart"} = $Starts[$c2];
								}
							## filling overlap2 with overlap1
							%{ $overlap2{$type}{$chr}{$Starts[$c2]}{"overlaps"} } = ();
							for my $i ($e..$#ends_1) {
								%{ $overlap2{$type}{$chr}{$Starts[$c2]}{"overlaps"} } = ( %{ $overlap2{$type}{$chr}{$Starts[$c2]}{"overlaps"} } , %{ ${$overlap1}{$type}{$chr}{$Starts[($c-1)]}{$ends_1[$i]} } );
								}
							foreach my $end2 (@ends_2) {
								%{ $overlap2{$type}{$chr}{$Starts[$c2]}{"overlaps"} } = ( %{ $overlap2{$type}{$chr}{$Starts[$c2]}{"overlaps"} } , %{ ${$overlap1}{$type}{$chr}{$Starts[$c2]}{$end2} } );
								}
							foreach my $c3 ($c..($c2-1)) {;
								my @ends_3 = sort{$a<=>$b}keys%{ ${$overlap1}{$type}{$chr}{$Starts[($c3)]} };
								for my $i (0..$#ends_3) {
									if ($Starts[$c2] <= $ends_3[$i]) {
										for my $j ($i..$#ends_3) {
											%{ $overlap2{$type}{$chr}{$Starts[$c2]}{"overlaps"} } = ( %{ $overlap2{$type}{$chr}{$Starts[$c2]}{"overlaps"} } , %{ ${$overlap1}{$type}{$chr}{$Starts[$c3]}{$ends_3[$j]} } );
											}
										last;
										}
									}
								if (exists $overlap2{$type}{$chr}{$Starts[$c3]}{"shortEnd"} && $Starts[$c2] <= $overlap2{$type}{$chr}{$Starts[$c3]}{"shortEnd"}) {
									%{ $overlap2{$type}{$chr}{$Starts[$c3]} } = ();
									}
								}
							}
						$c2++;
						}
					$e++;
					}
				$c++;
				}
			##last interval
			if (!exists $overlap2{$type}{$chr}{$Starts[-1]}) {
				my @ends_1 = sort{$a<=>$b}keys%{ ${$overlap1}{$type}{$chr}{$Starts[-1]} };
				$overlap2{$type}{$chr}{$Starts[-1]}{"shortEnd"} = $ends_1[0];
				$overlap2{$type}{$chr}{$Starts[-1]}{"printEnd"} = $ends_1[-1];
				$overlap2{$type}{$chr}{$Starts[-1]}{"printStart"} = $Starts[-1];
				%{ $overlap2{$type}{$chr}{$Starts[-1]}{"overlaps"} } = ();
				for my $i (0..$#ends_1) {
					%{ $overlap2{$type}{$chr}{$Starts[-1]}{"overlaps"} } = ( %{ $overlap2{$type}{$chr}{$Starts[-1]}{"overlaps"} } , %{ ${$overlap1}{$type}{$chr}{$Starts[-1]}{$ends_1[$i]} } );
					}
				}
			foreach my $start (keys%{ $overlap2{$type}{$chr} } ) {
				unless (%{ $overlap2{$type}{$chr}{$start} }) {
					delete $overlap2{$type}{$chr}{$start};
					}
				}
			}
		}

	}

if ($bedFile) {
	$Bed1based = readBed($bedFile);
	mergeIntervals($Bed1based);
	my $href = restrict2Intervals(\%overlap2,$Bed1based);
	%overlap2 = %{$href};
	}

## print all intersections
open($fh, ">", $outName.".txt") or die $!;
foreach my $type (sort(keys%overlap2)) {
	if ($CNVfmt eq "bed") { print $fh "\n##\n\n$type :\n\nchr\tstart\tend"; }
	else { print $fh "\n##\n\n$type :\n\nPos"; }
	foreach (@Samples) { print $fh "\t$_"; }
	print $fh "\n";
	foreach my $chr (sort(keys%{ $overlap2{$type} })) {
		foreach my $start (sort{$a<=>$b}keys%{ $overlap2{$type}{$chr} }) {
			if ($CNVfmt eq "bed") {
				print $fh $chromName{$chr}.":".($overlap2{$type}{$chr}{$start}{"printStart"} -1)."-".$overlap2{$type}{$chr}{$start}{"printEnd"};
				}
			else {
if (!exists $overlap2{$type}{$chr}{$start}{"printStart"}) { print STDERR "no start\t\t$chr:$start\n"; }
if (!exists $overlap2{$type}{$chr}{$start}{"printEnd"}) { print STDERR "no end\t\t$chr:$start\n"; }
				print $fh $chromName{$chr}.":".$overlap2{$type}{$chr}{$start}{"printStart"}."-".$overlap2{$type}{$chr}{$start}{"printEnd"};
				}
			foreach my $smpl (@Samples) {
				print $fh "\t";
				if (exists $overlap2{$type}{$chr}{$start}{"overlaps"}{$smpl}) {
					my $txt = join(";", keys%{ $overlap2{$type}{$chr}{$start}{"overlaps"}{$smpl} });
					print $fh $txt;
					}
				else { print $fh "."; }
				}
			print $fh "\n";
			}
		}
	}
close $fh;


my %overlap3;
if (%Fams) {
	foreach my $type (sort(keys%overlap2)) {
		foreach my $chr (sort(keys%{ $overlap2{$type} })) {
			foreach my $start (sort{$a<=>$b}keys%{ $overlap2{$type}{$chr} }) {
				foreach my $fam (keys%Fams) {
					my $keep = 0;
					foreach my $smpl (keys%{ $Fams{$fam}{"A"} }) {
						if (exists $overlap2{$type}{$chr}{$start}{"overlaps"}{$smpl}) { $keep++; }
						}
					if ($keep == scalar(keys%{ $Fams{$fam}{"A"} })) {
						foreach my $smpl (keys%{$Fams{$fam}{"H"}}) {
							if (exists $overlap2{$type}{$chr}{$start}{"overlaps"}{$smpl}) { $keep = 0; }
							}
						}
					else { $keep = 0; }
					if ($keep) {
						%{ $overlap3{$fam}{$type}{$chr}{$start} } = %{ $overlap2{$type}{$chr}{$start} };
						}
					}
				}
			}
		}
	}

foreach my $fam (keys%Fams) {
	if (exists $overlap3{$fam}) {
		open($fh, ">", $outName."_".$fam.".txt") or die $!;
		foreach my $type (sort(keys%{ $overlap3{$fam} })) {
			if (exists $overlap3{$fam}{$type}) {
				if ($CNVfmt eq "bed") { print $fh "\n##\n\n$type :\n\nchr\tstart\tend"; }
				else { print $fh "\n##\n\n$type :\n\nPos"; }
				foreach (@{ $Fams{$fam}{"A_all"} }) {
					if (exists $Fams{$fam}{"A"}{$_}) { print $fh "\t$_(A)"; }
					else { print $fh "\t$_(H)"; }
					}
				my $CNVothers = 0;
				foreach my $smpl (@Samples) {
					if (!exists $Fams{$fam}{"H_all"}{$smpl}) { $CNVothers++; }
					}
				print $fh "\tothers(n/$CNVothers)\n";
				foreach my $chr (sort(keys%{ $overlap3{$fam}{$type} })) {
					foreach my $start (sort{$a<=>$b}keys%{ $overlap3{$fam}{$type}{$chr} }) {
						if ($CNVfmt eq "bed") {
							print $fh $chromName{$chr}.":".($overlap2{$type}{$chr}{$start}{"printStart"} -1)."-".$overlap2{$type}{$chr}{$start}{"printEnd"};
							}
						else {
							print $fh $chromName{$chr}.":".$overlap2{$type}{$chr}{$start}{"printStart"}."-".$overlap2{$type}{$chr}{$start}{"printEnd"};
							}
						foreach my $smpl (@{ $Fams{$fam}{"A_all"} }) {
							print $fh "\t";
#							if (exists $Fams{$fam}{"A"}{$smpl}) {
							if (exists $overlap2{$type}{$chr}{$start}{"overlaps"}{$smpl}) {
								my $txt = join(";", keys%{ $overlap2{$type}{$chr}{$start}{"overlaps"}{$smpl} });
								print $fh $txt;
								}
							else { print $fh "."; }
							}
						$CNVothers = 0;
						foreach my $smpl (@Samples) {
							if (!exists $Fams{$fam}{"H_all"}{$smpl} && exists $overlap2{$type}{$chr}{$start}{"overlaps"}{$smpl}) {
								$CNVothers++;
								}
							}
						print $fh "\t$CNVothers\n";
						}
					}
				}
			else {
				print $fh "\n##\n\n$type :\tNONE\n";
				}
			}
		close $fh;
		}
	else { print STDERR "no CNV compatible with family $fam transmission hypotheses\n";}
	}




exit;




#################

#################

sub testMinR {
my ($minRatio, $choice, $S1start, $S1end, $S2start, $S2end) = @_;
my $significantOverlap = 2;
my ($minStart, $minEnd);
if ($S2start > $S1start) { $minStart = $S2start; }
else { $minStart = $S1start; }
if ($S2end < $S1end) { $minEnd = $S2end; }
else { $minEnd = $S1end; }
if ( ($minEnd-$minStart)/($S1end-$S1start) < $minRatio ) { $significantOverlap--; }
if ( ($minEnd-$minStart)/($S2end-$S2start) < $minRatio ) { $significantOverlap--; }
if ($choice eq "and" && $significantOverlap < 2) {
	$significantOverlap = 0;
	}
return($significantOverlap);
}


#################

sub fillOverlap1 {
my ($overlap1, $type, $chr, $smpl1, $smpl2, $S1start, $S1end, $S2start, $S2end, $longest) = @_;
my ($overStart,$overEnd);
if ($longest) {
	if ($S2start < $S1start) { $overStart = $S2start; }
	else { $overStart = $S1start; }
	if ($S2end > $S1end) { $overEnd = $S2end; }
	else { $overEnd = $S1end; }
	}
else {
	if ($S2start > $S1start) { $overStart = $S2start; }
	else { $overStart = $S1start; }
	if ($S2end < $S1end) { $overEnd = $S2end; }
	else { $overEnd = $S1end; }
	}
${$overlap1}{$type}{$chr}{$overStart}{$overEnd}{$smpl1}{$S1start."-".$S1end} = 1;
${$overlap1}{$type}{$chr}{$overStart}{$overEnd}{$smpl2}{$S2start."-".$S2end} = 1;
}


#################
sub readBed {

my ($bed) = @_;
my %Intervals = ();			#$covbed{$chr}{$start} = $end;	transform in hash with 1-based coord
my $nLines = 0;
open(my $fh, "<", $bed) || die "can't read file $bed ($!)\n";
print STDERR "reading $bed\n";
while (my $line = <$fh>) {
	if (($line =~ /^\w+\t\d+\t\d+/)&&($line !~ /^#/)) {
		$nLines++;
		chomp $line;
		my @tab = split(/\t/,$line);
		$tab[0] =~ s/^chr//i;
		#keep longest interval
		if ( exists $Intervals{$tab[0]}{($tab[1]+1)} ) {
			if ( $tab[2] > $Intervals{$tab[0]}{($tab[1]+1)}) {
				$Intervals{$tab[0]}{($tab[1]+1)} = $tab[2];
				}
			else { next; }
			}
		else {
			$Intervals{$tab[0]}{($tab[1]+1)} = $tab[2];
			}
		}
	}
if ($nLines) { print STDERR "\t$nLines line(s) read\n"; }
else { die "!!! err:\nno bed formatted line in $bed file";}
close $fh;
return(\%Intervals);

}


#################
sub mergeIntervals {

my ($Intervals) = @_;
print STDERR "merging bed intervals\n";
my %Interval2;
foreach my $chr (keys%{$Intervals}) {
	my @Starts = sort{$a<=>$b}(keys%{ ${$Intervals}{$chr} });
	my $start = $Starts[0];
	my $end = ${$Intervals}{$chr}{$start};
	for my $i (1..$#Starts) {
		if ($Starts[$i] <= $end) {
			if (${$Intervals}{$chr}{$Starts[$i]} > $end) {
				$end = ${$Intervals}{$chr}{$Starts[$i]};
				}
			else { next; }
			}
		else {
			$Interval2{$chr}{$start} = $end;
			$start = $Starts[$i];
			$end = ${$Intervals}{$chr}{$start};
			}
		}
	$Interval2{$chr}{$start} = $end;
	}
%{$Intervals} = %Interval2;
}


#################
sub restrict2Intervals {

my ($Overlaps,$Intervals) = @_;
print STDERR "restricting CNVs overlaps to bed intervals\n";
my %Overlap2;
my $nOread = 0;
my $nOkept = 0;
foreach my $type (keys%{$Overlaps}) {
	foreach my $chr (keys%{$Intervals}) {
		if (exists ${$Overlaps}{$type}{$chr}) {
			my @StartBed = sort{$a<=>$b}(keys%{ ${$Intervals}{$chr} });
			my @StartOverlap = sort{$a<=>$b}(keys%{ ${$Overlaps}{$type}{$chr} });
			my $b=0;	# idx in @StartBed
			my $o=0;	# idx in @StartOverlap
			my $endB;
			while ( $b < scalar@StartBed ) {	## 1st loop on bed
				$endB = ${$Intervals}{$chr}{$StartBed[$b]};
				if ( $endB < $StartOverlap[$o] ) {
					$b++;
					next;
					}
				while ( ($o < $#StartOverlap) && (${$Overlaps}{$type}{$chr}{$StartOverlap[$o]}{"shortEnd"} < $StartBed[$b]) ) {	## 2nd loop on overlap
					$o++; 
					}
				while ( ($o < scalar@StartOverlap) && ($StartOverlap[$o] <= $endB) ) {
					%{ $Overlap2{$type}{$chr}{$StartOverlap[$o]} } = %{ ${$Overlaps}{$type}{$chr}{$StartOverlap[$o]} };
					$o++;
					}
				if ( $StartBed[$b] > ${$Overlaps}{$type}{$chr}{$StartOverlap[-1]}{"shortEnd"} ) {
					last;
					}
				$b++;
				}
			}
		}
	}
#%{$Overlaps} = %Overlap2;
return(\%Overlap2);
}





