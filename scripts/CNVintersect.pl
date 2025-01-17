#!/usr/bin/perl


## print intersections of CNVs from all files provided or found


my $VERSION = "1.6";

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
	"mode|M=s",
	"info|i=s@",
	"header|H=s@",
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
    -b / --bed [file]: intervals file to restrict analysis
    -o / --outName [str]: prefix for output files (def: none in current working dir)
    -n / --name : [str]: list of sample names to select (comma separated or set several times) (def: all found),
                         names can be followed by a status: name:A/H for Affected/Healthy
    -N / --nameLs [str]: file with list of sample names to select, one per line (def: all found),
                         names can be followed by a status: name:A/H for Affected/Healthy
         --fam [file]: file with list of families, one family name per line followed by sample names and status,
                       famName<tab>sample1:A/H<tab>sample2:A/H...
    -m / --minRatio [0-1]: CNVs from 2 samples shared if overlapped by x fraction of their length (def: shared if overlapped by 1 base)
    -c / --choice [A/O]: CNVs shared if both samples are overlapped enough (And), or if only one is significantly overlapped (Or)
    -u / --union : overlap == union of intervals (def: intersection of intervals)
    -t / --transitivity : if A overlaps B and B overlaps C, then A overlaps C (def: no)
    -M / --mode [A/O]: report All CNVs, or only Overlapped ones (def: only overlaps)
    -i / --info [str][int]: add infos from input CNV lines, comma separated or set several times (def: none)
                            either col position, if input file has no header
                            or col name, if name is in a header line
	-H / --header : col_name for sup infos headers, in same order as provided in opt --info
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

my ($bedFile,$Bed1based);
if (exists $opt{bed}) {
	if (-f $opt{bed}) {
		$bedFile = $opt{bed};
	} else {
		die "!!! err: ".$opt{bed}." not found\n";
	}
}

##outName
my $outName = "CNVoverlap";
if (exists $opt{outName}) {
	if (-d $opt{outName}) {	#existing dir, no file name
		$opt{outName} =~ s/\/$//;
		$outName = $opt{outName}."/".$outName;
	} else {
		$outName = $opt{outName};
		my ($outFile,$outPath) = fileparse($outName);
		unless (-d $outPath) {
			$outPath =~ s/\/$//;
			my ($subdir,$dir) = fileparse($outPath);
			if (-d $dir) {
				mkdir $outPath;
			} else {
				die "!!! err: $outPath dir not found\n";
			}
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
		if ($status =~ /a/i) {
			$Fams{"fam1"}{"A"}{$smpl} = 1;
		} else {
			$Fams{"fam1"}{"H"}{$smpl} = 1;
		}
	} else {
		$selectedSmpl{$_} = 1;
	}
}
if (%Fams && !%selectedSmpl) { @selectedSmpl = (); }

if (exists $opt{fam}) {
	print STDERR "## reading ".$opt{fam}." family file\n";
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
					if ($status =~ /a/i) {
						$Fams{$fam}{"A"}{$smpl} = 1;
					} else {
						$Fams{$fam}{"H"}{$smpl} = 1;
					}
				} else {
					#die "!!! err: smpl format not recognized within ".$opt{fam}." file\n";
					$Fams{$fam}{"H_all"}{$_} = 1;
					push(@{ $Fams{$fam}{"A_all"} }, $_);
					$Fams{$fam}{"H"}{$_} = 1;
				}
			}
		}
	}
	close $fh;
}

my $minRatio = 0;
if (exists $opt{minRatio}) { $minRatio = $opt{minRatio}; }

my $choice = "and";
if (exists $opt{choice}) {
	if ($opt{choice} =~ /^a/i) {
		$choice = "and";
	} elsif ($opt{choice} =~ /^o/i) {
		$choice = "or";
	} else {
		die "!!! opt \"--choice\" not recognized (\"A\" or \"O\")\n";
	}
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

my $printMode = "onlyOvL";
if (exists $opt{mode}) {
	if ($opt{mode} =~ /^a/i) {
		$printMode = "allC";
	} elsif ($opt{mode} =~ /^o/i) {
		$printMode = "onlyOvL";
	} else {
		die "!!! opt \"--mode\" not recognized (\"A\" or \"O\")\n";
	}
}


my @supInfos = ();
if (exists $opt{info}) {
	@supInfos = split(/,/, join(",",@{ $opt{info} }));
}
my $infoType = "";
my %supInfoHeader = ();
if (@supInfos) {
	$infoType = "index";
	foreach (@supInfos) {
		$supInfoHeader{$_} = $_;
		if ($_ =~ /\D/) {
			$infoType = "header";
		}
	}
}
if (exists $opt{header}) {
	my $i = 0;
	foreach ( split(/,/, join(",",@{ $opt{header} })) ) {
		$supInfoHeader{$supInfos[$i]} = $_;
		$i++;
	}
	if (scalar@supInfos > $i) {
		warn "!!! warn: more sup. col. infos than header names\n";
	}
}




## reads all CNV files, filling %allCNVs{sample}{CNVtype}{chr}{start} = end
print STDERR "## reading CNV files:\n";
my (@Samples,%allCNVs,%chromName);
foreach my $file (@Files) {
	my $smplName = $File{$file}{"name"};
	$smplName =~ s/^CNV_//;
	if (!@selectedSmpl || exists $selectedSmpl{$smplName}) {
		push(@Samples,$smplName);
		if (exists $allCNVs{$smplName}) {
			die "!!! err: sample $smplName found several times\n";
		}
		print STDERR "\t$file\n";
		open($fh, "<", $file) or die $!;
		my $nLines = 0;
		my $nLinesOK = 0;
		my ($lineOk, @tab, %headerIdx, $chrName, $chr, $start, $end, $type);
		while (my $line = <$fh>) {
			$nLines++;
			chomp $line;
			@tab = split(/\t/,$line);
			if ($infoType eq "header" && $nLines == 1) {
				for (my $i = 0; $i < @tab; $i++) {
					$headerIdx{$tab[$i]} = $i;
				}
			} else {
				if ($line !~ /^#/ && $line !~ /^\s*$/) {
					$lineOk = 0;
					if ($CNVfmt eq "bed" && $line =~ m/^(\w+)\t(\d+)\t(\d+)/) {
						$lineOk = 1;
						$chrName = $1;
						$start = $2 + 1;
						$end = $3;
					} elsif ($CNVfmt eq "reg" && $tab[$idx{"pos"}] =~ m/^(\w+) *: *(\d+) *- *(\d+)/) {
						$lineOk = 1;
						$chrName = $1;
						$start = $2;
						$end = $3;
					}
					if ($lineOk) {
						$nLinesOK++;
						$chr = $chrName;
						$chr =~ s/^chr//i;
						$chromName{$chr} = $chrName;
						$type = lc($tab[$idx{"type"}]);
						$allCNVs{$smplName}{$type}{$chr}{$start}{"end"} = $end;
						foreach (@supInfos) {
							if ($infoType eq "header") {
								$allCNVs{$smplName}{$type}{$chr}{$start}{$_} = $tab[$headerIdx{$_}];
							} else {
								$allCNVs{$smplName}{$type}{$chr}{$start}{$_} = $tab[($_-1)];
							}
						}
					}
				}
			}
		}
		close $fh;
		unless ($nLinesOK) {
			print STDERR "no correctly formated line found in $file (format = $CNVfmt)\n";
		}
	}
}
if (@selectedSmpl) {
	foreach (@selectedSmpl) {
		unless (exists $allCNVs{$_}) {
			die "!!! err: sample $_ not found\n";
		}
	}
}
if (%Fams) {
	foreach my $fam (keys%Fams) {
		foreach (@{ $Fams{$fam}{"A_all"} }) {
			unless (exists $allCNVs{$_}) {
				die "!!! err: sample $_ not found\n";
			}
		}
	}
}

my $allStarts = {};
if ($bedFile) {
	$Bed1based = readBed($bedFile);
	mergeIntervals($Bed1based);
	restrict2Intervals(\%allCNVs,$allStarts,$Bed1based);
} else {
	## just order starts : so it will not be done several times
	foreach my $smpl (keys%allCNVs) {
		foreach my $type (keys%{ $allCNVs{$smpl} }) {
			foreach my $chr (keys%{ $allCNVs{$smpl}{$type} }) {
				@{ ${$allStarts}{$smpl}{$type}{$chr} } = sort{$a<=>$b}keys%{ $allCNVs{$smpl}{$type}{$chr} };
			}
		}
	}
}

## fill Overlap1 : same intervals
print STDERR "## overlapping pairs of CNV files\n";
if ($longest) { print STDERR "\tkeeping union of coordinates\n";
} else { print STDERR "\tkeeping intersection of coordinates\n";
}
my %donePairs = ();
my $overlap1 = {};	# ${$overlap1}{$type}{$chr}{$overStart}{$overEnd}{$smpl1}{$S1start."-".$S1end} = 1
foreach my $smpl1 (keys%allCNVs) {
	foreach my $smpl2 (keys%allCNVs) {
		if ($smpl2 ne $smpl1 && !exists $donePairs{$smpl1}{$smpl2} && !exists $donePairs{$smpl2}{$smpl1}) {
			print STDERR "\t$smpl1 with $smpl2\n";
			foreach my $type (keys%{ $allCNVs{$smpl1} }) {
				foreach my $chr (keys%{ $allCNVs{$smpl1}{$type} }) {
					if (exists $allCNVs{$smpl2}{$type}{$chr}) {
						my $S1start = \@{ ${$allStarts}{$smpl1}{$type}{$chr} };
						my $c = 0;	#idx of @{$S1start}
						foreach my $S2start (@{ ${$allStarts}{$smpl2}{$type}{$chr} }) {
							if ($S2start > $allCNVs{$smpl1}{$type}{$chr}{${$S1start}[-1]}{"end"}) { last; }
							while ( ($c < (scalar@{$S1start}-1)) && ($S2start > $allCNVs{$smpl1}{$type}{$chr}{${$S1start}[$c]}{"end"}) ) {
								$c++;
							}
							my $c2 = $c;
							my ($shortestStart,$shortestEnd);
							while ( ($c2 < scalar@{$S1start}) && ($allCNVs{$smpl2}{$type}{$chr}{$S2start}{"end"} >= ${$S1start}[$c2]) ) {
								my $significantOverlap = 1;
								if ($minRatio) {
									$significantOverlap = testMinR($minRatio,$choice,${$S1start}[$c2],$allCNVs{$smpl1}{$type}{$chr}{${$S1start}[$c2]}{"end"},
									                               $S2start,$allCNVs{$smpl2}{$type}{$chr}{$S2start}{"end"});
								}
								if ($significantOverlap) {
									fillOverlap1($overlap1,$type,$chr,$smpl1,$smpl2,${$S1start}[$c2],$S2start,\%allCNVs, $longest);
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


## merging overlapping overlaps
print STDERR "## Merging overlapping overlaps\n";
my %overlap2;	# $overlap2{$type}{$chr}{$start}{"overlaps"}{$smpl1}{$S1start."-".$S1end} = 1

if ($transitivity) {
	print STDERR "\tapplying transitivity accros overlaps\n";
	foreach my $type (sort(keys%{$overlap1})) {
		foreach my $chr (sort(keys%{ ${$overlap1}{$type} })) {
			# ${$overlap1}{$type}{$chr}{$overStart}{$overEnd}{$smpl1}{$S1start."-".$S1end} = 1
			my @Starts = sort{$a<=>$b}keys%{ ${$overlap1}{$type}{$chr} };
			my $c = 0;	#idx of @Starts
			my $currStart = $Starts[$c];
			my @ends_1 = sort{$a<=>$b}keys%{ ${$overlap1}{$type}{$chr}{$Starts[$c]} };
			my $currEnd = $ends_1[-1];
			$overlap2{$type}{$chr}{$currStart}{"printStart"} = $currStart;
			$overlap2{$type}{$chr}{$currStart}{"printEnd"} = $currEnd;
			%{ $overlap2{$type}{$chr}{$currStart}{"overlaps"} } = ();
			for my $e (0..$#ends_1) {
				%{ $overlap2{$type}{$chr}{$currStart}{"overlaps"} } = ( %{ $overlap2{$type}{$chr}{$currStart}{"overlaps"} },
				                                                        %{ ${$overlap1}{$type}{$chr}{$Starts[$c]}{$ends_1[$e]} } );
			}
			while ($c < $#Starts) {
				$c++;
				if ($Starts[$c] < $currEnd) {
					foreach my $end2 (keys%{ ${$overlap1}{$type}{$chr}{$Starts[$c]} }) {
						%{ $overlap2{$type}{$chr}{$currStart}{"overlaps"} } = ( %{ $overlap2{$type}{$chr}{$currStart}{"overlaps"} },
						                                                        %{ ${$overlap1}{$type}{$chr}{$Starts[$c]}{$end2} } );
					}
				} else {
					$currStart = $Starts[$c];
					@ends_1 = sort{$a<=>$b}keys%{ ${$overlap1}{$type}{$chr}{$Starts[$c]} };
					$currEnd = $ends_1[-1];
					$overlap2{$type}{$chr}{$currStart}{"printStart"} = $currStart;
					$overlap2{$type}{$chr}{$currStart}{"printEnd"} = $currEnd;
					%{ $overlap2{$type}{$chr}{$currStart}{"overlaps"} } = ();
					for my $e (0..$#ends_1) {
						%{ $overlap2{$type}{$chr}{$currStart}{"overlaps"} } = ( %{ $overlap2{$type}{$chr}{$currStart}{"overlaps"} },
						                                                        %{ ${$overlap1}{$type}{$chr}{$Starts[$c]}{$ends_1[$e]} } );
					}
				}
			}
		}
	}

} else {
	print STDERR "\twithout looking for transitivity accros overlaps\n";
	foreach my $type (sort(keys%{$overlap1})) {
		foreach my $chr (sort(keys%{ ${$overlap1}{$type} })) {
			# ${$overlap1}{$type}{$chr}{$overStart}{$overEnd}{$smpl1}{$S1start."-".$S1end} = 1
			my @Starts = sort{$a<=>$b}keys%{ ${$overlap1}{$type}{$chr} };
			my $c = 1;	#idx of @Starts
			my (@ends_1,$e);
			while ($c < scalar@Starts) {				#print STDERR "\t\t$c\n";print STDERR "\t\t\t".$ends_1[0]."\n";print STDERR "\t\t\t".$Starts[$c]."\n";
				@ends_1 = sort{$a<=>$b}keys%{ ${$overlap1}{$type}{$chr}{$Starts[($c-1)]} };
				$e = 0;	#idx of @ends_1
				##reaching next overlap, filling %overlap2 with overlap1 if no overlap, but merging hashs of different ends for same start
				while ( ($e < scalar@ends_1) && ($Starts[$c] > $ends_1[$e]) ) {
					if (!exists $overlap2{$type}{$chr}{$Starts[($c-1)]}) {
						$overlap2{$type}{$chr}{$Starts[($c-1)]}{"shortEnd"} = $ends_1[0];
						if ($longest) {
							$overlap2{$type}{$chr}{$Starts[($c-1)]}{"printEnd"} = $ends_1[-1];
						} else {
							$overlap2{$type}{$chr}{$Starts[($c-1)]}{"printEnd"} = $ends_1[0];
						}
						$overlap2{$type}{$chr}{$Starts[($c-1)]}{"printStart"} = $Starts[($c-1)];
						%{ $overlap2{$type}{$chr}{$Starts[($c-1)]}{"overlaps"} } = ();
						for my $i (0..$#ends_1) {
							%{ $overlap2{$type}{$chr}{$Starts[($c-1)]}{"overlaps"} } = ( %{ $overlap2{$type}{$chr}{$Starts[($c-1)]}{"overlaps"} },
							                                                             %{ ${$overlap1}{$type}{$chr}{$Starts[($c-1)]}{$ends_1[$i]} } );
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
							if ($ends_2[0] < $ends_1[$e]) {
								$overEnd = $ends_2[0];
							} else {
								$overEnd = $ends_1[$e];
							}
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
								if ($ends_2[-1] > $ends_1[-1]) {
									$longEnd = $ends_2[-1];
								} else {
									$longEnd = $ends_1[-1];
								}
								foreach my $c3 ($c..($c2-1)) {
									my @ends_3 = sort{$a<=>$b}keys%{ ${$overlap1}{$type}{$chr}{$Starts[($c3)]} };
									if ( ($Starts[$c2] <= $ends_3[-1]) && ($ends_3[-1] > $longEnd) ) {
										$longEnd = $ends_3[-1];
										}
									}
								$overlap2{$type}{$chr}{$Starts[$c2]}{"printEnd"} = $longEnd;
								## longest overlap start
								$overlap2{$type}{$chr}{$Starts[$c2]}{"printStart"} = $Starts[($c-1)];
							} else {
								$overlap2{$type}{$chr}{$Starts[$c2]}{"printEnd"} = $overEnd;
								## shortest overlap start
								$overlap2{$type}{$chr}{$Starts[$c2]}{"printStart"} = $Starts[$c2];
							}
							## filling overlap2 with overlap1
							%{ $overlap2{$type}{$chr}{$Starts[$c2]}{"overlaps"} } = ();
							for my $i ($e..$#ends_1) {
								%{ $overlap2{$type}{$chr}{$Starts[$c2]}{"overlaps"} } = ( %{ $overlap2{$type}{$chr}{$Starts[$c2]}{"overlaps"} },
								                                                          %{ ${$overlap1}{$type}{$chr}{$Starts[($c-1)]}{$ends_1[$i]} } );
							}
							foreach my $end2 (@ends_2) {
								%{ $overlap2{$type}{$chr}{$Starts[$c2]}{"overlaps"} } = ( %{ $overlap2{$type}{$chr}{$Starts[$c2]}{"overlaps"} },
								                                                          %{ ${$overlap1}{$type}{$chr}{$Starts[$c2]}{$end2} } );
							}
							foreach my $c3 ($c..($c2-1)) {;
								my @ends_3 = sort{$a<=>$b}keys%{ ${$overlap1}{$type}{$chr}{$Starts[($c3)]} };
								for my $i (0..$#ends_3) {
									if ($Starts[$c2] <= $ends_3[$i]) {
										for my $j ($i..$#ends_3) {
											%{ $overlap2{$type}{$chr}{$Starts[$c2]}{"overlaps"} } = ( %{ $overlap2{$type}{$chr}{$Starts[$c2]}{"overlaps"} },
											                                                          %{ ${$overlap1}{$type}{$chr}{$Starts[$c3]}{$ends_3[$j]} } );
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
					%{ $overlap2{$type}{$chr}{$Starts[-1]}{"overlaps"} } = ( %{ $overlap2{$type}{$chr}{$Starts[-1]}{"overlaps"} },
					                                                         %{ ${$overlap1}{$type}{$chr}{$Starts[-1]}{$ends_1[$i]} } );
				}
			}
			foreach my $start (keys%{ $overlap2{$type}{$chr} } ) {
				if (scalar(keys%{ $overlap2{$type}{$chr}{$start} }) == 0) {
					delete $overlap2{$type}{$chr}{$start};
				}
			}
		}
	}

}
if (%Fams) {
	foreach my $type (keys%overlap2) {
		foreach my $chr (keys%{ $overlap2{$type} }) {
			foreach my $overlapStart (keys%{ $overlap2{$type}{$chr} } ) {
				# $overlap2{$type}{$chr}{$start}{"overlaps"}{$smpl1}{$S1start} = $S1end
				# $allCNVs{$smpl1}{$type}{$chr}{$start}{"smplOv"} = $smpl2;
				foreach my $smpl (keys%{ $overlap2{$type}{$chr}{$overlapStart}{"overlaps"} }) {
					foreach my $smplStart (keys%{ $overlap2{$type}{$chr}{$overlapStart}{"overlaps"}{$smpl} }) {
						$allCNVs{$smpl}{$type}{$chr}{$smplStart}{"startOv"} = $overlapStart;
					}
				}
			}
		}
	}
}


## print all intersections
print STDERR "## printing all intersections in $outName.txt\n";
open($fh, ">", $outName.".txt") or die $!;
my $headers = "";
if ($CNVfmt eq "bed") {
	$headers .= "#chr\tstart\tend\ttype";
} else {
	$headers .= "#Pos\ttype";
}
foreach (@Samples) {
	$headers .= "\t$_";
	foreach my $info (@supInfos) {
		$headers .= "\t$_\_".$supInfoHeader{$info};
	}
}
my %mergeTypes = ();
if ($printMode eq "onlyOvL") {
	foreach my $type (sort(keys%overlap2)) {
		foreach my $chr (keys%{ $overlap2{$type} }) {
			foreach my $start (keys%{ $overlap2{$type}{$chr} }) {
				push(@{ $mergeTypes{$chr}{$start} } , makeLine($CNVfmt,$chromName{$chr},$type,\@Samples,\@supInfos,\%{ $overlap2{$type}{$chr}{$start} }));
			}
		}
	}
	if (%mergeTypes) {
		print $fh "$headers\n";
		foreach my $chr (sort(keys%mergeTypes)) {
			foreach my $start (sort{$a<=>$b}keys%{ $mergeTypes{$chr} }) {
				foreach (@{ $mergeTypes{$chr}{$start} }) {
					print $fh "$_\n";
				}
			}
		}
	}
} else {
	foreach my $smpl (@Samples) {
		foreach my $type (keys%{ $allCNVs{$smpl} }) {
			foreach my $chr (keys%{ $allCNVs{$smpl}{$type} }) {
				foreach my $start (keys%{ $allCNVs{$smpl}{$type}{$chr} }) {
					if (exists $allCNVs{$smpl}{$type}{$chr}{$start}{"smplOv"}) {
						my $overlapStart = $allCNVs{$smpl}{$type}{$chr}{$start}{"startOv"};
						if (!exists $mergeTypes{$chr}{$overlapStart}) {
							if (!exists $mergeTypes{$chr}{$overlapStart}{$type}) {
								$mergeTypes{$chr}{$overlapStart}{$type} = makeLine($CNVfmt,$chromName{$chr},$type,\@Samples,\@supInfos,
								                                                   \%{ $overlap2{$type}{$chr}{$overlapStart} });
							}
						}
					} else {
						my %miniHash;
						$miniHash{"printStart"} = $start;
						$miniHash{"printEnd"} = $allCNVs{$smpl}{$type}{$chr}{$start}{"end"};
						$miniHash{"overlaps"}{$smpl}{$start}{"end"} = $allCNVs{$smpl}{$type}{$chr}{$start}{"end"};
						foreach (@supInfos) {
							$miniHash{"overlaps"}{$smpl}{$start}{$_} = $allCNVs{$smpl}{$type}{$chr}{$start}{$_};
						}
						$mergeTypes{$chr}{$start}{$type} = makeLine($CNVfmt,$chromName{$chr},$type,\@Samples,\@supInfos,\%miniHash);
					}
				}
			}
		}
	}
	if (%mergeTypes) {
		print $fh "$headers\n";
		foreach my $chr (sort(keys%mergeTypes)) {
			foreach my $start (sort{$a<=>$b}keys%{ $mergeTypes{$chr} }) {
				foreach my $type (sort(keys%{ $mergeTypes{$chr}{$start} })) {
					print $fh $mergeTypes{$chr}{$start}{$type}."\n";
				}
			}
		}
	}
}
close $fh;


## family analysis
my %CNVbyFam;
if (%Fams) {

	print STDERR "## applying family filter\n";
	foreach my $fam (keys%Fams) {
		if ($printMode eq "onlyOvL") {
			my @affected = keys%{ $Fams{$fam}{"A"} };
			if (scalar@affected == 0) {
				foreach my $type (keys%overlap2) {
					foreach my $chr (keys%{ $overlap2{$type} }) {
						foreach my $start (keys%{ $overlap2{$type}{$chr} }) {
							my $keep = 0;
							foreach my $smpl (@{ $Fams{$fam}{"A_all"} }) {
								if (exists $overlap2{$type}{$chr}{$start}{"overlaps"}{$smpl}) {
									$keep ++;
								}
							}
							if ($keep > 1) {
								%{ $CNVbyFam{$fam}{$type}{$chr}{$start} } = %{ $overlap2{$type}{$chr}{$start} };
							}
						}
					}
				}
			} elsif (scalar@affected == 1) {
				# $allCNVs{$smpl1}{$type}{$chr}{$S1start}{"smplOv"} = 1
				my $affSmpl = $affected[0];
				foreach my $type (keys%{ $allCNVs{$affSmpl} }) {
					foreach my $chr (keys%{ $allCNVs{$affSmpl}{$type} }) {
						foreach my $start (keys%{ $allCNVs{$affSmpl}{$type}{$chr} }) {
							my $keep = 1;
							if (exists $allCNVs{$affSmpl}{$type}{$chr}{$start}{"smplOv"}) {
								foreach my $hSmpl (keys%{$Fams{$fam}{"H"}}) {
									if (exists $allCNVs{$affSmpl}{$type}{$chr}{$start}{"smplOv"}{$hSmpl}) {
										$keep = 0;
										last;
									}
								}
							}
							if ($keep) {
								if (exists $allCNVs{$affSmpl}{$type}{$chr}{$start}{"smplOv"}) {
									my $overlapStart = $allCNVs{$affSmpl}{$type}{$chr}{$start}{"startOv"};
									%{ $CNVbyFam{$fam}{$type}{$chr}{$overlapStart} } = %{ $overlap2{$type}{$chr}{$overlapStart} };
								} else {
									my $end = $allCNVs{$affSmpl}{$type}{$chr}{$start}{"end"};
									$CNVbyFam{$fam}{$type}{$chr}{$start}{"printStart"} = $start;
									$CNVbyFam{$fam}{$type}{$chr}{$start}{"printEnd"} = $end;
									$CNVbyFam{$fam}{$type}{$chr}{$start}{"overlaps"}{$affSmpl}{$start}{"end"} = $end;
									foreach (@supInfos) {
										$CNVbyFam{$fam}{$type}{$chr}{$start}{"overlaps"}{$affSmpl}{$start}{$_} = $allCNVs{$affSmpl}{$type}{$chr}{$start}{$_};
									}
								}
							}
						}
					}
				}
			} else {
				# $overlap2{$type}{$chr}{$start}{"overlaps"}{$smpl1}{$S1start} = $S1end
				foreach my $type (keys%overlap2) {
					foreach my $chr (keys%{ $overlap2{$type} }) {
						foreach my $start (keys%{ $overlap2{$type}{$chr} }) {
							my $keep = 0;
							foreach my $smpl (@affected) {
								if (exists $overlap2{$type}{$chr}{$start}{"overlaps"}{$smpl}) { $keep++; }
							}
							if ($keep == scalar@affected) {
								foreach my $smpl (keys%{$Fams{$fam}{"H"}}) {
									if (exists $overlap2{$type}{$chr}{$start}{"overlaps"}{$smpl}) { $keep = 0; }
								}
							} else {
								$keep = 0;
							}
							if ($keep) {
								%{ $CNVbyFam{$fam}{$type}{$chr}{$start} } = %{ $overlap2{$type}{$chr}{$start} };
							}
						}
					}
				}
			}
		} else {
			foreach my $smpl (@{ $Fams{$fam}{"A_all"} }) {
				foreach my $type (keys%{ $allCNVs{$smpl} }) {
					foreach my $chr (keys%{ $allCNVs{$smpl}{$type} }) {
						foreach my $start (keys%{ $allCNVs{$smpl}{$type}{$chr} }) {
							if (exists $allCNVs{$smpl}{$type}{$chr}{$start}{"smplOv"}) {
								my $overlapStart = $allCNVs{$smpl}{$type}{$chr}{$start}{"startOv"};
								%{ $CNVbyFam{$fam}{$type}{$chr}{$overlapStart} } = %{ $overlap2{$type}{$chr}{$overlapStart} };
							} else {
								my $end = $allCNVs{$smpl}{$type}{$chr}{$start}{"end"};
								$CNVbyFam{$fam}{$type}{$chr}{$start}{"printStart"} = $start;
								$CNVbyFam{$fam}{$type}{$chr}{$start}{"printEnd"} = $end;
								$CNVbyFam{$fam}{$type}{$chr}{$start}{"overlaps"}{$smpl}{$start}{"end"} = $end;
								foreach (@supInfos) {
									$CNVbyFam{$fam}{$type}{$chr}{$start}{"overlaps"}{$smpl}{$start}{$_} = $allCNVs{$smpl}{$type}{$chr}{$start}{$_};
								}
							}
						}
					}
				}
			}
		}
	}

	print STDERR "## printing intersections by family\n";
	foreach my $fam (keys%Fams) {
		print STDERR "\t$fam, -> ".$outName."_".$fam.".txt file\n";
		open($fh, ">", $outName."_".$fam.".txt") or die $!;
		if (exists $CNVbyFam{$fam}) {
			## headers :
			my $headers = "";
			if ($CNVfmt eq "bed") {
				$headers .= "#chr\tstart\tend\ttype";
			} else {
				$headers .= "#Pos\ttype";
			}
			foreach (@{ $Fams{$fam}{"A_all"} }) {
				if (exists $Fams{$fam}{"A"}{$_}) {
					$headers .= "\tboundaries_$_(A)";
				} else {
					$headers .= "\tboundaries_$_(H)";
				}
			}
			my $otherGroups = scalar(keys%Fams) - 1;
			foreach my $smpl (@Samples) {
				my $belong2Fam = 0;
				foreach my $f2 (keys%Fams) {
					if (exists $Fams{$f2}{"H_all"}{$smpl}) { 
						$belong2Fam = 1;
						last;
					}
				}
				unless ($belong2Fam) { $otherGroups++; }
			}
			$headers .= "\tin_other_groups(tot=$otherGroups)";
			foreach my $info (@supInfos) {
				foreach my $smpl (@{ $Fams{$fam}{"A_all"} }) {
					$headers .= "\t$smpl\_".$supInfoHeader{$info};
				}
			}
			## CNVs :
			my %mergeTypes = ();
			foreach my $type (sort(keys%{ $CNVbyFam{$fam} })) {
				foreach my $chr (keys%{ $CNVbyFam{$fam}{$type} }) {
					foreach my $start (keys%{ $CNVbyFam{$fam}{$type}{$chr} }) {
						my $line = "";
						## pos type
						if ($CNVfmt eq "bed") {
							$line .= $chromName{$chr}.":".($CNVbyFam{$fam}{$type}{$chr}{$start}{"printStart"} -1)."-".$CNVbyFam{$fam}{$type}{$chr}{$start}{"printEnd"};
						} else {
							$line .= $chromName{$chr}.":".$CNVbyFam{$fam}{$type}{$chr}{$start}{"printStart"}."-".$CNVbyFam{$fam}{$type}{$chr}{$start}{"printEnd"};
						}
						$line .= "\t$type";
						## included intervals
						my (%smplStarts);
						foreach my $smpl (@{ $Fams{$fam}{"A_all"} }) {
							$line .=  "\t";
							if (exists $CNVbyFam{$fam}{$type}{$chr}{$start}{"overlaps"}{$smpl}) {
								my $txt = "";
								foreach (sort{$a<=>$b}keys%{ $CNVbyFam{$fam}{$type}{$chr}{$start}{"overlaps"}{$smpl} }) {
									$txt .= $_."-".$CNVbyFam{$fam}{$type}{$chr}{$start}{"overlaps"}{$smpl}{$_}{"end"}.";";
									push(@{ $smplStarts{$smpl} }, $_);
								}
								chop $txt;
								$line .=  $txt;
							} else {
								$line .=  ".";
							}
						}
						## occurence in other families
						my $CNVothers = 0;
						my %occByFam = ();
						foreach my $smpl (@Samples) {
							if (!exists $Fams{$fam}{"H_all"}{$smpl} && exists $overlap2{$type}{$chr}{$start}{"overlaps"}{$smpl}) {
								my $belong2Fam = 0;
								foreach my $f2 (keys%Fams) {
									if (exists $Fams{$f2}{"H_all"}{$smpl}) {
										$belong2Fam = 1;
										if (!exists $occByFam{$f2}) {
											$CNVothers++;
											$occByFam{$f2} = 1;
											last;
										}
									}
								}
								unless ($belong2Fam) { $CNVothers++; }
							}
						}
						$line .= "\t$CNVothers";
						## add sup Infos
						foreach my $info (@supInfos) {
							foreach my $smpl (@{ $Fams{$fam}{"A_all"} }) {
								my $smplInfo = "";
								foreach my $smplStart (@{ $smplStarts{$smpl} }) {
									$smplInfo .= $CNVbyFam{$fam}{$type}{$chr}{$start}{"overlaps"}{$smpl}{$smplStart}{$info}.";";
								}
								chop $smplInfo;
								if ($smplInfo  eq "") {
									$smplInfo = ".";
								}
								$line .= "\t".$smplInfo;
							}
						}
						push(@{ $mergeTypes{$chr}{$start} } , $line);
					}
				}
			}
			if (%mergeTypes) {
				print $fh "$headers\n";
				foreach my $chr (sort(keys%mergeTypes)) {
					foreach my $start (sort{$a<=>$b}keys%{ $mergeTypes{$chr} }) {
						foreach (@{ $mergeTypes{$chr}{$start} }) {
							print $fh "$_\n";
						}
					}
				}
			}
		} else {
			my $line = "no CNV compatible with family $fam transmission hypotheses\n";
			print $fh $line;
			print STDERR $line;
		}
		close $fh;
	}

}



exit;




#################

#################

sub testMinR {

	my ($minRatio, $choice, $S1start, $S1end, $S2start, $S2end) = @_;
	my $significantOverlap = 2;
	my ($minStart, $minEnd);
	if ($S2start > $S1start) {
		$minStart = $S2start;
	} else {
		$minStart = $S1start;
	}
	if ($S2end < $S1end) {
		$minEnd = $S2end;
	} else { $minEnd = $S1end;
	}
	if ( ($minEnd-$minStart)/($S1end-$S1start) < $minRatio ) {
		$significantOverlap--;
	}
	if ( ($minEnd-$minStart)/($S2end-$S2start) < $minRatio ) {
		$significantOverlap--;
	}
	if ($choice eq "and" && $significantOverlap < 2) {
		$significantOverlap = 0;
		}
	return($significantOverlap);

}


#################

sub fillOverlap1 {

	my ($overlap1, $type, $chr, $smpl1, $smpl2, $S1start, $S2start, $allCNVs, $longest) = @_;
	my ($S1end, $S2end) = (${$allCNVs}{$smpl1}{$type}{$chr}{$S1start}{"end"}, ${$allCNVs}{$smpl2}{$type}{$chr}{$S2start}{"end"});
	my ($overStart,$overEnd);
	if ($longest) {
		if ($S2start < $S1start) { $overStart = $S2start; } else { $overStart = $S1start; }
		if ($S2end > $S1end) { $overEnd = $S2end; } else { $overEnd = $S1end; }
	} else {
		if ($S2start > $S1start) { $overStart = $S2start; } else { $overStart = $S1start; }
		if ($S2end < $S1end) { $overEnd = $S2end; } else { $overEnd = $S1end; }
	}
	%{ ${$overlap1}{$type}{$chr}{$overStart}{$overEnd}{$smpl1}{$S1start} } = %{ ${$allCNVs}{$smpl1}{$type}{$chr}{$S1start} };
	%{ ${$overlap1}{$type}{$chr}{$overStart}{$overEnd}{$smpl2}{$S2start} } = %{ ${$allCNVs}{$smpl2}{$type}{$chr}{$S2start} };

	${$allCNVs}{$smpl1}{$type}{$chr}{$S1start}{"smplOv"}{$smpl2} = 1;
	${$allCNVs}{$smpl2}{$type}{$chr}{$S2start}{"smplOv"}{$smpl1} = 1;

}


#################

sub readBed {

	my ($bed) = @_;
	my %Intervals = ();			#$covbed{$chr}{$start} = $end;	transform in hash with 1-based coord
	my $nLines = 0;
	open(my $fh, "<", $bed) || die "can't read file $bed ($!)\n";
	print STDERR "## reading $bed\n";
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
				} else {
					next;
				}
			} else {
				$Intervals{$tab[0]}{($tab[1]+1)} = $tab[2];
			}
		}
	}
	if ($nLines) {
		print STDERR "\t$nLines line(s) read\n";
	} else {
		die "!!! err:\nno bed formatted line in $bed file";
	}
	close $fh;
	return(\%Intervals);

}


#################

sub mergeIntervals {

	my ($Intervals) = @_;
	print STDERR "## merging bed intervals\n";
	my %Interval2;
	foreach my $chr (keys%{$Intervals}) {
		my @Starts = sort{$a<=>$b}(keys%{ ${$Intervals}{$chr} });
		my $start = $Starts[0];
		my $end = ${$Intervals}{$chr}{$start};
		for my $i (1..$#Starts) {
			if ($Starts[$i] <= $end) {
				if (${$Intervals}{$chr}{$Starts[$i]} > $end) {
					$end = ${$Intervals}{$chr}{$Starts[$i]};
				} else {
					next;
				}
			} else {
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

	my ($CNV2cut,$orderedStarts,$Intervals) = @_;
	# $allCNVs{$smpl}{$type}{$chr}{$start}{"end"} = $end
	print STDERR "## restricting CNVs to bed intervals\n";
	foreach my $smpl (keys%{$CNV2cut}) {
		foreach my $type (keys%{ ${$CNV2cut}{$smpl} }) {
			foreach my $chr (keys%{ ${$CNV2cut}{$smpl}{$type} }) {
				if (exists ${$Intervals}{$chr}) {
					my %keptCNV = ();
					my @StartBed = sort{$a<=>$b}(keys%{ ${$Intervals}{$chr} });
					my @StartCNV = sort{$a<=>$b}(keys%{ ${$CNV2cut}{$smpl}{$type}{$chr} });
					my $b = 0;	# idx in @StartBed
					my $c = 0;	# idx in @StartCNV
					my $endB;
					while ( $b < scalar@StartBed ) {
						## 1st loop on bed
						$endB = ${$Intervals}{$chr}{$StartBed[$b]};
						if ( $endB < $StartCNV[$c] ) {
							$b++;
							next;
						}
						while ( ($c < $#StartCNV) && (${$CNV2cut}{$smpl}{$type}{$chr}{$StartCNV[$c]}{"end"} < $StartBed[$b]) ) {
							## 2nd loop on overlap
							$c++; 
						}
						while ( ($c < scalar@StartCNV) && ($StartCNV[$c] <= $endB) ) {
							%{ $keptCNV{$StartCNV[$c]} } = %{ ${$CNV2cut}{$smpl}{$type}{$chr}{$StartCNV[$c]} };
							push(@{ ${$orderedStarts}{$smpl}{$type}{$chr} }, $StartCNV[$c]);
							$c++;
						}
						if ( $StartBed[$b] > ${$CNV2cut}{$smpl}{$type}{$chr}{$StartCNV[-1]}{"end"} ) {
							last;
						}
						$b++;
					}
					if (%keptCNV) {
						%{ ${$CNV2cut}{$smpl}{$type}{$chr} } = %keptCNV;
					} else {
						delete ${$CNV2cut}{$smpl}{$type}{$chr};
					}
				} else {
					delete ${$CNV2cut}{$smpl}{$type}{$chr};
				}
			}
		}
	}

}

#################

sub makeLine {
	my ($CNVfmt,$chromName,$type,$Samples,$supInfos,$overlap) = @_;
	my $line = "";
	## Pos and Type
	if ($CNVfmt eq "bed") {
		$line .= $chromName.":".(${$overlap}{"printStart"} -1)."-".${$overlap}{"printEnd"};
	} else {
		$line .= $chromName.":".${$overlap}{"printStart"}."-".${$overlap}{"printEnd"};
	}
	$line .= "\t$type";
	foreach my $smpl (@{$Samples}) {
		## CNVs, if any
		$line .= "\t";
		my @CNVstarts;
		if (exists ${$overlap}{"overlaps"}{$smpl}) {
			my $txt = "";
			foreach (sort{$a<=>$b}keys%{ ${$overlap}{"overlaps"}{$smpl} }) {
				$txt .= $_."-".${$overlap}{"overlaps"}{$smpl}{$_}{"end"}.";";
				push(@CNVstarts,$_);
			}
			chop $txt;
			$line .= $txt;
		} else {
			$line .= ".";
		}
		## add sup Infos
		foreach my $info (@{$supInfos}) {
			my $txt = "";
			foreach my $start (@CNVstarts) {
				$txt .= ${$overlap}{"overlaps"}{$smpl}{$start}{$info}.";";
			}
			chop $txt;
			if ($txt  eq "") {
				$txt = ".";
			}
			$line .= "\t$txt";
		}
	}
	return $line;
}





