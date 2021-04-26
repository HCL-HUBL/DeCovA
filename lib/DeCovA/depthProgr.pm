package DeCovA::depthProgr;


use strict;
use warnings;


####################

sub samT_filter {

my ($samT,$bedFile,$threads,$gatk,$mmq,$bam,$outBam) = @_;

my $cmd = "$samT view -bh -L $bedFile";
if ($threads) { $cmd .= " -@ ".($threads-1); }
if (!$gatk && $mmq) { $cmd .= " -q $mmq"; }
$cmd .= " $bam > $outBam";
print "$cmd\n";
system "$cmd";
if ($?) { die "samtools not happy\n"; }
$cmd = "$samT index $outBam";
print "$cmd\n";
system "$cmd";

}


####################

sub bedToolsCov {

my ($intervalFile,$Files,$sName,$fName,$bam2Use,$depthFile,$tmpDir,$withChr,$chromLength,$bedT,$bedTversion,$keepCov) = @_;

my %covFiles;
my %isChr;
foreach my $f (@{$Files}) { 
	if (${$withChr}{$f}) {
		$isChr{$f} = ${$withChr}{$f};
	} else {
		$isChr{$f} = "_0Chr.bed";
	}
}
foreach my $f1 (@{$Files}) {
	$covFiles{$f1} = "$tmpDir/".${$fName}{$f1}."_cov.txt";
	bedToolCmd($intervalFile.$isChr{$f1}, ${$bam2Use}{$f1}, $covFiles{$f1}, $chromLength, $bedT, $bedTversion);
	unless (${$withChr}{$f1}) {	#test if cov == 0 at each base
		open(my $fh, "<", $covFiles{$f1}) || die "could not read ".$covFiles{$f1}." ($!)\n";
		my $ok=0;
		while (my $line = <$fh>) {
			chomp $line;
			my @tab = split(/\t/,$line);
			if ($tab[-1] != 0) {
				$ok=1;
				last;
			}
		}
		close $fh;
		unless ($ok) {
			print "try other reference genome:\n";
			foreach my $f2 (@{$Files}) { ##change all $isChr{$f} if undef ${$withChr}{$f}
				unless (${$withChr}{$f2}) { 
					if ($isChr{$f2}eq"_0Chr.bed") {
						$isChr{$f2} = "_wChr.bed";
					} else {
						$isChr{$f2} = "_0Chr.bed";
					}
				} 
			}
			bedToolCmd($intervalFile.$isChr{$f1}, ${$bam2Use}{$f1}, $covFiles{$f1}, $chromLength, $bedT, $bedTversion);
		}
		${$withChr}{$f1} = $isChr{$f1};
	}
	my $cmd = "sort -k 1,1 -k 2n,2n -k 4n,4n -T $tmpDir -o $covFiles{$f1} $covFiles{$f1}";
	print "$cmd\n";
	system "$cmd";
	if ($?) { die "sort $covFiles{$f1} failed\n"; }
}

##make unique file, as the gatk one
my $allSamplesCov = "$tmpDir/allSamples";
my $headers = "Locus\tTotal_Depth\tAverage_Depth_sample\tDepth_for_".${$sName}{${$Files}[0]};
my %smplIdx = ( ${$Files}[0]=>3 );
open(my $fhNew, ">", "$allSamplesCov.0") || die "could not create $allSamplesCov.0 ($!)\n";
open(my $fhSmpl, "<", $covFiles{${$Files}[0]}) || die "could not read ".$covFiles{${$Files}[0]}.".cov ($!)\n";
while (my $line = <$fhSmpl>) {
	chomp $line;
	my @tab = split(/\t/,$line);
	print $fhNew  $tab[0].":".($tab[1]+$tab[3])."\t.\t.\t".$tab[-1]."\n";
}
close $fhNew; close $fhSmpl;
unlink $covFiles{${$Files}[0]};
for my $i (1..$#{$Files}) {
	$headers .= "\tDepth_for_".${$sName}{${$Files}[$i]};
	$smplIdx{${$Files}[$i]} = 3+$i;
	open(my $fhOld, "<", "$allSamplesCov.".($i-1)) || die "could not write into $allSamplesCov.".($i-1)." ($!)\n";
	open($fhNew, ">", "$allSamplesCov.$i") || die "could not write into $allSamplesCov.$i ($!)\n";
	open($fhSmpl, "<", $covFiles{${$Files}[$i]}) || die "could not read ".$covFiles{${$Files}[$i]}." ($!)\n";
#	while (! eof $fhOld and ! eof $fhSmpl) {
	while (my $line1 = <$fhOld>) {
		chomp $line1;
		my $line2 = <$fhSmpl>;
		chomp $line2;
		my @tab = split(/\t/,$line2);
		print $fhNew $line1."\t".$tab[-1]."\n";
	}
	close $fhOld; close $fhNew; close $fhSmpl;
	unlink "$allSamplesCov.".($i-1);
	unlink $covFiles{${$Files}[$i]};
}
use File::Copy qw(move);
move("$allSamplesCov.".$#{$Files} , "$allSamplesCov") or die "move failed: $!\n";

##make 1 file / chr
open($fhNew, ">", $depthFile) || die "could not create $depthFile ($!)\n";
print $fhNew "$headers\n";
my $currentChr = "";
my $fhChr;
my %depthFilePerChr;
open(my $fhIn, "<", $allSamplesCov) || die "could not read $allSamplesCov ($!)\n";
while (my $line = <$fhIn>) {
	print $fhNew $line;
	$line =~ /^(\w+):\d+\t/;
	my $chr = $1; $chr =~ s/^chr//i;
	if ($chr ne $currentChr) {
		if ($currentChr ne "") { close $fhChr; }
		$currentChr = $chr;
		$depthFilePerChr{"raw"}{$chr} = "$tmpDir/$chr"."_cov.txt";
		open($fhChr, ">", $depthFilePerChr{"raw"}{$chr}) || die "could not create ".$depthFilePerChr{"raw"}{$chr}." ($!)\n";
	}
	if (defined $fhChr) {
		$line =~ s/^\w+://;
		print $fhChr $line;
	}
}
close $fhIn;
close $fhNew;
unlink $allSamplesCov;

if ($keepCov) {
	use IO::Compress::Gzip qw(gzip $GzipError) ;
	gzip $depthFile => "$depthFile.gz" or die "gzip failed: $GzipError\n";
}
unlink $depthFile;

return(\%smplIdx, \%depthFilePerChr);

}

##

sub bedToolCmd {

my ($bedFile,$bamFile,$outFile,$chromLength,$bedT,$bedTversion) = @_;
print "perform bedtools coverageBed on $bamFile\n";
my $cmd = "$bedT coverage";
if ($bedTversion == 1) {
	$cmd .= " -abam $bamFile -b $bedFile -d";
} else {
	$cmd .= " -a $bedFile -b $bamFile -d -sorted -g $chromLength";
}
$cmd .= " > $outFile"; 
print "$cmd\n";
system "$cmd";
if ($?) { die "bedtools not happy\n"; }

}


####################


#($smplIdx, $depthFilePerChr) = DeCovA::depthProgr::gatkCov($intervalName,\@Files,\%sName,\%bam2Use,$depthFile,$tmpDir,$withChr{"all"},$mmq,$mbq,$dedup,$threads,\%gatk,$genom,$jobs{"keepCov"});

sub gatkCov {

my ($intervalFile,$Files,$sName,$bam2Use,$depthFile,$tmpDir,$withChr,$mmq,$mbq,$dedup,$threads,$gatk,$genom,$keepCov) = @_;

my $allSamplesCov = "$tmpDir/allSamples";
my $cmd = ${$gatk}{"exec"};
if (${$gatk}{"version"} == 3) {
	$cmd .= " -T DepthOfCoverage -R $genom -L $intervalFile$withChr";
	foreach my $f (@{$Files}) { $cmd .= " -I ".${$bam2Use}{$f}; }
	unless ($dedup) { $cmd .= " -drf DuplicateRead"; }
	if ($threads) { $cmd .= " -nt $threads"; }
	$cmd .= " -mbq ";
	if ($mbq) {
		$cmd .= "$mbq";
	} else {
		$cmd .= "0";
	}
	$cmd .= " -mmq ";
	if ($mmq) {
		$cmd .= "$mmq";
	} else {
		$cmd .= "0";
	}
	$cmd .= " --countType COUNT_FRAGMENTS -omitLocusTable -omitSampleSummary -omitIntervals -o $allSamplesCov";
} elsif (${$gatk}{"version"} == 4) {
	$cmd .= " DepthOfCoverage -R $genom -L $intervalFile$withChr";
	foreach my $f (@{$Files}) { $cmd .= " -I ".${$bam2Use}{$f}; }
	unless ($dedup) { $cmd .= " --disable-read-filter NotDuplicateReadFilter"; }
	if ($mbq) { $cmd .= " --min-base-quality $mbq"; }	# def 0
	if ($mmq) { $cmd .= " --read-filter MappingQualityReadFilter --minimum-mapping-quality $mmq"; }
	$cmd .= " --omit-locus-table --omit-per-sample-statistics --omit-interval-statistics --output-format TABLE -O $allSamplesCov";	#--count-type COUNT_FRAGMENTS : "java.lang.UnsupportedOperationException: Fragment based counting is currently unsupported"
}
print "$cmd\n";
system "$cmd";
if ($?) { die "!!! GATK not happy\n"; }

##get smplIdx, make 1 file / chr
my (%gatkIdx,%depthFilePerChr);
open(my $fhIn, "<", $allSamplesCov) or die "could not read $allSamplesCov ($!)\n";
my $firstLine = <$fhIn>;
print $firstLine;
chomp $firstLine;
my @header_cov = split(/\t/,$firstLine);
foreach my $f (@{$Files}) {
	for (my $i=0;$i<scalar@header_cov;$i++) {
		if ($header_cov[$i] eq "Depth_for_".${$sName}{$f}) {
			$gatkIdx{$f} = $i;
			last;
		}
	}
}
my $currentChr = "";
my $fhChr;
while (my $line = <$fhIn>) {
	$line =~ /^(\w+):\d+\t/;
	my $chr = $1; $chr =~ s/^chr//i;
	if ($chr ne $currentChr) {
		if ($currentChr ne "") { close $fhChr; }
		$currentChr = $chr;
		$depthFilePerChr{"raw"}{$chr} = "$tmpDir/$chr"."_cov.txt";
		open($fhChr, ">", $depthFilePerChr{"raw"}{$chr}) || die "could not create ".$depthFilePerChr{"raw"}{$chr}." ($!)\n";
	}
	if (defined $fhChr) {
		$line =~ s/^\w+://;
		print $fhChr $line;
	}
}
close $fhChr; close $fhIn;

if ($keepCov) {
	use IO::Compress::Gzip qw(gzip $GzipError) ;
	gzip $allSamplesCov => "$depthFile.gz" or die "gzip failed: $GzipError\n";
}
unlink $allSamplesCov;

return(\%gatkIdx, \%depthFilePerChr);

}

####################

##GATK cov:
#Locus	Total_Depth	Average_Depth_sample	Depth_for_L1266a_S5
#chr9:36882001	3084	3084.00	3084
#chr9:36882002	3084	3084.00	3084

#bedtools cov:
#chr9	36882000	36882102	1	5159
#chr9	36882000	36882102	2	5159





1;

