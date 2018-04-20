package Bio::NGS::HCL::DeCovA::depthProgr;


use strict;
use warnings;


####################

sub samT_filter {

my($samT,$bedFile,$threads,$gatk,$mmq,$bam,$outBam) = @_;

my$cmd = "$samT view -bh -L $bedFile";
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

my($intervalFile,$Files,$sName,$fName,$bam2Use,$outdir,$outName,$withChr,$chromLength,$bedT,$bedTversion) = @_;

my%isChr;
foreach my$f (@{$Files}) { 
	if (${$withChr}{$f}) { $isChr{$f} = ${$withChr}{$f}; }
	else { $isChr{$f} = "_0Chr.bed"; }
	}
foreach my$f1 (@{$Files}) {
	my$outCovFile = "$outdir/".${$fName}{$f1}.".cov";
	bedToolCmd("$outdir/$intervalFile".$isChr{$f1}, ${$bam2Use}{$f1}, $outCovFile, $chromLength, $bedT, $bedTversion);
	unless (${$withChr}{$f1}) {	#test if cov == 0 at each base
		open(my$fh, "<", "$outCovFile") || die "could not read $outCovFile ($!)\n";
		my$ok=0;
		while (my$line = <$fh>) {
			chomp $line;
			my@tab = split(/\t/,$line);
			if ($tab[-1] != 0)
				{ $ok=1; last; }
			}
		close $fh;
		unless ($ok) {
			print "try other reference genome:\n";
			foreach my$f2 (@{$Files}) { ##change all $isChr{$f} if undef ${$withChr}{$f}
				unless (${$withChr}{$f2}) { 
					if ($isChr{$f2}eq"_0Chr.bed") { $isChr{$f2} = "_wChr.bed"; }
					else { $isChr{$f2} = "_0Chr.bed"; }
					} 
				}
			bedToolCmd("$outdir/$intervalFile".$isChr{$f1}, ${$bam2Use}{$f1}, $outCovFile, $chromLength, $bedT, $bedTversion);
			}
		${$withChr}{$f1} = $isChr{$f1};
		}
	}

##make unique file, as the gatk one

my$allSamplesCov = "$outdir/$outName";
my$headers = "Locus\tTotal_Depth\tAverage_Depth_sample\tDepth_for_".${$sName}{${$Files}[0]};
my%smplIdx = ( ${$Files}[0]=>3 );
open(my$fhNew, ">", "$allSamplesCov.0") || die "could not create $allSamplesCov.0 ($!)\n";
open(my$fhSmpl, "<", "$outdir/".${$fName}{${$Files}[0]}.".cov") || die "could not read $outdir/".${$fName}{${$Files}[0]}.".cov ($!)\n";
while (my $line = <$fhSmpl>) {
	chomp $line;
	my@tab = split(/\t/,$line);
	print $fhNew  $tab[0].":".($tab[1]+$tab[3])."\t.\t.\t".$tab[-1]."\n";
	}
close $fhNew; close $fhSmpl;
for my$i (1..$#{$Files}) {
	$headers .= "\tDepth_for_".${$sName}{${$Files}[$i]};
	$smplIdx{${$Files}[$i]} = 3+$i;
	my$smplCovFile = "$outdir/".${$fName}{${$Files}[$i]}.".cov";
	open(my$fhOld, "<", "$allSamplesCov.".($i-1)) || die "could not write into $allSamplesCov.".($i-1)." ($!)\n";
	open($fhNew, ">", "$allSamplesCov.$i") || die "could not write into $allSamplesCov.$i ($!)\n";
	open($fhSmpl, "<", "$smplCovFile") || die "could not read $smplCovFile ($!)\n";
#	while (! eof $fhOld and ! eof $fhSmpl) {
	while (my$line1 = <$fhOld>) {
		chomp $line1;
		my$line2 = <$fhSmpl>;
		chomp $line2;
		my@tab = split(/\t/,$line2);
		print $fhNew $line1."\t".$tab[-1]."\n";
		}
	close $fhOld; close $fhNew; close $fhSmpl;
	unlink "$allSamplesCov.".($i-1);
	}
use File::Copy qw/move/;
move("$allSamplesCov.".$#{$Files} , "$allSamplesCov") or die "move failed: $!\n";

my$cmd = "sort -k 1,1 -k 2n,2n -k 4n,4n -T $outdir -o $allSamplesCov $allSamplesCov";
print "$cmd\n";
system "$cmd";

open($fhNew, ">", "$outdir/$outName.headers.txt") || die "could not create $outdir/$outName.headers.txt ($!)\n";
print $fhNew "$headers\n";
close $fhNew;

return(\%smplIdx);

}

##

sub bedToolCmd {

my($bedFile,$bamFile,$outFile,$chromLength,$bedT,$bedTversion) = @_;
print "perform bedtools coverageBed on $bamFile\n";
my$cmd = "$bedT coverage";
if ($bedTversion == 1) { $cmd .= " -abam $bamFile -b $bedFile -d"; }
else	{ $cmd .= " -a $bedFile -b $bamFile -d -sorted -g $chromLength"; }
$cmd .= " > $outFile"; 
print "$cmd\n";
system "$cmd";
if ($?) { die "bedtools not happy\n"; }

}


####################

#%gatkIdx = gatkCov($intervalName,\%path,$extenS,\@Files,\%fName,$outdir,"all.cov",$withChr{"all"},$mmq,$mbq,$dedup,$gatk,$picard,$genom);
sub gatkCov {

my($intervalFile,$Files,$sName,$bam2Use,$outdir,$outName,$withChr,$mmq,$mbq,$dedup,$threads,$gatk,$genom) = @_;

my$cmd = "$gatk -T DepthOfCoverage -R $genom -L $outdir/$intervalFile$withChr";
foreach my$f (@{$Files}) { $cmd .= " -I ".${$bam2Use}{$f}; }
unless ($dedup) { $cmd .= " -drf DuplicateRead"; }
if ($threads) { $cmd .= " -nt $threads"; }
$cmd .= " -mbq ";
if ($mbq) { $cmd .= "$mbq"; }
else { $cmd .= "0"; }
$cmd .= " -mmq ";
if ($mmq) { $cmd .= "$mmq"; }
else { $cmd .= "0"; }
$cmd .= " --countType COUNT_FRAGMENTS -omitLocusTable -omitSampleSummary -omitIntervals -o $outdir/$outName";
print "$cmd\n";
system "$cmd";
if ($?) { die "GATK not happy\n"; }

my%gatkIdx;
open(my$fh, "<","$outdir/$outName") or die "could not read $outdir/$outName ($?)\n";
my$firstLine = <$fh>;
print $firstLine;
chomp $firstLine;
my@header_cov = split(/\t/,$firstLine);
foreach my$f (@{$Files}) {
	for (my$i=0;$i<scalar@header_cov;$i++) {
		if ($header_cov[$i] eq "Depth_for_".${$sName}{$f})
			{ $gatkIdx{$f} = $i; last; }
		}
	}
close $fh;

$cmd = "head -n 1 $outdir/$outName > $outdir/$outName.headers.txt";
print "$cmd\n";
system "$cmd";
$cmd = "tail -n +2 $outdir/$outName > $outdir/tmpTail; mv $outdir/tmpTail $outdir/$outName";
print "$cmd\n";
system "$cmd";
return(\%gatkIdx);

}

####################

#GATK cov:
#Locus	Total_Depth	Average_Depth_sample	Depth_for_L1266a_S5
#chr9:36882001	3084	3084.00	3084
#chr9:36882002	3084	3084.00	3084

#bedtools cov:
#chr9	36882000	36882102	1	5159
#chr9	36882000	36882102	2	5159





1;

