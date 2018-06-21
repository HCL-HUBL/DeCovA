#!/usr/bin/perl

## -> mean alt freq and homoZ% per CNV interval, from a vcf variant file
## vcf file needs to be sorted

my $VERSION = "1.4";

use warnings;
use strict;
use Getopt::Long;

Getopt::Long::Configure ("bundling");
my %opt = ();
GetOptions (\%opt, 
	"cnv|c=s",
	"vcf|v=s",
	"smpl|s=s",
	"typeC|C=i",
	"posC|P=i",
	"help|h") or usage();
if ($opt{help}) { usage(); }
unless(%opt) { usage(); }
sub usage {
die "usage: perl /path/to/script [option] > outfile
options:
	-c / --cnv : CNVs, in bed format (chr<tab>start(0-based)<tab>end) or region format (chr:start(1-based)-end)
	-v / --vcf : to get AD infos from
	-s / --smpl : sample name (to be searched in vcf)
	-C / --typeC : index of CNV type info field (1-based)
	-P / --posC : index of position field, if not a bed file and not the first one (1-based)
	-h / --help\n";
}


my $vcf = $opt{vcf};
if ($vcf) {
	unless (-f $vcf) { die "$vcf not found\n"; }
	}
else { die "vcf file required\n"; }

my $cnvFile =  $opt{cnv};
if ($cnvFile) {
	unless (-f $cnvFile) { die "$cnvFile not found\n"; }
	}
else { die "CNV file required\n"; }

my $sample = "";
if ($opt{smpl}) { $sample = $opt{smpl}; }

my $typeCol = "";
if ($opt{typeC}) {
	$typeCol = $opt{typeC}-1;
	}

my $CNVfmt = "bed";
my $posCol = 0;
if ($opt{posC}) {
	$CNVfmt = "reg";
	$posCol = $opt{posC}-1;
	}


## read bed-CNV
print STDERR "## reading $cnvFile file\n";
my (%CNVs,@chrOrder,%chrFound,%chromName);
my $l = 0;	## line nber
open(my $fh, "<", $cnvFile) or die "cannot open $cnvFile file ($!)\n";
while (my $line = <$fh>) {
	if ($line !~ /^#/ && $line !~ /^\s*$/) {
		chomp $line;
		my @tab = split(/\t/,$line);
		my($chrName,$start,$end);
		if ($CNVfmt eq "bed" && $line =~ m/^(\w+)\t(\d+)\t(\d+)/) {
			$chrName = $1;
			$start = $2 + 1;
			$end = $3;
			}
		elsif ($CNVfmt eq "reg" && $tab[$posCol] =~ m/^(\w+) *: *(\d+) *- *(\d+)/) {
			$chrName = $1;
			$start = $2;
			$end = $3;
			}
		if (defined $chrName) {
			my $chr = $chrName;
			$chr =~ s/^chr//i;
			$chromName{$chr} = $chrName;
			unless (exists $chrFound{$chr}) {
				push(@chrOrder,$chr);
				$chrFound{$chr} = 1;
				}
			$CNVs{"coord"}{$chr}{$start}{$end} = $l;
			if ($typeCol) {
				$CNVs{"idx"}{$l}{"type"} = $tab[$typeCol];
				}
			}
		}
	$l++;
	}
close ($fh);


##get idx for sample within vcf, and variantCaller
print STDERR "## get idx of samples in $vcf\n";
my ($idx0,%idx,%sample2idx,$caller);
open($fh, "<", "$vcf") || die "can't open file $vcf ($!)\n";
while (my $line=<$fh>) {
	if($line =~ /^##/) { next; }
	elsif($line =~ /^#[^#]/) {
		chomp($line);
		my @tab = split(/\t/,$line);
		my @sampleFound = ();
		for my $i (0..$#tab) {
			$idx{$tab[$i]} = $i;
			if ($tab[$i] eq "FORMAT") { $idx0 = ($i + 1); }
			if ($idx0 && $i>=$idx0) {
				$sample2idx{$i} = $tab[$i];
				if ($sample && $sample =~ /$tab[$i]/) { push(@sampleFound,$tab[$i]); } 
				}
			}
		if ($sample) {
			if (@sampleFound) {
				if (scalar@sampleFound == 1) { $sample = $sampleFound[0]; }
				else { die "too many sample ids match $sample in $vcf\n"; }
				}
			else { die "$sample not found in $vcf\n"; }
			}
		else {
			if (keys%sample2idx > 1) { die "need to provide a sample name as there are multiple samples in $vcf\n"; }
			}
		}
	else {
		chomp($line);
		my @tab = split(/\t/,$line);
		my @Format=split(/:/,$tab[($idx0-1)]);
		##freeB	GT:DP:RO:QR:AO:QA:GL
		##gatk	GT:AD:DP:GQ:PL
		for my $i (0..$#Format) { 
			if ($Format[$i] eq "AD") { $caller="gatk"; last; }	#"AD" field in "gatk" and "vardict"
			elsif ($Format[$i] eq "AO") { $caller="freeB"; last; }
			#elsif ($Format[$i] eq "ALD") { $caller="vardict"; }
			}
		if ($caller) { last; }
		}
	}
close ($fh);

## intersect with variants from vcf
print STDERR "## intersect CNVs from $cnvFile with variants from $vcf\n";
open($fh, "<", $vcf) or die "cannot read $vcf file ($!)\n";
my $currChrom = "";
my $doneChrom = 0;
my (@Starts,%Ends,%keptVar);
my $c = 0;		## idx in @Starts
my($chr,$pos,$prevPos);
while (my $line = <$fh>) {
	if ($line !~ /^#/ && $line !~ /^\s*$/) {
#		CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2	...
		$line =~ m/^(\w+)\t(\d+)\t/;
		$chr = $1;
		$pos = $2;
		$chr =~ s/^chr//i;
		if (exists $CNVs{"coord"}{$chr}) {
			if ($chr ne $currChrom) {
				print STDERR "\t$chr\n";
				@Starts = sort{$a<=>$b}keys(%{ $CNVs{"coord"}{$chr} });
				$c = 0;
				%Ends = ();
				foreach my $start (@Starts) { 
					foreach (sort{$a<=>$b}keys(%{ $CNVs{"coord"}{$chr}{$start} })) { push(@{ $Ends{$start} }, $_); }
					}
				$doneChrom++;
				#else { print STDERR "$chr absent in $cnvFile\n"; }
				$currChrom = $chr;
				}
			else {
				if ($pos < $prevPos) { die "$vcf file not sorted? ($chr:$prevPos followed by $chr:$pos)\n";}
				}
			while ( ($c < $#Starts) && ($pos > $Ends{$Starts[$c]}[-1]) ) { $c++; }
			my $c2 = $c;			## c2 = @Starts iteration within while loop
			while ( ($c2 < scalar@Starts) && ($pos >= $Starts[$c2]) ) {
				foreach my $end (@{ $Ends{$Starts[$c2]} }) {
					if ($pos <= $end) {
						$keptVar{"$chr:$pos"} = $line;
						push(@{ $CNVs{"idx"}{$CNVs{"coord"}{$chr}{$Starts[$c2]}{$end}}{"vcf"} }, \$keptVar{"$chr:$pos"});
						}
				 	}
				$c2++;
				}
			if ( $pos > $Ends{$Starts[-1]}[-1] ) {
				if (scalar(@chrOrder) == $doneChrom) { last; }
				}
			$prevPos = $pos;
			}
		}
	}
close ($fh);
@Starts = (); %Ends = ();


## print bed-CNV with AF infos
print STDERR "## print allele frequencies associated with CNVs\n";
open($fh, "<", $cnvFile) or die "cannot open $cnvFile file ($!)\n";
$l = 0;
my $txt = "";
while (my $line = <$fh>) {
	chomp $line;
	print STDOUT "$line";
	if (exists $CNVs{"idx"}{$l}) {
		if (exists $CNVs{"idx"}{$l}{"vcf"}) {
			$txt = varFreq($idx0,\%idx,\%sample2idx,$caller,$sample,\@{ $CNVs{"idx"}{$l}{"vcf"} },$CNVs{"idx"}{$l}{"type"});
			}
	else { $txt = "\tNA"; }
		print STDOUT "$txt\n";
		}
	else 	{ print STDOUT "\n"; }
	$l++;
	}
close ($fh);



exit;


###############
sub varFreq {
my ($idx0,$idx_r,$sample2idx_r,$caller,$sample,$vcfLines_r,$CNVtype) = @_;
## initialize %counts of genotypes
my %count;
for (my $s=$idx0;$s<scalar(keys%{$idx_r});$s++) {
	$count{$s}{"het"} = 0; $count{$s}{"hom"} = 0;
	}
## extract hetC homC altF from selected vcf lines
foreach my $ref (@{$vcfLines_r}) {	
	my $line = ${$ref};
	chomp($line);
	my @fields = split(/\t/,$line);
	my @alleles = split(/,/,$fields[${$idx_r}{"ALT"}]);
	my @Format=split(/:/,$fields[($idx0-1)]);
	my (%format2idx, %idx2Format);
	for my $i (0..$#Format) {
		$format2idx{$Format[$i]} = $i;
		$idx2Format{$i} = $Format[$i];
		}
	for (my $s=$idx0;$s<@fields;$s++) {				# $s: sample idx
		if ($fields[$s] !~ /^\./) {
			my @sampleInfo = split(/:/,$fields[$s]);
			my @GT = split(/\//,$sampleInfo[$format2idx{"GT"}]);
			if ($GT[0]==0) {
				if ($GT[1]!=0) {
					$count{$s}{"het"}++;
					my $allelF = allelFreq1(\@sampleInfo,\%format2idx,$caller);
					push(@{ $count{$s}{"allelFreq"}}, $allelF);
					}
				}
			else {
				if ($GT[0]==$GT[1]) {
					$count{$s}{"hom"}++;
					}
				else {
					$count{$s}{"het"}++;
					my $allelF  = allelFreq2(\@sampleInfo,\%format2idx,$caller);
					push(@{ $count{$s}{"allelFreq"} }, $allelF);
					}
				}
			}
		}
	}

my $txt = "";
if ($sample) {
	##for case sample
	$txt .= "\t$sample:";
	my $s = ${$idx_r}{$sample};
	my ($homFrq,$allFqMoy,$allFqSD);
	my $nGen = ($count{$s}{"hom"}+$count{$s}{"het"});
	if ($nGen > 0)
		{ $homFrq = $count{$s}{"hom"} / $nGen; }
	if (exists $count{$s}{"allelFreq"}) {
		my $nAllFq = scalar@{ $count{$s}{"allelFreq"} };
		foreach (@{ $count{$s}{"allelFreq"} }) { $allFqMoy += $_; }
		$allFqMoy /= $nAllFq;
		if ($nAllFq > 1) {
			my $sqtotal;
			foreach (@{ $count{$s}{"allelFreq"} }) { $sqtotal += ($allFqMoy-$_)**2; }
			$allFqSD = ($sqtotal / ($nAllFq-1))**0.5;
			}
		else { $allFqSD = 0; }
		}
	if ($CNVtype) {
		if ($CNVtype =~ /del/i) {
			$txt .= " homF = ";
			if (defined $homFrq) { $txt .= sprintf("%.2f",$homFrq)." (n=$nGen)"; }
			else { $txt .= "na"}
			}
		elsif ($CNVtype =~ /dup/i) {
			$txt .= " allelF = ";
			if ($allFqMoy) { $txt .= sprintf("%.3f",$allFqMoy)." (SD=".sprintf("%.3f",$allFqSD).")"; }
			else { $txt .= "na"}
			}
		}
	else {
		$txt .= " homF = ";
		if (defined $homFrq) { $txt .= sprintf("%.2f",$homFrq)." (n=$nGen)"; }
		else { $txt .= "na"}
		$txt .= " ; allelF = ";
		if ($allFqMoy) { $txt .= sprintf("%.3f",$allFqMoy)." (SD=".sprintf("%.3f",$allFqSD).")"; }
		else { $txt .= "na"}
		}
	##for other samples
	$txt .= "\tother samples:";
	my@sampleList = ();
	for (my $s=$idx0;$s<scalar(keys%idx);$s++) {
		if ($s ne ${$idx_r}{$sample}) { push(@sampleList,$s); }
		}
	$homFrq=0,$allFqMoy=0,$allFqSD=0;
	my $totGT = 0;
	my @allMaf = ();
	foreach my $s (@sampleList) {
		$homFrq += $count{$s}{"hom"};
		$totGT += $count{$s}{"hom"}+$count{$s}{"het"};
		if (exists $count{$s}{"allelFreq"}) {
			foreach (@{ $count{$s}{"allelFreq"} }) { push(@allMaf, $_); }
			}
		}
	if (@allMaf) {
		my $nAllFq = scalar@allMaf;
		foreach (@allMaf) { $allFqMoy += $_; }
		$allFqMoy /= $nAllFq;
		if ($nAllFq > 1) {
			my $sqtotal;
			foreach (@allMaf) { $sqtotal += ($allFqMoy-$_)**2; }
			$allFqSD = ($sqtotal / ($nAllFq-1))**0.5;
			}
		else { $allFqSD = 0; }
		}
	if ($CNVtype) {
		if ($CNVtype =~ /del/i) {
			$txt .= " homF = ";
			if ($totGT) { $txt .= sprintf("%.2f",($homFrq/$totGT))." (n=$totGT);"; }
			else { $txt .= "na"}
			}
		elsif ($CNVtype =~ /dup/i) {
			$txt .= " allelF = ";
			if ($allFqMoy) { $txt .= sprintf("%.3f",$allFqMoy)." (SD=".sprintf("%.3f",$allFqSD).")"; }
			else { $txt .= "na"}
			}
		}
	else {
		$txt .= " homF = ";
		if ($totGT) { $txt .= sprintf("%.2f",($homFrq/$totGT))." (n=$totGT);"; }
		else { $txt .= "na"}
		$txt .= " ; allelF = ";
		if ($allFqMoy) { $txt .= sprintf("%.3f",$allFqMoy)." (SD=".sprintf("%.3f",$allFqSD).")"; }
		else { $txt .= "na"}
		}

	}
else {
	for (my $s=$idx0;$s<scalar(keys%idx);$s++) {
		$txt .= "\t".${$sample2idx_r}{$s}.":";
		my ($homFrq,$allFqMoy,$allFqSD);
		my $nGen = ($count{$s}{"hom"}+$count{$s}{"het"});
		if ($nGen > 0)
			{ $homFrq = $count{$s}{"hom"} / $nGen; }
		if (exists $count{$s}{"allelFreq"}) {
			my $AllFq = scalar@{ $count{$s}{"allelFreq"} };
			foreach (@{ $count{$s}{"allelFreq"} }) { $allFqMoy += $_; }
			$allFqMoy /= $AllFq;
			if ($AllFq > 1) {
				my $sqtotal;
				foreach (@{ $count{$s}{"allelFreq"} }) { $sqtotal += ($allFqMoy-$_)**2; }
				$allFqSD = ($sqtotal / ($AllFq-1))**0.5;
				}
			else { $allFqSD = 0; }
			}
		$txt .= " homF = ";
		if (defined $homFrq) { $txt .= sprintf("%.2f",$homFrq)." (n=$nGen)"; }
		else { $txt .= "na"}
		$txt .= " ; allelF = ";
		if ($allFqMoy) { $txt .= sprintf("%.3f",$allFqMoy)." (SD=".sprintf("%.3f",$allFqSD).")"; }
		else { $txt .= "na"}
		}
	}
return($txt);
}


######################
## get the higher allele freq among all alleles, including ref
sub allelFreq1 {
my ($sampleInfo_r,$format2idx_r,$caller) = @_;
my ($tot,@reads,$alt);
if ($caller eq "freeB") {
	@reads = split(/,/,${$sampleInfo_r}[${$format2idx_r}{"AO"}]);
	push(@reads,${$sampleInfo_r}[${$format2idx_r}{"RO"}]);
	}
else {
	@reads = split(/,/,${$sampleInfo_r}[${$format2idx_r}{"AD"}]);
	}
foreach (@reads) { $tot += $_; }
if ($tot) {
	@reads = sort{$b<=>$a}(@reads);
	$alt = $reads[0]/$tot;
	}
return($alt);
}

######################
sub allelFreq2 {
my ($sampleInfo_r,$format2idx_r,$caller) = @_;
my ($tot,@reads,$alt);
if ($caller eq "freeB") {
	@reads = split(/,/,${$sampleInfo_r}[${$format2idx_r}{"AO"}]);
	$tot = ${$sampleInfo_r}[${$format2idx_r}{"RO"}];
	foreach (@reads) { $tot += $_; }
	}
else {
	@reads = split(/,/,${$sampleInfo_r}[${$format2idx_r}{"AD"}]);
	foreach (@reads) { $tot += $_; }
	shift(@reads);
	}
if ($tot) {
	@reads = sort{$b<=>$a}(@reads);
	$alt = $reads[0]/$tot;
	}
return($alt);
}


