#!/usr/bin/perl

## use a tab file with a header first line, within which a "gene" field must be specified

use warnings;
use strict;
use Getopt::Long;
use File::Basename;

Getopt::Long::Configure ("bundling");
my %opt = ();
GetOptions (\%opt,
	"tab|t=s",			#tab file used to annotate
	"query|q=s",		#field name in tab file used to query input file ; def: "gene"
	"target|a=s",		#field name in input file used to look for query ; def: "ANN.*GENE"
	"fields|f=s@",		#field names in tab file used to annotate
	"col|c=s@",			#columns order in tab output
	"input|i=s",		#input file to annotate
	"output|o=s",		#output file name
	"help|h")
or die "!!! err in command line arguments\n$!\n";

unless(%opt) { usage(); }
if (exists $opt{help}) { usage(); }
sub usage {
die "usage: perl script.pl [options]
    options:
    -t / --tab [file]: tab file used for annotation
    -i / --input [file]: input file to be annotated
    -o / --output: output file name (def : STDOUT)
    -q / --query [str]: field name in tab file used to query input file ; def: \"gene\"
    -a / --target [str]: field name in input file used to look for query ; def: \"ANN.\*GENE\"
    -f / --fields [str]: field names in tab file used to annotate input file (set several times or comma separated list)(def: all)
    -c / --col [str]: columns order in tab output (set several times or comma separated list)
    -h / help: this help\n";
}

my $tabFile;
if (exists $opt{"tab"}) { 
	if ( -e $opt{"tab"} ) {
		$tabFile = $opt{"tab"};
	} else {
		die $opt{"tab"}." tab file not found\n";
	}
} else {
	die "tab file used to annotate with required (opt -t)\n";
}

my $inFile;
if (exists $opt{"input"}) { 
	if ( -e $opt{"input"} ) {
		$inFile = $opt{"input"};
	} else {
		die $opt{"input"}." input file not found\n";
	}
} else {
	die "input file to be annotated required (opt -i)\n";
}

my $outFile;
if (exists $opt{"output"}) {
	if (-d dirname($opt{"outFile"})) {
		$outFile = $opt{"output"};
	} else {
		die "!!! err: ".dirname($opt{"outFile"})." dir not found\n";
	}

}

my $query = "gene";
if (exists $opt{"query"}) {
	$query = $opt{"query"};
	$query =~ s/^#+//;
}

my $target = "ANN.*GENE";
if (exists $opt{target}) {
	$target = $opt{"target"};
	}

my @fieldOut = ();
if (exists $opt{"fields"}) {
	@fieldOut = split(/,/, join(',',@{$opt{"fields"}}));
	}

my@colOrder = ();
if (exists $opt{col}) {
	@colOrder = split(/,/, join(',',@{$opt{"col"}}));
	}


##tab file in hash:
##get idx, get values
my (%idxTab, %outTab);
my $l = 0;
open(my $fhTab, "<", $tabFile) or die "$!\n";
while (my $line = <$fhTab>) {
	chomp($line);
	my @tab = split(/\t/,$line);
	## headers
	if ($l == 0) {
		$tab[0] =~ s/^#+//;
		for (my $i=0; $i<scalar@tab; $i++) {
			$idxTab{$tab[$i]} = $i; 
		}
		unless (exists $idxTab{$query}) { die "\"$query\" field not found in $tabFile\n"; }
		if (@fieldOut) {
			foreach my$f (@fieldOut) {
				unless (exists $idxTab{$f}) { die "\"$f\" field not found in $tabFile\n"; }
			}
		} else {
			foreach (@tab) { 
				unless ($_ eq $query) { push(@fieldOut, $_); }
			}
		}
	## following lines
	} else {
		foreach my$f (@fieldOut) {
			$outTab{$tab[$idxTab{$query}]}{$f} = $tab[$idxTab{$f}];
		}
	}
	$l++;
}
close($fhTab);


my $fhOut;
if ($outFile) {
	open($fhOut, '>', $outFile) or die "$!\n";
} else {
	$fhOut = \*STDOUT;
}
my ($idxTarget, @changeCol, @whichCol, @whichField);
@changeCol = ("0");
$l = 0;
open(my $fhIn, "<", $inFile) or die "$!\n";
while (my $line = <$fhIn>) {
	chomp($line);
	my @tab = split(/\t/,$line);
	#get and print headers
	if ($l == 0) {
		$tab[0] =~ s/^#+//;
		my $ok = 0;
		for (my $i=0; $i<scalar@tab; $i++) {
			if ( $tab[$i] =~ /$target/ ) {
				$idxTarget=$i;
				$ok = 1;
				last;
			}
		}
		unless ($ok) { die "no match found for $target in headers\n"; }
		my $c = 0;
		my $newline = "#";
		foreach (@colOrder) {
			$ok = 0;
			for (my $i=0; $i<scalar@tab; $i++) {
				if ($_ eq $tab[$i]) {
					$newline .= $tab[$i]."\t";
					push(@whichCol,$i);
					$ok = 1;
					$c++;
					last;
				}
			}
			unless ($ok) {
				for (my $f=0; $f<scalar@fieldOut; $f++) {
					if ($_ eq $fieldOut[$f]) {
						$newline .= $fieldOut[$f]."\t";
						push(@changeCol,$c);
						push(@whichField,$f);
					}
				}
			}
		}
		if ($ok) { push(@changeCol,$c); }
		chop $newline;
		print $fhOut "$newline\n";
	## following lines
	} else {
		my $newline = "";
		my $c = 1;
		while ($c < scalar@changeCol) {
			for (my $i=$changeCol[$c-1]; $i<$changeCol[$c]; $i++) {
				$newline .=  $tab[$whichCol[$i]]."\t";
			}
			if ( ($c-1) < scalar@whichField ) {
				if (exists $outTab{$tab[$idxTarget]}) {
					$newline .= $outTab{$tab[$idxTarget]}{$fieldOut[$whichField[$c-1]]}."\t";
				} else {
					if ($target =~ /ANN.*FEATUREID/) {
						$tab[$idxTarget] =~ s/\.\d+$//;
						if (exists $outTab{$tab[$idxTarget]}) {
							$newline .= $outTab{$tab[$idxTarget]}{$fieldOut[$whichField[$c-1]]}."\t";
						} else {
							$newline .= "NA\t";
						}
					} else {
						$newline .= "NA\t";
					}
				}
			}
			$c++;
		}
		for (my $i=$changeCol[$c-1]; $i<scalar@whichCol; $i++) {
			$newline .=  $tab[$whichCol[$i]]."\t";
		}
		chop $newline;
		print $fhOut "$newline\n";
	}
	$l++;
}
close $fhIn; 
if ($outFile) { close $fhOut; }


exit;


#CHROM	POS	ID	REF	ALT	QUAL	DP	SAF	SAR	EFF[*].GENE	EFF[*].TRID	EFF[*].EFFECT	EFF[*].IMPACT	EFF[*].FUNCLASS	EFF[*].CODON	EFF[*].AA	EFF[*].AA_LEN	dbNSFP_Uniprot_acc	CAF	dbNSFP_1000Gp1_AF	dbNSFP_ESP6500_AA_AF	dbNSFP_ESP6500_EA_AF	dbNSFP_GERP_NR	dbNSFP_GERP_RS	dbNSFP_MutationTaster_pred	dbNSFP_SIFT_pred	dbNSFP_Polyphen2_HDIV_pred	dbNSFP_Polyphen2_HVAR_pred	dbNSFP_LRT_pred	dbNSFP_phastCons100way_vertebrate	CLNDSDBID	CLNDBN	CLNSIG	CLNORIGIN	OM	MUT	GEN[8].GT	GEN[9].GT	GEN[10].GT	GEN[11].GT	GEN[0].DP	GEN[1].DP	GEN[2].DP	GEN[3].DP	GEN[4].DP	GEN[5].DP	GEN[6].DP	GEN[7].DP	GEN[8].DP	GEN[9].DP	GEN[10].DP	GEN[11].DP	GEN[12].DP	GEN[13].DP	GEN[0].AO	GEN[1].AO	GEN[2].AO	GEN[3].AO	GEN[4].AO	GEN[5].AO	GEN[6].AO	GEN[7].AO	GEN[8].AO	GEN[9].AO	GEN[10].AO	GEN[11].AO	GEN[12].AO	GEN[13].AO
#1	949608	rs1921	G	A	17939.2	1407	433	494	ISG15	ENST00000379389	NON_SYNONYMOUS_CODING	LOW	MISSENSE	aGc/aAc	S83N	165	P05161	[0.6593,0.3407,.]	0.33974358974358976	0.414662	0.394884	3.99	-2.36	P	D	B	B	N	0.000000	.	.	.	.	.	.	1/1	1/1	0/1	1/1	51	55	53	45	27	19	127	113	193	169	123	182	140	110	21	25	0	23	0	7	0	0	192	167	100	181	140	71

