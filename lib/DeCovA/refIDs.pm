package DeCovA::refIDs;


use strict;
use warnings;


##############################

#print lines of RefSeq matching each $id of @IDs

sub UCSC2Ids {

my($refFile,$IDs,$chromName,$wNonCod)=@_;

my $fhIn;
if ($refFile =~ /.gz$/) {
	use IO::Zlib;
	$fhIn = new IO::Zlib;
	$fhIn->open($refFile, "rb")
	}
else {
	open($fhIn, "<", $refFile) or die "could not read $refFile ($!)\n";
	}

print "\treading $refFile\n";
my $allRefs = UCSCfile2hash($fhIn,$chromName);
close($fhIn);

## selecting IDs from hash:
my %targets;
foreach my $id (@{$IDs}) {
	$id = uc($id);
	if (exists ${$allRefs}{"transcript"}{$id}) { %{ $targets{$id} } = %{ ${$allRefs}{"transcript"}{$id} }; }
	elsif (exists ${$allRefs}{"gene"}{$id}) {
		foreach my $nm (keys%{ ${$allRefs}{"gene"}{$id} }) {
			if ($wNonCod || exists ${$allRefs}{"transcript"}{$nm}{"CDS_start"}) {
				%{ $targets{$nm} } = %{ ${$allRefs}{"transcript"}{$nm} };
				}
			}
		}
	elsif (exists ${$allRefs}{"Id_noExt"}{$id}) {
		foreach my $nm (keys%{ ${$allRefs}{"Id_noExt"}{$id} }) {
			if ($wNonCod || exists ${$allRefs}{"transcript"}{$nm}{"CDS_start"}) {
				%{ $targets{$nm} } = %{ ${$allRefs}{"transcript"}{$nm} };
				}
			}
		}
	else { die "$id not found in $refFile file\n"; }
	}

return(\%targets);

}

####

sub UCSCfile2hash {

my($fhIn,$chromName) = @_;

my(%allRefs,%idx);
$idx{"ID"} = 1;
$idx{"chr"} = 2;
$idx{"strand"} = 3;
$idx{"cdsStart"} = 6;
$idx{"cdsEnd"} = 7;
$idx{"exStart"} = 9;
$idx{"exEnd"} = 10;
$idx{"gene"} = 12;

while (my $line = <$fhIn>) {
	unless($line =~ /^#/) {
		chomp $line;
		my @tab = split(/\t/,$line);
		my $transcript_id = uc($tab[$idx{"ID"}]);
		my $gene_name = uc($tab[$idx{"gene"}]);
		my $Id_noExt = $transcript_id;
		$Id_noExt =~ s/\.(\d+)$//;
		my $chr = $tab[$idx{"chr"}];
		$chr =~ s/^chr//i;
		## for synonyms: 
		if (exists $allRefs{"transcript"}{$transcript_id}) {
			$allRefs{"transcript"}{$transcript_id}{"nbr"}++;
			$transcript_id = $transcript_id."-".$allRefs{"transcript"}{$transcript_id}{"nbr"};
			}
		${$chromName}{"refId"}{$chr} = $tab[$idx{"chr"}];
		$allRefs{"transcript"}{$transcript_id}{"chr"} = $chr;
		$allRefs{"transcript"}{$transcript_id}{"strand"} = $tab[$idx{"strand"}];
		if ($tab[$idx{"cdsStart"}] != $tab[$idx{"cdsEnd"}]) {
			$allRefs{"transcript"}{$transcript_id}{"CDS_start"} = $tab[$idx{"cdsStart"}] + 1;
			$allRefs{"transcript"}{$transcript_id}{"CDS_end"} = $tab[$idx{"cdsEnd"}];
			}
		my @starts = split(/,/, $tab[$idx{"exStart"}]);
		for my $i (0..$#starts) { $allRefs{"transcript"}{$transcript_id}{"ex_starts"}[$i] = $starts[$i] + 1; }		# -> 1-based
		@{ $allRefs{"transcript"}{$transcript_id}{"ex_ends"} } = split(/,/, $tab[$idx{"exEnd"}]);
		$allRefs{"transcript"}{$transcript_id}{"gene"} = $gene_name;
		}
	}

return(\%allRefs);

}


##############################

sub GTF2Ids {

my($refFile,$IDs,$chromName,$wNonCod)=@_;

## reading $refFile, in hash %allRefs:

# GTF (General Transfer Format) is identical to GFF version 2
#tab-separated fields:
#1	seqname - name of the chromosome or scaffold
#2	source - name of the program that generated this feature, or the data source (database or project name)
#3	feature - feature type name, e.g. Gene, Variation, Similarity
#4	start - Start position of the feature, with sequence numbering starting at 1.
#5	end - End position of the feature, with sequence numbering starting at 1.
#6	score - A floating point value.
#7	strand - defined as + (forward) or - (reverse).
#8	frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
#9	attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.

#features:
#	CDS
#	Selenocysteine
#	exon
#	five_prime_utr
#	gene
#	start_codon
#	stop_codon
#	three_prime_utr
#	transcript

#transcript_biotypes:
#	unprocessed_pseudogene
#	IG_J_pseudogene
#	miRNA
#	pseudogene
#	sense_intronic
#	lincRNA
#	TR_V_pseudogene
#	protein_coding
#	transcribed_processed_pseudogene
#	TR_J_pseudogene
#	processed_transcript
#	Mt_tRNA
#	rRNA
#	transcribed_unprocessed_pseudogene
#	non_stop_decay
#	IG_V_pseudogene
#	polymorphic_pseudogene
#	nonsense_mediated_decay
#	IG_C_pseudogene
#	TR_C_gene
#	misc_RNA
#	TR_D_gene
#	processed_pseudogene
#	IG_C_gene
#	Mt_rRNA
#	snRNA
#	TR_J_gene
#	unitary_pseudogene
#	translated_processed_pseudogene
#	retained_intron
#	TR_V_gene
#	3prime_overlapping_ncrna
#	snoRNA
#	IG_D_gene
#	IG_J_gene
#	sense_overlapping
#	antisense
#	IG_V_gene


#ex:

#1       ensembl_havana  gene    955503  991496  .       +       .       gene_id "ENSG00000188157"; gene_version "9"; gene_name "AGRN"; gene_source "ensembl_havana"; gene_biotype "protein_coding";
#1       ensembl_havana  transcript      955503  991496  .       +       .       gene_id "ENSG00000188157"; gene_version "9"; transcript_id "ENST00000379370"; transcript_version "2"; gene_name "AGRN"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "AGRN-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS30551"; havana_transcript "OTTHUMT00000097990"; havana_transcript_version "2"; tag "basic";
#1       ensembl_havana  exon    955503  955753  .       +       .       gene_id "ENSG00000188157"; gene_version "9"; transcript_id "ENST00000379370"; transcript_version "2"; exon_number "1"; gene_name "AGRN"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "AGRN-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS30551"; havana_transcript "OTTHUMT00000097990"; havana_transcript_version "2"; exon_id "ENSE00002277560"; exon_version "1"; tag "basic";
#1       ensembl_havana  CDS     955553  955753  .       +       0       gene_id "ENSG00000188157"; gene_version "9"; transcript_id "ENST00000379370"; transcript_version "2"; exon_number "1"; gene_name "AGRN"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "AGRN-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS30551"; havana_transcript "OTTHUMT00000097990"; havana_transcript_version "2"; protein_id "ENSP00000368678"; protein_version "2"; tag "basic";
#1       ensembl_havana  start_codon     955553  955555  .       +       0       gene_id "ENSG00000188157"; gene_version "9"; transcript_id "ENST00000379370"; transcript_version "2"; exon_number "1"; gene_name "AGRN"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "AGRN-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS30551"; havana_transcript "OTTHUMT00000097990"; havana_transcript_version "2"; tag "basic";
#1       ensembl_havana  exon    957581  957842  .       +       .       gene_id "ENSG00000188157"; gene_version "9"; transcript_id "ENST00000379370"; transcript_version "2"; exon_number "2"; gene_name "AGRN"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "AGRN-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS30551"; havana_transcript "OTTHUMT00000097990"; havana_transcript_version "2"; exon_id "ENSE00001673004"; exon_version "1"; tag "basic";

#...

#1       ensembl_havana  exon    990204  991496  .       +       .       gene_id "ENSG00000188157"; gene_version "9"; transcript_id "ENST00000379370"; transcript_version "2"; exon_number "36"; gene_name "AGRN"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "AGRN-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS30551"; havana_transcript "OTTHUMT00000097990"; havana_transcript_version "2"; exon_id "ENSE00001883271"; exon_version "1"; tag "basic";
#1       ensembl_havana  CDS     990204  990358  .       +       2       gene_id "ENSG00000188157"; gene_version "9"; transcript_id "ENST00000379370"; transcript_version "2"; exon_number "36"; gene_name "AGRN"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "AGRN-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS30551"; havana_transcript "OTTHUMT00000097990"; havana_transcript_version "2"; protein_id "ENSP00000368678"; protein_version "2"; tag "basic";
#1       ensembl_havana  stop_codon      990359  990361  .       +       0       gene_id "ENSG00000188157"; gene_version "9"; transcript_id "ENST00000379370"; transcript_version "2"; exon_number "36"; gene_name "AGRN"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "AGRN-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS30551"; havana_transcript "OTTHUMT00000097990"; havana_transcript_version "2"; tag "basic";
#1       ensembl_havana  five_prime_utr  955503  955552  .       +       .       gene_id "ENSG00000188157"; gene_version "9"; transcript_id "ENST00000379370"; transcript_version "2"; gene_name "AGRN"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "AGRN-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS30551"; havana_transcript "OTTHUMT00000097990"; havana_transcript_version "2"; tag "basic";
#1       ensembl_havana  three_prime_utr 990362  991496  .       +       .       gene_id "ENSG00000188157"; gene_version "9"; transcript_id "ENST00000379370"; transcript_version "2"; gene_name "AGRN"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "AGRN-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS30551"; havana_transcript "OTTHUMT00000097990"; havana_transcript_version "2"; tag "basic";

print "\treading $refFile\n";
my $fhIn;
if ($refFile =~ /.gz$/) {
	use IO::Zlib;
	$fhIn = new IO::Zlib;
	$fhIn->open($refFile, "rb")
	}
else {
	open($fhIn, "<", $refFile) or die "could not read $refFile ($!)\n";
	}
my $allRefs = GTFfile2hash($fhIn,$chromName);
close($fhIn);

## selecting IDs from hash:
my %targets;
foreach my $id (@{$IDs}) {
	$id = uc($id);
	if (exists ${$allRefs}{"transcript"}{$id}) { %{ $targets{$id} } = %{ ${$allRefs}{"transcript"}{$id} }; }
	elsif (exists ${$allRefs}{"gene"}{$id}) {
		foreach my $nm (keys%{ ${$allRefs}{"gene"}{$id} }) {
			if ($wNonCod || exists ${$allRefs}{"transcript"}{$nm}{"CDS_starts"}) {
				%{ $targets{$nm} } = %{ ${$allRefs}{"transcript"}{$nm} };
				}
			}
		}
	elsif (exists ${$allRefs}{"Id_noExt"}{$id}) {
		foreach my $nm (keys%{ ${$allRefs}{"Id_noExt"}{$id} }) {
			if ($wNonCod || exists ${$allRefs}{"transcript"}{$nm}{"CDS_starts"}) {
				%{ $targets{$nm} } = %{ ${$allRefs}{"transcript"}{$nm} };
				}
			}
		}
	else { die "$id not found in $refFile file\n"; }
	}
foreach my $id (keys%targets) {
	delete $targets{$id}{"ex_starts"};
	@{ $targets{$id}{"ex_starts"} } = sort{$a<=>$b}keys%{ ${$allRefs}{"transcript"}{$id}{"ex_starts"} };
	delete $targets{$id}{"ex_ends"};
	@{ $targets{$id}{"ex_ends"} } = sort{$a<=>$b}keys%{ ${$allRefs}{"transcript"}{$id}{"ex_ends"} };
	if (exists $targets{$id}{"CDS_starts"}) {
		delete $targets{$id}{"CDS_starts"};
		my @pos = sort{$a<=>$b}keys%{ ${$allRefs}{"transcript"}{$id}{"CDS_starts"} };
		$targets{$id}{"CDS_start"} = $pos[0];
		delete $targets{$id}{"CDS_ends"};
		@pos = sort{$a<=>$b}keys%{ ${$allRefs}{"transcript"}{$id}{"CDS_ends"} };
		$targets{$id}{"CDS_end"} = $pos[-1];
		}
	}

return(\%targets);

}

####

sub GTFfile2hash {

my($fhIn,$chromName) = @_;

my(%allRefs,%idx);
$idx{"chr"} = 0;
$idx{"type"} = 2;
$idx{"start"} = 3;
$idx{"end"} = 4;
$idx{"strand"} = 6;
$idx{"info"} = 8;
while (my $line = <$fhIn>) {
	unless($line =~ /^#/) {
		chomp $line;
		my @tab = split(/\t/,$line);
		my $chr = $tab[$idx{"chr"}];
		$chr =~ s/^chr//i;
		${$chromName}{"refId"}{$chr} = $tab[$idx{"chr"}];
		my($gene_id,$transcript_id,$gene_name,$biotype);
		if ($tab[$idx{"type"}] eq "transcript") {
			my @infos = split(/;/, $tab[$idx{"info"}]);
#GTF transcript attributes:
#gene_id "ENSG00000188157"; gene_version "9"; transcript_id "ENST00000379370"; transcript_version "2"; gene_name "AGRN"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "AGRN-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS30551"; havana_transcript "OTTHUMT00000097990"; havana_transcript_version "2"; tag "basic";
			foreach (@infos) {
				$_ =~ s/^\s+//;
				if ($_ =~ /^(gene_id )/) { $gene_id = $_; $gene_id =~ s/^$1//; $gene_id =~ s/"//g; $gene_id = uc($gene_id); }
				elsif ($_ =~ /^(gene_name )/) { $gene_name = $_; $gene_name =~ s/^$1//; $gene_name =~ s/"//g; $gene_name = uc($gene_name); }
				elsif ($_ =~ /^(transcript_id )/) { $transcript_id = $_; $transcript_id =~ s/^$1//; $transcript_id =~ s/"//g; $transcript_id = uc($transcript_id); }
				elsif ($_ =~ /^(transcript_biotype )/) { $biotype = $_; $biotype =~ s/^$1//; $biotype =~ s/"//g; }
				}
			#if (exists $allRefs{"transcript"}{$transcript_id}) {
			#	print "$transcript_id...$chr...$gene_name\n";
			#	}
			$allRefs{"transcript"}{$transcript_id}{"gene"} = $gene_name;
			$allRefs{"transcript"}{$transcript_id}{"biotype"} = $biotype;
			$allRefs{"transcript"}{$transcript_id}{"chr"} = $chr;
			$allRefs{"transcript"}{$transcript_id}{"strand"} = $tab[$idx{"strand"}];
			my $Id_noExt = $transcript_id;
			$Id_noExt =~ s/\.(\d+)$//;
			$allRefs{"Id_noExt"}{$Id_noExt}{$transcript_id} = 1;
			$allRefs{"gene"}{$gene_name}{$transcript_id} = 1;
			$allRefs{"gene"}{$gene_id}{$transcript_id} = 1;
			}
		elsif ($tab[$idx{"type"}] =~ /^(exon|CDS)$/) {
			my @infos = split(/;/, $tab[$idx{"info"}]);
			foreach (@infos) {
				$_ =~ s/^\s+//;
				if ($_ =~ /^(transcript_id )/) { $transcript_id = $_; $transcript_id =~ s/^$1//; $transcript_id =~ s/"//g; $transcript_id = uc($transcript_id); last; }
				}
#GTF:	chr	source	feature	start	end	score	strand	frame	attribute
#Ref:	bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames
			if ($tab[$idx{"type"}] eq "exon") {
				$allRefs{"transcript"}{$transcript_id}{"ex_starts"}{$tab[$idx{"start"}]} = 1;
				$allRefs{"transcript"}{$transcript_id}{"ex_ends"}{$tab[$idx{"end"}]} = 1;
				}
			else {
				$allRefs{"transcript"}{$transcript_id}{"CDS_starts"}{$tab[$idx{"start"}]} = 1;
				$allRefs{"transcript"}{$transcript_id}{"CDS_ends"}{$tab[$idx{"end"}]} = 1;
				}
			}

		}
	}

return(\%allRefs);

}


##############################

sub GFF2Ids {

my($refFile,$IDs,$chromName,$wNonCod)=@_;

# GFF (General Feature Format)  version 3
#tab-separated fields:
#1	seqid - name of the chromosome or scaffold
#2	source - name of the program that generated this feature, or the data source (database or project name)
#3	type - type of feature. Must be a term or accession from the SOFA sequence ontology
#4	start - Start position of the feature, with sequence numbering starting at 1.
#5	end - End position of the feature, with sequence numbering starting at 1.
#6	score - A floating point value.
#7	strand - defined as + (forward) or - (reverse).
#8	phase - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
#9	attributes - A semicolon-separated list of tag-value pairs, providing additional information about each feature. Some of these tags are predefined, e.g. ID, Name, Alias, Parent

#types:
#	CDS
#	C_gene_segment
#	D_gene_segment
#	J_gene_segment
#	V_gene_segment
#	biological_region
#	chromosome
#	exon
#	five_prime_UTR
#	gene
#	lnc_RNA
#	mRNA
#	miRNA
#	ncRNA
#	ncRNA_gene
#	pseudogene
#	pseudogenic_transcript
#	rRNA
#	scRNA
#	snRNA
#	snoRNA
#	supercontig
#	tRNA
#	three_prime_UTR
#	transcript
#	vaultRNA_primary_transcript

#biotypes:
#	processed_transcript
#	Mt_rRNA
#	rRNA
#	retained_intron
#	polymorphic_pseudogene
#	processed_pseudogene
#	IG_C_pseudogene
#	snRNA
#	transcribed_processed_pseudogene
#	antisense
#	sense_intronic
#	TR_J_gene
#	TR_D_gene
#	TR_C_gene
#	3prime_overlapping_ncrna
#	non_stop_decay
#	protein_coding
#	IG_V_pseudogene
#	unitary_pseudogene
#	Mt_tRNA
#	pseudogene
#	snoRNA
#	translated_processed_pseudogene
#	IG_C_gene
#	misc_RNA
#	sense_overlapping
#	unprocessed_pseudogene
#	transcribed_unprocessed_pseudogene
#	TR_J_pseudogene
#	IG_J_gene
#	nonsense_mediated_decay
#	TR_V_gene
#	lincRNA
#	miRNA
#	IG_V_gene
#	IG_J_pseudogene
#	IG_D_gene
#	TR_V_pseudogene

#ex:
#1       ensembl_havana  gene    955503  991496  .       +       .       ID=gene:ENSG00000188157;Name=AGRN;biotype=protein_coding;description=agrin [Source:HGNC Symbol%3BAcc:329];gene_id=ENSG00000188157;logic_name=ensembl_havana_gene;version=9
#1       ensembl_havana  mRNA    955503  991496  .       +       .       ID=transcript:ENST00000379370;Parent=gene:ENSG00000188157;Name=AGRN-001;biotype=protein_coding;ccdsid=CCDS30551.1;havana_transcript=OTTHUMT00000097990;havana_version=2;tag=basic;transcript_id=ENST00000379370;version=2
#1       ensembl_havana  five_prime_UTR  955503  955552  .       +       .       Parent=transcript:ENST00000379370
#1       ensembl_havana  exon    955503  955753  .       +       .       Parent=transcript:ENST00000379370;Name=ENSE00002277560;constitutive=0;ensembl_end_phase=0;ensembl_phase=-1;exon_id=ENSE00002277560;rank=1;version=1
#1       ensembl_havana  CDS     955553  955753  .       +       0       ID=CDS:ENSP00000368678;Parent=transcript:ENST00000379370;protein_id=ENSP00000368678
#1       ensembl_havana  exon    957581  957842  .       +       .       Parent=transcript:ENST00000379370;Name=ENSE00001673004;constitutive=0;ensembl_end_phase=1;ensembl_phase=0;exon_id=ENSE00001673004;rank=2;version=1
#1       ensembl_havana  CDS     957581  957842  .       +       0       ID=CDS:ENSP00000368678;Parent=transcript:ENST00000379370;protein_id=ENSP00000368678

#...

#1       ensembl_havana  exon    989828  989931  .       +       .       Parent=transcript:ENST00000379370;Name=ENSE00003619422;constitutive=0;ensembl_end_phase=1;ensembl_phase=2;exon_id=ENSE00003619422;rank=35;version=1
#1       ensembl_havana  CDS     989828  989931  .       +       1       ID=CDS:ENSP00000368678;Parent=transcript:ENST00000379370;protein_id=ENSP00000368678
#1       ensembl_havana  CDS     990204  990361  .       +       2       ID=CDS:ENSP00000368678;Parent=transcript:ENST00000379370;protein_id=ENSP00000368678
#1       ensembl_havana  exon    990204  991496  .       +       .       Parent=transcript:ENST00000379370;Name=ENSE00001883271;constitutive=0;ensembl_end_phase=-1;ensembl_phase=1;exon_id=ENSE00001883271;rank=36;version=1
#1       ensembl_havana  three_prime_UTR 990362  991496  .       +       .       Parent=transcript:ENST00000379370

#1       havana  processed_transcript    969486  976105  .       +       .       ID=transcript:ENST00000477585;Parent=gene:ENSG00000188157;Name=AGRN-002;biotype=processed_transcript;havana_transcript=OTTHUMT00000097991;havana_version=1;transcript_id=ENST00000477585;version=1

#12	ensembl	miRNA_gene	105721684	105721770	.	-	.	ID=gene:ENSG00000264749;Name=AC011313.1;biotype=miRNA;gene_id=ENSG00000264749;logic_name=ncrna;version=1
#12	ensembl	miRNA	105721684	105721770	.	-	.	ID=transcript:ENST00000584487;Parent=gene:ENSG00000264749;Name=AC011313.1-201;biotype=miRNA;tag=basic;transcript_id=ENST00000584487;version=1

#10	havana	lincRNA	73813052	73813795	.	-	.	ID=transcript:ENST00000609403;Parent=gene:ENSG00000272988;Name=RP11-150D20.5-001;biotype=lincRNA;havana_transcript=OTTHUMT00000473110;havana_version=1;tag=basic;transcript_id=ENST00000609403;version=1

print "\treading $refFile\n";
my $fhIn;
if ($refFile =~ /.gz$/) {
	use IO::Zlib;
	$fhIn = new IO::Zlib;
	$fhIn->open($refFile, "rb")
	}
else {
	open($fhIn, "<", $refFile) or die "could not read $refFile ($!)\n";
	}
my $allRefs = GFFfile2hash($fhIn,$chromName);
close($fhIn);

## selecting IDs from hash:
my %targets;
foreach my $id (@{$IDs}) {
	$id = uc($id);
	if (exists ${$allRefs}{"transcript"}{$id}) { %{ $targets{$id} } = %{ ${$allRefs}{"transcript"}{$id} }; }
	elsif (exists ${$allRefs}{"gene"}{$id}) {
		if (exists ${$allRefs}{"gene"}{$id}{"2gene_id"}) {
			my $gene_id = ${$allRefs}{"gene"}{$id}{"2gene_id"};
			foreach my $nm (keys%{ ${$allRefs}{"gene"}{$gene_id}{"transcript"} }) {
				if ($wNonCod || exists ${$allRefs}{"transcript"}{$nm}{"CDS_starts"}) {
					%{ $targets{$nm} } = %{ ${$allRefs}{"transcript"}{$nm} };
					}
				}
			}
		else {
			foreach my $nm (keys%{ ${$allRefs}{"gene"}{$id}{"transcript"} }) {
				if ($wNonCod || exists ${$allRefs}{"transcript"}{$nm}{"CDS_starts"}) {
					%{ $targets{$nm} } = %{ ${$allRefs}{"transcript"}{$nm} };
					}
				}
			}
		}
	elsif (exists ${$allRefs}{"Id_noExt"}{$id}) {
		foreach my $nm (keys%{ ${$allRefs}{"Id_noExt"}{$id} }) {
			if ($wNonCod || exists ${$allRefs}{"transcript"}{$nm}{"CDS_starts"}) {
				%{ $targets{$nm} } = %{ ${$allRefs}{"transcript"}{$nm} };
				}
			}
		}
	else { die "can't find $id in $refFile file\n"; }
	}

foreach my $id (keys%targets) {
	delete $targets{$id}{"ex_starts"};
	@{ $targets{$id}{"ex_starts"} } = sort{$a<=>$b}keys%{ ${$allRefs}{"transcript"}{$id}{"ex_starts"} };
	delete $targets{$id}{"ex_ends"};
	@{ $targets{$id}{"ex_ends"} } = sort{$a<=>$b}keys%{ ${$allRefs}{"transcript"}{$id}{"ex_ends"} };
	if (exists $targets{$id}{"CDS_starts"}) {
		delete $targets{$id}{"CDS_starts"};
		my @pos = sort{$a<=>$b}keys%{ ${$allRefs}{"transcript"}{$id}{"CDS_starts"} };
		$targets{$id}{"CDS_start"} = $pos[0];
		delete $targets{$id}{"CDS_ends"};
		@pos = sort{$a<=>$b}keys%{ ${$allRefs}{"transcript"}{$id}{"CDS_ends"} };
		$targets{$id}{"CDS_end"} = $pos[-1];
		}
	if (exists $targets{$id}{"gene_id"}) {
		my $gene_id = $targets{$id}{"gene_id"};
		if (exists ${$allRefs}{"gene"}{$gene_id}{"2gene_name"}) { $targets{$id}{"gene"} = ${$allRefs}{"gene"}{$gene_id}{"2gene_name"}; }
		else { $targets{$id}{"gene"} = $gene_id; }
		}
	}

return(\%targets);

}

####

sub GFFfile2hash {

my($fhIn,$chromName) = @_;

my(%allRefs,%idx);
$idx{"chr"} = 0;
$idx{"type"} = 2;
$idx{"start"} = 3;
$idx{"end"} = 4;
$idx{"strand"} = 6;
$idx{"info"} = 8;
while (my $line = <$fhIn>) {
	if ($line !~ /^#/) {
		chomp $line;
		my @tab = split(/\t/,$line);
		my $chr = $tab[$idx{"chr"}];
		$chr =~ s/^chr//i;
		${$chromName}{"refId"}{$chr} = $tab[$idx{"chr"}];
		my($gene_id,$transcript_id,$gene_name);
		if ($tab[$idx{"type"}] eq "gene") {
#GFF attributes
#ID=gene:ENSG00000188157;Name=AGRN;biotype=protein_coding;description=agrin [Source:HGNC Symbol%3BAcc:329];gene_id=ENSG00000188157;logic_name=ensembl_havana_gene;version=9
			if ($tab[$idx{"info"}] =~ /ID=gene:(\w+)/) {
				$gene_id = $1; $gene_id = uc($gene_id);
				if ($tab[$idx{"info"}] =~ /(|;)Name=(\w+)/) {
					$gene_name = $2; $gene_name = uc($gene_name);
					$allRefs{"gene"}{$gene_id}{"2gene_name"} = $gene_name;
					$allRefs{"gene"}{$gene_name}{"2gene_id"} = $gene_id;
					}
				}
			}
#		elsif ($tab[$idx{"type"}] =~ /^mRNA$/) {
#1       ensembl_havana  mRNA    955503  991496  .       +       .       ID=transcript:ENST00000379370;Parent=gene:ENSG00000188157;Name=AGRN-001;biotype=protein_coding;ccdsid=CCDS30551.1;havana_transcript=OTTHUMT00000097990;havana_version=2;tag=basic;transcript_id=ENST00000379370;version=2
#			if ($tab[$idx{"info"}] =~ /ID=transcript:(\w+);/) {
#				$transcript_id = $1; $transcript_id = uc($transcript_id);
#				if (exists $allRefs{"transcript"}{$transcript_id}) {
#					print "$transcript_id...$chr...$gene_name\n";
#					}
#				}
#			}
		elsif ($tab[$idx{"type"}] =~ /^(exon|CDS)$/) {
#GFF attributes
#exon :	Parent=transcript:ENST00000379370;Name=ENSE00002277560;constitutive=0;ensembl_end_phase=0;ensembl_phase=-1;exon_id=ENSE00002277560;rank=1;version=1
#CDS :	ID=CDS:ENSP00000368678;Parent=transcript:ENST00000379370;protein_id=ENSP00000368678
			if ($tab[$idx{"info"}] =~ /Parent=transcript:(\w+)/) {
				$transcript_id = $1; $transcript_id = uc($transcript_id);
				if ($tab[$idx{"type"}] eq "exon") {
					$allRefs{"transcript"}{$transcript_id}{"ex_starts"}{$tab[$idx{"start"}]} = 1;
					$allRefs{"transcript"}{$transcript_id}{"ex_ends"}{$tab[$idx{"end"}]} = 1;
					}
				else {
					$allRefs{"transcript"}{$transcript_id}{"CDS_starts"}{$tab[$idx{"start"}]} = 1;
					$allRefs{"transcript"}{$transcript_id}{"CDS_ends"}{$tab[$idx{"end"}]} = 1;
					}
				}
			}
		else {
			if ($tab[$idx{"info"}] =~ /ID=transcript:(\w+)/) {
#ID=transcript:ENST00000379370;Parent=gene:ENSG00000188157;Name=AGRN-001;biotype=protein_coding;ccdsid=CCDS30551.1;havana_transcript=OTTHUMT00000097990;havana_version=2;tag=basic;transcript_id=ENST00000379370;version=2
				$transcript_id = $1; $transcript_id = uc($transcript_id);
				my $Id_noExt = $transcript_id;
				$Id_noExt =~ s/\.(\d+)$//;;
				$allRefs{"Id_noExt"}{$Id_noExt}{$transcript_id} = 1;
				if ($tab[$idx{"info"}] =~ /Parent=gene:(\w+)/) {
					$gene_id = $1; $gene_id = uc($gene_id);
					$allRefs{"transcript"}{$transcript_id}{"gene_id"} = $gene_id;
					$allRefs{"gene"}{$gene_id}{"transcript"}{$transcript_id} = 1;
					}
				if ($tab[$idx{"info"}] =~ /biotype=(\w+)/) {
					$allRefs{"transcript"}{$transcript_id}{"biotype"} = $1;
					}
				$allRefs{"transcript"}{$transcript_id}{"chr"} = $chr;
				$allRefs{"transcript"}{$transcript_id}{"strand"} = $tab[$idx{"strand"}];
				}
			}

		}
	}
return(\%allRefs);

}


##############################

sub bed2IDs {

my($allRefs,$refFile,$refFmt,$interval_r,$chromName,$wNonCod) = @_;
#my %interval ->	#$Bed{$chr}{$start} = $end , in 1-based

print "\n\tlooking for genes/transcripts overlapping bed file intervals\n";

#@{ $RefCoord{$chr}{$start}{$end} } = [$NM1,$NM2,...]
my($fhIn,%RefCoord);
if (!$allRefs) {
	print "\t\treading $refFile\n";
	if ($refFile =~ /.gz$/) {
		use IO::Zlib;
		$fhIn = new IO::Zlib;
		$fhIn->open($refFile, "rb")
		}
	else {
		open($fhIn, "<", $refFile) or die "could not read $refFile ($!)\n";
		}
	}
if ($refFmt eq "ucsc") {
	unless ($allRefs) { $allRefs = UCSCfile2hash($fhIn,$chromName); }
	foreach my $nm (keys%{ ${$allRefs}{"transcript"} }) {
		if ($wNonCod || exists ${$allRefs}{"transcript"}{$nm}{"CDS_start"}) {
			push(@{ $RefCoord{${$allRefs}{"transcript"}{$nm}{"chr"}}{${$allRefs}{"transcript"}{$nm}{"ex_starts"}[0]}{${$allRefs}{"transcript"}{$nm}{"ex_ends"}[-1]} }, $nm);
			}
		}
	}
elsif ($refFmt eq "gtf") {
	unless ($allRefs) { $allRefs = GTFfile2hash($fhIn,$chromName); }
	foreach my $nm (keys%{ ${$allRefs}{"transcript"} }) {
		if ($wNonCod || exists ${$allRefs}{"transcript"}{$nm}{"CDS_starts"}) {
			my @ex_starts = sort{$a<=>$b}keys%{ ${$allRefs}{"transcript"}{$nm}{"ex_starts"} };
			delete ${$allRefs}{"transcript"}{$nm}{"ex_starts"};
			@{ ${$allRefs}{"transcript"}{$nm}{"ex_starts"} } = @ex_starts;
			my @ex_ends = sort{$a<=>$b}keys%{ ${$allRefs}{"transcript"}{$nm}{"ex_ends"} };
			delete ${$allRefs}{"transcript"}{$nm}{"ex_ends"};
			@{ ${$allRefs}{"transcript"}{$nm}{"ex_ends"} } = @ex_ends;
			push(@{ $RefCoord{${$allRefs}{"transcript"}{$nm}{"chr"}}{$ex_starts[0]}{$ex_ends[-1]} }, $nm);
			}
		}
	}
elsif ($refFmt eq "gff3") {
	unless ($allRefs) { $allRefs = GFFfile2hash($fhIn,$chromName); }
	foreach my $nm (keys%{ ${$allRefs}{"transcript"} }) {
		if ($wNonCod || exists ${$allRefs}{"transcript"}{$nm}{"CDS_starts"}) {
			my @ex_starts = sort{$a<=>$b}keys%{ ${$allRefs}{"transcript"}{$nm}{"ex_starts"} };
			delete ${$allRefs}{"transcript"}{$nm}{"ex_starts"};
			@{ ${$allRefs}{"transcript"}{$nm}{"ex_starts"} } = @ex_starts;
			my @ex_ends = sort{$a<=>$b}keys%{ ${$allRefs}{"transcript"}{$nm}{"ex_ends"} };
			delete ${$allRefs}{"transcript"}{$nm}{"ex_ends"};
			@{ ${$allRefs}{"transcript"}{$nm}{"ex_ends"} } = @ex_ends;
			push(@{ $RefCoord{${$allRefs}{"transcript"}{$nm}{"chr"}}{$ex_starts[0]}{$ex_ends[-1]} }, $nm);
			}
		}
	}
if ($fhIn) { close($fhIn); }

my %targets;
foreach my $chr (sort(keys%{$interval_r})) {
	my @bedStarts = sort{$a<=>$b}keys%{ ${$interval_r}{$chr} };
	my @refStarts = sort{$a<=>$b}keys%{ $RefCoord{$chr} };
	my $c = 0;		#idx of @refStarts
	my @RefEnds = sort{$b<=>$a}keys%{ $RefCoord{$chr}{$refStarts[$c]} };
	foreach my $bedstart(@bedStarts) {
		while ( ($c < $#refStarts) && ($bedstart > $RefEnds[0]) ) { 
			$c++; 
			@RefEnds = sort{$b<=>$a}keys%{ $RefCoord{$chr}{$refStarts[$c]} };
			}
		my $c2 = $c;
		while (($c2 < $#refStarts) && (${$interval_r}{$chr}{$bedstart} >= $refStarts[$c2])) {
			foreach my $end (@RefEnds) {
				if ($bedstart <= $end) { 
					foreach (@{ $RefCoord{$chr}{$refStarts[$c2]}{$end} }) {
						%{ $targets{$_} } = %{ ${$allRefs}{"transcript"}{$_} };
						}
					}
				else { last; }
				}
			$c2++;
			@RefEnds = sort{$b<=>$a}keys%{ $RefCoord{$chr}{$refStarts[$c2]} };
			}
		if (${$interval_r}{$chr}{$bedstart} >= $refStarts[$c2]) {
			foreach my $end (@RefEnds) {
				if ($bedstart <= $end) { 
					foreach (@{ $RefCoord{$chr}{$refStarts[$c2]}{$end} }) {
						%{ $targets{$_} } = %{ ${$allRefs}{"transcript"}{$_} };
						}
					}
				else { last; }
				}
			}
		}
	}
if ($refFmt eq "gtf" || $refFmt eq "gff3") {
	foreach my $id (keys%targets) {
		if (exists $targets{$id}{"CDS_starts"}) {
			my @pos = sort{$a<=>$b}keys%{ $targets{$id}{"CDS_starts"} };
			delete $targets{$id}{"CDS_starts"};
			$targets{$id}{"CDS_start"}= $pos[0];
			@pos = sort{$a<=>$b}keys%{ $targets{$id}{"CDS_ends"} };
			delete $targets{$id}{"CDS_ends"};
			$targets{$id}{"CDS_end"} = $pos[-1];
			}
		if ($refFmt eq "gff3") {
			if (exists $targets{$id}{"gene_id"}) {
				my $gene_id = $targets{$id}{"gene_id"};
				if (exists ${$allRefs}{"gene"}{$gene_id}{"2gene_name"}) { $targets{$id}{"gene"} = ${$allRefs}{"gene"}{$gene_id}{"2gene_name"}; }
				else { $targets{$id}{"gene"} = $gene_id; }
				}
			}
		}
	}

return(\%targets);

}


##############################

sub ref2Hash {

my($refFile,$refFmt,$chromName) = @_;

print "\t\treading $refFile\n";

my($fhIn,$allRefs);
if ($refFile =~ /.gz$/) {
	use IO::Zlib;
	$fhIn = new IO::Zlib;
	$fhIn->open($refFile, "rb")
	}
else {
	open($fhIn, "<", $refFile) or die "could not read $refFile ($!)\n";
	}

if ($refFmt eq "ucsc") {
	$allRefs = UCSCfile2hash($fhIn,$chromName);
	}
elsif ($refFmt eq "gtf") {
	$allRefs = GTFfile2hash($fhIn,$chromName);
	}
elsif ($refFmt eq "gff3") {
	$allRefs = GFFfile2hash($fhIn,$chromName);
	}
close($fhIn);

return($allRefs);

}

##############################

sub refHash2Transcripts {

##put back in 0-based : easier to compare with bed file!

my($allRefs,$refFmt) = @_;

my(%targets);
if ($refFmt eq "ucsc") {
	foreach my $nm (keys%{ ${$allRefs}{"transcript"} }) {
		%{ $targets{$nm} } = %{ ${$allRefs}{"transcript"}{$nm} };
		delete $targets{$nm}{"ex_starts"};
		foreach (@{ ${$allRefs}{"transcript"}{$nm}{"ex_starts"} }) { push(@{ $targets{$nm}{"ex_starts"} }, ($_-1)); }
		if (exists $targets{$nm}{"CDS_start"}) { $targets{$nm}{"CDS_start"} -= 1; }
#		else {
#			$targets{$nm}{"CDS_start"} = $targets{$nm}{"ex_ends"}[-1];
#			$targets{$nm}{"CDS_end"} = $targets{$nm}{"ex_ends"}[-1];
#			}
		}
	}
elsif ($refFmt eq "gtf") {
	foreach my $nm (keys%{ ${$allRefs}{"transcript"} }) {
		%{ $targets{$nm} } = %{ ${$allRefs}{"transcript"}{$nm} };
		my @ex_starts = sort{$a<=>$b}keys%{ ${$allRefs}{"transcript"}{$nm}{"ex_starts"} };
		delete $targets{$nm}{"ex_starts"};
		foreach (@ex_starts) { push(@{ $targets{$nm}{"ex_starts"} }, ($_-1)); }
		my @ex_ends = sort{$a<=>$b}keys%{ ${$allRefs}{"transcript"}{$nm}{"ex_ends"} };
		delete $targets{$nm}{"ex_ends"};
		@{ $targets{$nm}{"ex_ends"} } = @ex_ends;
		if (exists $targets{$nm}{"CDS_starts"}) {
			my @pos = sort{$a<=>$b}keys%{ ${$allRefs}{"transcript"}{$nm}{"CDS_starts"} };
			delete $targets{$nm}{"CDS_starts"};
			$targets{$nm}{"CDS_start"} = $pos[0]-1;
			@pos = sort{$a<=>$b}keys%{ ${$allRefs}{"transcript"}{$nm}{"CDS_ends"} };
			delete $targets{$nm}{"CDS_ends"};
			$targets{$nm}{"CDS_end"} = $pos[-1];
			}
#		else {
#			$targets{$nm}{"CDS_start"} = $ex_ends[-1];
#			$targets{$nm}{"CDS_end"} = $ex_ends[-1];
#			}
		}
	}
elsif ($refFmt eq "gff3") {
	foreach my $nm (keys%{ ${$allRefs}{"transcript"} }) {
		%{ $targets{$nm} } = %{ ${$allRefs}{"transcript"}{$nm} };
		my @ex_starts = sort{$a<=>$b}keys%{ ${$allRefs}{"transcript"}{$nm}{"ex_starts"} };
		delete $targets{$nm}{"ex_starts"};
		foreach (@ex_starts) { push(@{ $targets{$nm}{"ex_starts"} }, ($_-1)); }
		my @ex_ends = sort{$a<=>$b}keys%{ ${$allRefs}{"transcript"}{$nm}{"ex_ends"} };
		delete $targets{$nm}{"ex_ends"};
		@{ $targets{$nm}{"ex_ends"} } = @ex_ends;
		if (exists $targets{$nm}{"CDS_starts"}) {
			delete $targets{$nm}{"CDS_starts"};
			my @pos = sort{$a<=>$b}keys%{ ${$allRefs}{"transcript"}{$nm}{"CDS_starts"} };
			$targets{$nm}{"CDS_start"} = $pos[0]-1;
			delete $targets{$nm}{"CDS_ends"};
			@pos = sort{$a<=>$b}keys%{ ${$allRefs}{"transcript"}{$nm}{"CDS_ends"} };
			$targets{$nm}{"CDS_end"} = $pos[-1];
			}
#		else {
#			$targets{$nm}{"CDS_start"} = $ex_ends[-1];
#			$targets{$nm}{"CDS_end"} = $ex_ends[-1];
#			}
		if (exists ${$allRefs}{"transcript"}{$nm}{"gene_id"}) {
			my $gene_id = ${$allRefs}{"transcript"}{$nm}{"gene_id"};
			if (exists ${$allRefs}{"gene"}{$gene_id}{"2gene_name"}) { $targets{$nm}{"gene"} = ${$allRefs}{"gene"}{$gene_id}{"2gene_name"}; }
			else { $targets{$nm}{"gene"} = $gene_id; }
			}
		}
	}

return(\%targets);

}

##############################

#RefSeq:
##bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames
#0	NM_032291	chr1	+	66999824	67210768	67000041	67208778	25	66999824,67091529,67098752,67101626,67105459,67108492,67109226,67126195,67133212,67136677,67137626,67138963,67142686,67145360,67147551,67154830,67155872,67161116,67184976,67194946,67199430,67205017,67206340,67206954,67208755,	67000051,67091593,67098777,67101698,67105516,67108547,67109402,67126207,67133224,67136702,67137678,67139049,67142779,67145435,67148052,67154958,67155999,67161176,67185088,67195102,67199563,67205220,67206405,67207119,67210768,	0	SGIP1	cmpl	cmpl	0,1,2,0,0,0,1,0,0,0,1,2,1,1,1,1,0,1,1,2,2,0,2,1,1,
## @hashSub = Id2Coord($idFile,$len5,$len3,$upstream,$downstream,$splitBedFromId,$mergeBedFromId,$wNonCod,$wUTR,$id2Bed,\@chromOrder);

sub Id2Coord {

my($IDinRef,$len5,$len3,$upstream,$downstream,$noOverlap,$mergeBedFromId,$wUTR,$id2Bed,$chromOrder_r,$chromName_r) = @_;

print "\n\tgetting positions from gene IDs\n";

my(%NMgene,%geneNM,%NMchr,%NMsens,%NMstartCod,%NMendCod,%NM_Ex,%Regions,%gOK,%hashBed,%printBed);

foreach my $NM (keys%{$IDinRef}) {
	my $chr = ${$IDinRef}{$NM}{"chr"};
	my $gene = ${$IDinRef}{$NM}{"gene"};

	## if same gene name on different chr? add -chr. at the end of this gene
#	if ( (exists $geneNM{$gene}) && ($chr ne $NMchr{$geneNM{$gene}[0]}) ) { $tab[$geneIdx] .= "-$chr"; }
	##if same transcript on different chr ? skip this one
#	if (exists $NMchr{$NM}) {
#		if ($chr =~ /^$NMchr{$NM}.+/) { $ok = 0; }
#		unless ($NMchr{$NM} =~ /$chr.+/) { $ok = 0; }
#		}

	$NMgene{$NM} = $gene;
	push(@{ $geneNM{$gene} }, $NM);
	$NMchr{$NM} = $chr;
	$NMsens{$NM} = ${$IDinRef}{$NM}{"strand"};
	my @Starts = @{ ${$IDinRef}{$NM}{"ex_starts"} };
	my @Ends = @{ ${$IDinRef}{$NM}{"ex_ends"} };
	if (exists ${$IDinRef}{$NM}{"CDS_start"}) {
		$NMstartCod{$NM} = ${$IDinRef}{$NM}{"CDS_start"};
		$NMendCod{$NM} = ${$IDinRef}{$NM}{"CDS_end"};
		}
	else {
		if (!$wUTR) { next; }		#next if non coding transcript and $wUTR not wanted
		$NMstartCod{$NM} = $Ends[-1];
		$NMendCod{$NM} = $Ends[-1];
		}

	my(%interval,$firstCodingEx,$lastCodingEx,$minLen,$plusLen,$exonNb); # %interval in 1-based coord
	if ($NMsens{$NM} eq "+") { $minLen = $len5; $plusLen = $len3; }
	else { $minLen = $len3; $plusLen = $len5; }

	#expand exons +/-length:
	if ($wUTR) {

		##intervals
		if ($upstream) { $interval{($Starts[0]-$upstream)} = $Ends[0]+$plusLen; }
		else { $interval{($Starts[0]-$minLen)} = $Ends[0]+$plusLen; }
		for (my $i=1;$i<$#Starts;$i++) {
			$interval{($Starts[$i]-$minLen)} = $Ends[$i]+$plusLen;
			}
		if ($downstream) { $interval{($Starts[-1]-$minLen)} = $Ends[-1]+$downstream; }
		else { $interval{($Starts[-1]-$minLen)} = $Ends[-1]+$plusLen; }

		##bed
		if ($id2Bed) {
			if ($NMsens{$NM} eq "+" && $upstream) {
				$printBed{$chr}{($Starts[0]-1-$upstream)}{($Starts[0]-1-$minLen)}{$gene}{$NM} = "upstream";
				}
			elsif ($NMsens{$NM} eq "-" && $downstream) {
				$printBed{$chr}{($Starts[0]-1-$downstream)}{($Starts[0]-1-$minLen)}{$gene}{$NM} = "downstream";
				}
			for (my $i=0;$i<scalar(@Starts);$i++) {
				if ($noOverlap && $Starts[$i+1] && ($Ends[$i]+$plusLen) > ($Starts[$i+1]-1-$minLen)) {
					my $ecart  = ($Starts[$i+1]-1-$Ends[$i]) / 2;
					$Ends[$i] = int($Ends[$i]+$ecart-$plusLen);
					$Starts[$i+1] = int($Starts[$i+1]-1-$ecart+$minLen);
					}
				if ($NMsens{$NM} eq "+") { $exonNb = ($i+1); }
				else { $exonNb = (scalar(@Starts)-$i); }
				$printBed{$chr}{($Starts[$i]-1-$minLen)}{($Ends[$i]+$plusLen)}{$gene}{$NM} = "exon$exonNb";
				if (($Ends[$i] < $NMstartCod{$NM}) || (($Starts[$i]-1) > $NMendCod{$NM})) {
					$printBed{$chr}{($Starts[$i]-1-$minLen)}{($Ends[$i]+$plusLen)}{$gene}{$NM} .= ".NC";
					}
				}
			if ($NMsens{$NM} eq "+" && $downstream) {
				$printBed{$chr}{($Ends[-1]+$plusLen)}{($Ends[-1]+$downstream)}{$gene}{$NM} = "downstream";
				}
			elsif ($NMsens{$NM} eq "-" && $upstream) {
				$printBed{$chr}{($Ends[-1]+$plusLen)}{($Ends[-1]+$upstream)}{$gene}{$NM} = "upstream";
				}
			}
		}

	else {

		##intervals
		for (my $i=0;$i<scalar(@Starts);$i++) {
			if ($Ends[$i] < $NMstartCod{$NM}) { next; }
			elsif ( ($Ends[$i] >= $NMstartCod{$NM}) && ($Starts[$i] < $NMstartCod{$NM}) ) {
				$firstCodingEx = $i;
				if ($Ends[$i] > $NMendCod{$NM}) {
					$interval{($NMstartCod{$NM}-$minLen)} = $NMendCod{$NM}+$plusLen;
					}
				else {
					$interval{($NMstartCod{$NM}-$minLen)} = $Ends[$i]+$plusLen;
					}
				}
			elsif ( ($Starts[$i] >= $NMstartCod{$NM}) && ($Starts[$i] <= $NMendCod{$NM}) ) {
				if ($Ends[$i] < $NMendCod{$NM}) {
					$interval{($Starts[$i]-$minLen)} = $Ends[$i]+$plusLen;
					}
				else {
					$lastCodingEx = $i;
					$interval{($Starts[$i]-$minLen)} = $NMendCod{$NM}+$plusLen;
					}
				}
			else { last; }
			}

		##bed
		if ($id2Bed) {
			for (my $i=0;$i<scalar(@Starts);$i++) {
				if ($noOverlap && $Starts[$i+1] && ($Ends[$i]+$plusLen) > ($Starts[$i+1]-1-$minLen)) {
					my $ecart  = ($Starts[$i+1]-1-$Ends[$i]) / 2;
					$Ends[$i] = int($Ends[$i]+$ecart-$plusLen);
					$Starts[$i+1] = int($Starts[$i+1]-1-$ecart+$minLen);
					}
				if ($NMsens{$NM} eq "+") { $exonNb = ($i+1); }
				else { $exonNb = (scalar(@Starts)-$i); }
				if ($Ends[$i] < $NMstartCod{$NM}) { next; }
				elsif ( ($Ends[$i] >= $NMstartCod{$NM}) && ($Starts[$i] < $NMstartCod{$NM}) ) {
					$firstCodingEx = $i;
					if ($Ends[$i] > $NMendCod{$NM}) {
						$printBed{$chr}{($NMstartCod{$NM}-1-$minLen)}{($NMendCod{$NM}+$plusLen)}{$gene}{$NM} = "exon$exonNb";
						}
					else {
						$printBed{$chr}{($NMstartCod{$NM}-1-$minLen)}{($Ends[$i]+$plusLen)}{$gene}{$NM} = "exon$exonNb";
						}
					}
				elsif ( ($Starts[$i] >= $NMstartCod{$NM}) && ($Starts[$i] <= $NMendCod{$NM}) ) {
					if ($Ends[$i] < $NMendCod{$NM}) {
						$printBed{$chr}{($Starts[$i]-1-$minLen)}{($Ends[$i]+$plusLen)}{$gene}{$NM} = "exon$exonNb";
						}
					else {
						$lastCodingEx = $i;
						$printBed{$chr}{($Starts[$i]-1-$minLen)}{($NMendCod{$NM}+$plusLen)}{$gene}{$NM} = "exon$exonNb";
						}
					}
				else { last; }
				}
			}
		}

	#eventually merges overlapping intervals of regions
	#associates real exons to intervals
	my @Starts2 = sort{$a<=>$b}(keys%interval);
	my %interval2;
	my $start = $Starts2[0];
	my $end = $interval{$start};
	if ($wUTR) {
		$NM_Ex{$NM}{$start}{$Starts[0]} = $Ends[0];
		for (my $i=1;$i<scalar(@Starts2);$i++) {
			if ($Starts2[$i] <= $end) { 
				$end = $interval{$Starts2[$i]}; 
				$NM_Ex{$NM}{$start}{$Starts[$i]} = $Ends[$i];
				}
			else { 
				$interval2{$start} = $end;
				$start = $Starts2[$i];
				$end = $interval{$start};
				$NM_Ex{$NM}{$start}{$Starts[$i]} = $Ends[$i];
				}
			}
		$interval2{$start} = $end;
		$NM_Ex{$NM}{$start}{$Starts[-1]} = $Ends[-1];
		}
	else {
		$NM_Ex{$NM}{$start}{$NMstartCod{$NM}} = $Ends[$firstCodingEx];
		for (my $i=1;$i<scalar(@Starts2);$i++) {
			if ($Starts2[$i] <= $end) { 
				$end = $interval{$Starts2[$i]}; 
				$NM_Ex{$NM}{$start}{$Starts[$i+$firstCodingEx]} = $Ends[$i+$firstCodingEx];
				}
			else { 
				$interval2{$start} = $end;
				$start = $Starts2[$i];
				$end = $interval{$start};
				$NM_Ex{$NM}{$start}{$Starts[$i+$firstCodingEx]} = $Ends[$i+$firstCodingEx];
				}
			}
		$interval2{$start} = $end;
		if ($lastCodingEx) {
			$NM_Ex{$NM}{$start}{$Starts[$lastCodingEx]} = $NMendCod{$NM};
			}
		else { print "$NM\t $firstCodingEx\n"; exit;}
		}
	%{ $Regions{$chr}{$NM} } = %interval2;
	foreach (keys%interval2) {
		if ( exists $hashBed{$chr}{$_} ) {
			if ( $interval2{$_} > $hashBed{$chr}{$_}) {
				$hashBed{$chr}{$_} = $interval2{$_};
				}
			}
		else {
			$hashBed{$chr}{$_} = $interval2{$_};
			}
		}
	}


if ($id2Bed) {
	my @allChrom = ();
	if (@{$chromOrder_r}) {
		foreach my $chr (@{$chromOrder_r}) {
			if (exists $printBed{$chr}) { push(@allChrom,$chr); }
			}
		}
	else {
		foreach my $chr (sort(keys%printBed)) { push(@allChrom,$chr); }
		}

	if (exists ${$chromName_r}{"fasta"}) { %{ ${$chromName_r}{"bed"} } = %{ ${$chromName_r}{"fasta"} }; }
	else { %{ ${$chromName_r}{"bed"} } = %{ ${$chromName_r}{"refId"} }; }

	open(my $fh, ">", "$id2Bed") || die "can't create $id2Bed ($!)\n";
	if ($mergeBedFromId) {
		foreach my $chr (@allChrom) {
			#exons with same start
			my %printBed2 = ();
			foreach my $start (keys%{ $printBed{$chr} }) {
				my @ends = sort{$a<=>$b}keys%{ $printBed{$chr}{$start} };
				$printBed2{$start}{"end"} = $ends[-1];
				foreach my $end (keys%{ $printBed{$chr}{$start} }) {
					foreach my $gene (keys%{ $printBed{$chr}{$start}{$end} }) {
						foreach my $NM (keys%{ $printBed{$chr}{$start}{$end}{$gene} }) {
							push(@{ $printBed2{$start}{"info"}{$gene}{$NM} }, $printBed{$chr}{$start}{$end}{$gene}{$NM});
							}
						}
					}	
				}
			#overlapping exons
			my @Starts = sort{$a<=>$b}(keys%printBed2);
			my $start = $Starts[0];
			my $end = $printBed2{$start}{"end"};
			my %printBed3 = ();
			foreach my $gene (keys%{ $printBed2{$start}{"info"} }) {
				foreach my $NM (keys%{ $printBed2{$start}{"info"}{$gene} }) {
					push(@{ $printBed3{$gene}{$NM} }, @{ $printBed2{$start}{"info"}{$gene}{$NM} });
					}
				}
			for (my $i=1;$i<scalar(@Starts);$i++) {
				if ($Starts[$i] <= $end) {
					foreach my $gene (keys%{ $printBed2{$Starts[$i]}{"info"} }) {
						foreach my $NM (keys%{ $printBed2{$Starts[$i]}{"info"}{$gene} }) {
							push(@{ $printBed3{$gene}{$NM} }, @{ $printBed2{$Starts[$i]}{"info"}{$gene}{$NM} });
							}
						}
					if ($printBed2{$Starts[$i]}{"end"} > $end)
						{ $end = $printBed2{$Starts[$i]}{"end"}; }
					else { next; }
					}
				else {
					printID2BedMerge($fh,${$chromName_r}{"bed"}{$chr},$start,$end,\%printBed3);
					%printBed3 = ();
					foreach my $gene (keys%{ $printBed2{$Starts[$i]}{"info"} }) {
						foreach my $NM (keys%{ $printBed2{$Starts[$i]}{"info"}{$gene} }) {
							push(@{ $printBed3{$gene}{$NM} }, @{ $printBed2{$Starts[$i]}{"info"}{$gene}{$NM} });
							}
						}
					$start = $Starts[$i];
					$end = $printBed2{$start}{"end"};
					}
				}
			printID2BedMerge($fh,${$chromName_r}{"bed"}{$chr},$start,$end,\%printBed3);
			}
		}

	else {
		foreach my $chr (@allChrom) {
			foreach my $start (sort{$a<=>$b}keys%{ $printBed{$chr} }) {
				foreach my $end (sort{$a<=>$b}keys%{ $printBed{$chr}{$start} }) {
					my $info="";
					foreach my $gene (sort(keys%{ $printBed{$chr}{$start}{$end} })) {
						$info .= "$gene:";
						foreach my $NM (sort(keys%{ $printBed{$chr}{$start}{$end}{$gene} })) {
							$info .= "$NM:".$printBed{$chr}{$start}{$end}{$gene}{$NM}.",";
							}
						chop $info;
						$info .= ";";
						}
					chop $info;
					print $fh ${$chromName_r}{"bed"}{$chr}."\t$start\t$end\t$info\n";
					}	
				}
			}
		}

	close $fh;
	}

my %Genes;
foreach (keys%geneNM) { 
	$Genes{$_} = 1;
	@{ $geneNM{$_} } = sort@{ $geneNM{$_} };
	}
return (\%geneNM,\%NMgene,\%NMchr,\%NMsens,\%NMstartCod,\%NMendCod,\%Regions,\%NM_Ex,\%Genes,\%hashBed);
}

####
sub printID2BedMerge {
my($fh,$chr,$start,$end,$infoRef) = @_;
my $info = "";
foreach my $gene (sort(keys%{$infoRef})) {
	$info .= "$gene:";
	foreach my $NM (sort(keys%{ ${$infoRef}{$gene} })) {
		$info .= "$NM:";
		foreach my $ex (@{ ${$infoRef}{$gene}{$NM} }) {
			$info .= "$ex-";
			}
		chop $info;
		$info .= ",";
		}
	chop $info;
	$info .= ";";
	}
chop $info;
print $fh "$chr\t$start\t$end\t$info\n";
}


##############################


1;



