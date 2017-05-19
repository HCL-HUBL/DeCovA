package Bio::NGS::HCL::DeCovA::CNV_tool;


use strict;
use warnings;


#@hashSub = Bio::NGS::HCL::DeCovA::CNV_tool::CNV_detect(\@Files,\%sName2,"$outdir/CNV_analysis",$nGraf,$fichier_sexe,\%CNV_opt,$mbq,$mmq,$refBedLines,\@ChromOrder,\%fai);
sub CNV_detect
{
my($h1,$h2,$outdir,$nGraf,$fichier_sexe,$h3,$mbq,$mmq,$h4,$h5,$h6)=@_;
my@Files = @$h1;
my%sampleName = %$h2;
my%CNV_opt = %$h3;
my@Regions = @$h4;
my@ChromOrder = @$h5;
my%fai = %$h6;
my$norm = $CNV_opt{"norm"};
my$normByGender = $CNV_opt{"normByGender"};
my$Ref = $CNV_opt{"RefDepth"};
my$RefByGender = $CNV_opt{"RefByGender"};
my$seuil_region = $CNV_opt{"seuil_region"};
if ($seuil_region < (1/scalar@Files)) { $seuil_region = (1/scalar@Files); }
my$seuil_patient = $CNV_opt{"seuil_patient"};
if ($seuil_patient < (1/scalar@Regions)) { $seuil_patient = (1/scalar@Regions); }
my$minCov = $CNV_opt{"seuil_cov"};
my$seuil_deletion = $CNV_opt{"seuil_deletion"};
my$seuil_duplication = $CNV_opt{"seuil_duplication"};
my$minCNV = $CNV_opt{"min_following_CNV"};
my$minDP = $CNV_opt{"min_DP"} ;
my$maxNonCNV = $CNV_opt{"max_Non_CNV"};
my$maxNonCNVrate = $CNV_opt{"max_Non_CNV_rate"};
my$CNVgraph = $CNV_opt{"chromGraph"};
my@cnvFields = @{ $CNV_opt{"fields"} };
my$ok = 0;
foreach my$i (0..$#cnvFields) {
	if ($cnvFields[$i] eq $norm) { 
		$ok = 1;
		if ($norm eq "std") { @cnvFields = (@cnvFields[0..($i-1)],"moy",@cnvFields[$i..$#cnvFields]); }
		last;
		}
	}
my@cnvVal = @cnvFields;		# = @cnvFields plus $norm (if not already present)
unless ($ok) { 
	if ($norm eq "std") { push(@cnvVal, ("moy","std")); }
	else { push(@cnvVal, $norm); }
	}

my$maxDepthGraph = 10;

my$nCol = scalar(split(/\t/,$Regions[0]{"allLine"}));

my%Patients;

#read gender file
if ($fichier_sexe) {
	open(PATIENTS, "$fichier_sexe") or die "Fichier $fichier_sexe impossible a lire\n";
	foreach my $ligne (<PATIENTS>) {
		$ligne =~ m/(\S+)\s+(\S+)/;
		my$ok="";
		foreach my$file (@Files) {
			if ($1 eq $sampleName{$file}) {
				$Patients{$file}{"Sexe"} = $2;
				$ok=1; last;
				}
			}
		unless ($ok) { print "sample $1 not found in bam files\n"; }
		}
	close(PATIENTS);
	foreach my$file (keys%Patients) {
		if (exists $Patients{$file}{"Sexe"}) {
			if ($Patients{$file}{"Sexe"} =~ /^m|^h/i) { $Patients{$file}{"Sexe"} = "H"; }
			else { $Patients{$file}{"Sexe"} = "F"; }
			}
		}
	}


#CALCUL de LA PROFONDEUR TOTALE POUR CHAQUE PATIENT, EN TENANT UNIQUEMENT COMPTE DES REGIONS CORRECTEMENT SEQUENCEES
foreach my$file (@Files) {
	$Patients{$file}{"tot_Autosomes_depth"} = 0;
	$Patients{$file}{"tot_chrX_depth"} = 0;
	$Patients{$file}{"tot_chrY_depth"} = 0;
	$Patients{$file}{"Ref_depth"} = 0;
	}

my$n_Regions = 0;
my$autosom_Regions = 0;
my$gonosom_Regions = 0;

#$bedLines[$i]{"allLine"} = $line;
#$bedLines[$i]{"Chr"} = $tab[0]; $bedLines[$i]{"Start"} = ($tab[1]+1); $bedLines[$i]{"End"} = $tab[2];
#$bedLines[$i]{$file}{"Cov"}
#$bedLines[$i]{$file}{"Mean"}
for my$r (0..$#Regions) {
	my $region_a_conserver = 1;
	if ($minCov) {
		foreach my$file (@Files) {
			if ( $Regions[$r]{$file}{"Cov"} >= $minCov)
				{ $region_a_conserver = 0; last; }
			}
		}
	if ($minDP) {
		my@allDepth = ();
		foreach my$file (@Files) {
			push(@allDepth, $Regions[$r]{$file}{"Mean"});
			}
		my$normDepth = norm($norm,\@allDepth);
		if ( $normDepth < $minDP)
			{ $region_a_conserver = 0; }
		}
	# si region à conserver, calcul Profondeur_Autosomes/gonosomes
	if ($region_a_conserver) {
		foreach my$file (@Files) {
			if ($Regions[$r]{"Chrom"} !~ /^(chr[XY]|[XY])$/) {
				$autosom_Regions++;
				if ($Ref eq "mean") { $Patients{$file}{"tot_Autosomes_depth"} += $Regions[$r]{$file}{"Mean"}; }
				}
			else {
				if ($Regions[$r]{"Chrom"} =~ /^(chrX|X)$/) {
					$gonosom_Regions++;
					if ($Ref eq "mean") { $Patients{$file}{"tot_chrX_depth"} += $Regions[$r]{$file}{"Mean"}; }
					}
				else {
					if ($Ref eq "mean") { $Patients{$file}{"tot_chrY_depth"} += $Regions[$r]{$file}{"Mean"}; }
					}
				}
			}
		}
	# Si region mal couverte, on l'imprime dans un fichier POUBELLE
	else
		{ $Regions[$r]{"Appel"} = "Couverture_Faible"; }
	}

my$depthTxt="";
my$sexTxt="";
foreach my$file (@Files) {

	##if tot bases
	if ($Ref eq "tot") {
		#autosomes
		my$cmd = "samtools view -F 0x4";
		if ($mmq) { $cmd .= " -q $mmq"; }
		$cmd .= " $file";
		foreach my$chr (@{ $fai{$file} }) {
			unless ($chr =~ /^(chr[XY]|[XY])$/)
				{ $cmd .= " $chr"; }
			}
		$cmd .= " | perl -ne 'END { print \"\$cnt\" } \@fields = split (/\\t/,\$_); \@bases = split (//,\$fields[9]);";
		if ($mbq) { $cmd .= " \@quals = map { unpack(\"C*\", \$_ )} split(//,\$fields[10]); for \$i (0..\$#bases) { if ((\$bases[\$i] =~ /[ACGT]/i) && (\$quals[\$i] >= $mbq)) { \$cnt++; } }'"; }
		else { $cmd .= " for \$i (0..\$#bases) { if (\$bases[\$i] =~ /[ACGT]/i) { \$cnt++; } }'"; }
		print "$cmd\n";
		$Patients{$file}{"tot_Autosomes_depth"} = `$cmd`;
		chomp $Patients{$file}{"tot_Autosomes_depth"};
		#chrX
		$cmd = "samtools view -F 0x4";
		if ($mmq) { $cmd .= " -q $mmq"; }
		$cmd .= " $file";
		foreach my$chr (@{ $fai{$file} }) {
			if ($chr =~ /^(chrX|X)$/)
				{ $cmd .= " $chr"; }
			}
		$cmd .= " | perl -ne 'END { print \"\$cnt\" } \@fields = split (/\\t/,\$_); \@bases = split (//,\$fields[9]);";
		if ($mbq) { $cmd .= " \@quals = map { unpack(\"C*\", \$_ )} split(//,\$fields[10]); for \$i (0..\$#bases) { if ((\$bases[\$i] =~ /[ACGT]/i) && (\$quals[\$i] >= $mbq)) { \$cnt++; } }'"; }
		else { $cmd .= " for \$i (0..\$#bases) { if (\$bases[\$i] =~ /[ACGT]/i) { \$cnt++; } }'"; }
		print "$cmd\n";
		$Patients{$file}{"tot_chrX_depth"} = `$cmd`;
		chomp $Patients{$file}{"tot_chrX_depth"};
		#chrY
		$cmd = "samtools view -F 0x4";
		if ($mmq) { $cmd .= " -q $mmq"; }
		$cmd .= " $file";
		foreach my$chr (@{ $fai{$file} }) {
			if ($chr =~ /^(chrY|Y)$/)
				{ $cmd .= " $chr"; }
			}
		$cmd .= " | perl -ne 'END { print \"\$cnt\" } \@fields = split (/\\t/,\$_); \@bases = split (//,\$fields[9]);";
		if ($mbq) { $cmd .= " \@quals = map { unpack(\"C*\", \$_ )} split(//,\$fields[10]); for \$i (0..\$#bases) { if ((\$bases[\$i] =~ /[ACGT]/i) && (\$quals[\$i] >= $mbq)) { \$cnt++; } }'"; }
		else { $cmd .= " for \$i (0..\$#bases) { if (\$bases[\$i] =~ /[ACGT]/i) { \$cnt++; } }'"; }
		print "$cmd\n";
		$Patients{$file}{"tot_chrY_depth"} = `$cmd`;
		chomp $Patients{$file}{"tot_chrY_depth"};
		
		$depthTxt .= $sampleName{$file}.":\n\ttotal sequenced bases for $sampleName{$file} : \n\tautoZ:".$Patients{$file}{"tot_Autosomes_depth"}."\n\tchrX:".$Patients{$file}{"tot_chrX_depth"}."\n\tchrY:".$Patients{$file}{"tot_chrY_depth"}."\n";
		}

	##gender if no sex file
	$sexTxt .= $sampleName{$file}.":\n";
	unless (exists $Patients{$file}{"Sexe"}) {
		if ($Patients{$file}{"tot_chrX_depth"} && $gonosom_Regions && $Patients{$file}{"tot_Autosomes_depth"} && $autosom_Regions) {
			$sexTxt .= "\tautoZ: ".($Patients{$file}{"tot_Autosomes_depth"}/$autosom_Regions)."\n\tchrX: ".($Patients{$file}{"tot_chrX_depth"}/$gonosom_Regions)."\n";
			if ( ($Patients{$file}{"tot_chrX_depth"}/$gonosom_Regions) > (1.2*$Patients{$file}{"tot_Autosomes_depth"}/$autosom_Regions) ) { $sexTxt .= "\t-> sexe ambigu\n"; $Patients{$file}{"Sexe"} = "F"; }
			elsif ( (($Patients{$file}{"tot_chrX_depth"}/$gonosom_Regions) <= (1.2*$Patients{$file}{"tot_Autosomes_depth"}/$autosom_Regions)) && (($Patients{$file}{"tot_chrX_depth"}/$gonosom_Regions) >= (0.8*$Patients{$file}{"tot_Autosomes_depth"}/$autosom_Regions)) ) { $Patients{$file}{"Sexe"} = "F"; }
			elsif ( (($Patients{$file}{"tot_chrX_depth"}/$gonosom_Regions) < (0.8*$Patients{$file}{"tot_Autosomes_depth"}/$autosom_Regions)) && (($Patients{$file}{"tot_chrX_depth"}/$gonosom_Regions) > (1.2*0.5*$Patients{$file}{"tot_Autosomes_depth"}/$autosom_Regions)) ) { $sexTxt .= "\t-> sexe ambigu\n"; $Patients{$file}{"Sexe"} = "F"; }
			elsif ( (($Patients{$file}{"tot_chrX_depth"}/$gonosom_Regions) <= (1.2*0.5*$Patients{$file}{"tot_Autosomes_depth"}/$autosom_Regions)) && (($Patients{$file}{"tot_chrX_depth"}/$gonosom_Regions) >= (0.8*0.5*$Patients{$file}{"tot_Autosomes_depth"}/$autosom_Regions)) ) { $Patients{$file}{"Sexe"} = "H"; }
			else { $sexTxt .= "\t-> sexe ambigu\n"; $Patients{$file}{"Sexe"} = "H"; }
			}
		else {
			$sexTxt .= "\tno data available to determine the gender\n";
			$Patients{$file}{"Sexe"} = "F"; 
			}
		}
	$sexTxt .=  "\t-> $Patients{$file}{Sexe}\n";

	##Ref_Profondeur_Patient
	if ($RefByGender) {
		if ($Patients{$file}{"Sexe"} eq "H") {
			$Patients{$file}{"Ref_depth"} = ($Patients{$file}{"tot_chrX_depth"} * 2) + $Patients{$file}{"tot_Autosomes_depth"}; } 
		else {
			$Patients{$file}{"Ref_depth"} = $Patients{$file}{"tot_chrX_depth"} + $Patients{$file}{"tot_Autosomes_depth"}; }
		}
	else {
		$Patients{$file}{"Ref_depth"} = $Patients{$file}{"tot_chrX_depth"} + $Patients{$file}{"tot_chrY_depth"} + $Patients{$file}{"tot_Autosomes_depth"};
		}

	}
if ($depthTxt) { print "\nbases counts:\n$depthTxt\n"; }
print "\ngender:\n$sexTxt\n";

my$meanRef;
foreach my$file (@Files) {
	if ($Patients{$file}{"Ref_depth"})
		{ $meanRef += $Patients{$file}{"Ref_depth"}; }
	else { $Patients{$file}{"ecarte"} = 1; }
	}
$meanRef /= scalar@Files;
foreach my$file (@Files) {
	$Patients{$file}{"Ref_depth"} /= $meanRef; 
	print "normalization factor for $sampleName{$file} : ".$Patients{$file}{"Ref_depth"}."\n";
	}



##iterations
my $nb_parcours = 0;
my $continuer = 1;
my %Results;
while ($continuer == 1 || $nb_parcours <= 1) {

	my $iteration = $nb_parcours + 1;
	print "Iteration numero \: ".$iteration."\n";

	undef %Results;
	my %CNV;
	my $nb_regions_conservees = 0;
	my $regions_ecartees = 0;
	$continuer = 0;

	my $log = "$outdir/Logs_Iteration_".$iteration."\.txt";
	open(LOG,">",$log) or die("Pb lors de l'ecriture du fichier log $!\n");

	my $sortie = "$outdir/CNV_Iteration_".$iteration."\.tab";
	open(SORTIE,">",$sortie) or die("Pb lors de l'ecriture du fichier sortie $!\n");
	print SORTIE "Chrom"."\t";
	print SORTIE "Start"."\t";
	print SORTIE "End"."\t";
	if ($nCol>3) {
		print SORTIE "Infos";
		for (0..($nCol-4)) { print SORTIE "\t"; }
		}
	print SORTIE "Numero_Region"."\t";
	foreach my$file (@Files)
		{ print SORTIE $sampleName{$file}."\t"; }
	print SORTIE "Statut_Region\n";


	# PREMIER PARCOURS BIS #
	# RECALCUL DE LA PROFONDEUR TOTALE POUR LA PONDERATION INTRA
	# CE (RE)CALCUL A LIEU SI DES REGIONS ETAIENT "MOCHES" APRES L'APPEL DE CNV
	if ($nb_parcours > 0) {

		if ($Ref eq "mean") {
		# REINITIALISATION DE LA PROFONDEUR TOTALE PAR PATIENT
			foreach my$file (@Files) {
				$Patients{$file}{"tot_Autosomes_depth"} = 0;
				$Patients{$file}{"Ref_Gonosomes_Patient"} = 0;
				}
			for (my $r = 0 ; $r < $#Regions ; $r++) {
				if(!defined $Regions[$r]{"Appel"}) {
					foreach my$file (@Files) {
						if ($Regions[$r]{"Chrom"} !~ /^(chr[XY]|[XY])$/) {
							$Patients{$file}{"tot_Autosomes_depth"} += $Regions[$r]{$file}{"Mean"};
							}
						elsif ($Regions[$r]{"Chrom"} =~ /^(chrX|X)$/) {
							$Patients{$file}{"Ref_Gonosomes_Patient"} += $Regions[$r]{$file}{"Mean"};
							}
						else {
							if ((!$RefByGender) && ($Patients{$file}{"Sexe"} eq "H")) {
								$Patients{$file}{"Ref_Gonosomes_Patient"} += $Regions[$r]{$file}{"Mean"};
								}
							}
						}
					}
				}
			
			foreach my$file (@Files) {
				if ($RefByGender && ($Patients{$file}{"Sexe"} eq "H")) {
					$Patients{$file}{"Ref_depth"} = ($Patients{$file}{"Ref_Gonosomes_Patient"} * 2) + $Patients{$file}{"tot_Autosomes_depth"}; } 
				else {
					$Patients{$file}{"Ref_depth"} = $Patients{$file}{"Ref_Gonosomes_Patient"} + $Patients{$file}{"tot_Autosomes_depth"}; }
				}
			foreach my$file (@Files)
				{ $meanRef += $Patients{$file}{"Ref_depth"}; }
			$meanRef /= scalar@Files;
			foreach my$file (@Files) {
				$Patients{$file}{"Ref_depth"} /= $meanRef;
				}
			}
		}
	print LOG "sample normalization factors :\n";
	foreach my$file (@Files) {
		print LOG "\t$sampleName{$file} : $Patients{$file}{Ref_depth}\n";
	}
	print LOG "\n";


	# SECOND PARCOURS DES REGIONS #
	# PERMET DE PONDERER LA PROFONDEUR PAR LES AUTRES REGIONS (INTRA) ET ENTRE LES PATIENTS (INTER)
	for my$r (0..$#Regions) {

		print SORTIE $Regions[$r]{"allLine"}."\t";
		print SORTIE $r."\t";

		# SI LA REGION N'EST PAS "MOCHE"
		if (!(defined($Regions[$r]{"Appel"}))) {

			$nb_regions_conservees++;
			$Regions[$r]{"nb_CNV"}{"DEL"} = 0;
			$Regions[$r]{"nb_CNV"}{"DUP"} = 0;
			my $nb_evts = 0;
			my$prof_Moyenne_Inter = "";

			if ($normByGender) {

				@{ $Regions[$r]{"all_normByS_depths_fem"} } = ();
				@{ $Regions[$r]{"all_normByS_depths_males"} } = ();
				# Pour chaque patient, un premier parcours permet la pondération intra-patient (en divisant la profondeur de la region par la profondeur totale du patient)
				foreach my$file (@Files) {

					if ($Patients{$file}{"Sexe"} eq "F" && !($Patients{$file}{"ecarte"})) {

						# Le controle permet de diviser uniquement par le nombre de patients restant apres l filtre, et non par l'ensemble des patients
						# Si la region concerne un gonosome et le patient est une femme, la référence est celle qui correspond aux autosomes + au chromosome X
						$Regions[$r]{$file}{"normByS_depth"} = $Regions[$r]{$file}{"Mean"} / $Patients{$file}{"Ref_depth"};
						# On sauvegarde la profondeur ponderee intra patient, pour determiner ensuite la moyenne de ces profondeurs ponderees
						push( @{ $Regions[$r]{"all_normByS_depths_fem"} }, $Regions[$r]{$file}{"normByS_depth"} );


					} elsif ($Patients{$file}{"Sexe"} eq "H" && !($Patients{$file}{"ecarte"})) {

						# Le controle permet de divisier uniquement par le nombre de patients restant apres l filtre, et non par l'ensemble des patients
						# Si la region concerne un gonosome et le patient est un homme, la référence est celle qui correspond uniquement au chromosome X
						$Regions[$r]{$file}{"normByS_depth"} = $Regions[$r]{$file}{"Mean"} / $Patients{$file}{"Ref_depth"};
						# On sauvegarde la profondeur ponderee intra patient, pour determiner ensuite la moyenne de ces profondeurs ponderees
						push( @{ $Regions[$r]{"all_normByS_depths_males"} }, $Regions[$r]{$file}{"normByS_depth"} );

					}

				}

				# Nous calculons pour chaque region et chaque sexe la moyenne de profondeur (ponderee pour chaque patient par rapport aux autres regions)
				if ($normByGender eq "all") {
					if(@{ $Regions[$r]{"all_normByS_depths_fem"} }) {
						foreach (@cnvVal) 
							{ $Regions[$r]{"normByR_depth_fem"}{$_} = norm($_,\@{ $Regions[$r]{"all_normByS_depths_fem"} }); }
						if ($norm eq "std") {
							my$sqtotal = 0;
							foreach (@{ $Regions[$r]{"all_normByS_depths_fem"} })
								{ $sqtotal += ($Regions[$r]{"normByR_depth_fem"}{"moy"}-$_)**2; }
							$Regions[$r]{"normByR_depth_fem"}{"std"} = ($sqtotal / (scalar@{ $Regions[$r]{"all_normByS_depths_fem"} }-1))**0.5;
						}
					}
					if(@{ $Regions[$r]{"all_normByS_depths_males"} }) {
						foreach (@cnvVal)
							{ $Regions[$r]{"normByR_depth_males"}{$_} = norm({$_},\@{ $Regions[$r]{"all_normByS_depths_males"} }); }
						if ($norm eq "std") {
							my$sqtotal = 0;
							foreach (@{ $Regions[$r]{"all_normByS_depths_males"} })
								{ $sqtotal += ($Regions[$r]{"normByR_depth_males"}{"moy"}-$_)**2; }
							$Regions[$r]{"normByR_depth_males"}{"std"} = ($sqtotal / (scalar@{ $Regions[$r]{"all_normByS_depths_males"} }-1))**0.5;
						}
					}


				} else {
					if ($Regions[$r]{"Chrom"} !~ /^(chr[XY]|[XY])$/) {
						my@Autosomes=();
						if(@{ $Regions[$r]{"all_normByS_depths_fem"} }) {
							push(@Autosomes, @{ $Regions[$r]{"all_normByS_depths_fem"} });
						}
						if(@{ $Regions[$r]{"all_normByS_depths_males"} }) {
							push(@Autosomes, @{ $Regions[$r]{"all_normByS_depths_males"} });
						}
						if (@Autosomes) {
							foreach (@cnvVal) {
								$Regions[$r]{"normByR_depth"}{$_} = norm($_,\@Autosomes);
								}
							if ($norm eq "std") {
								my$sqtotal = 0;
								foreach (@Autosomes)
									{ $sqtotal += ($Regions[$r]{"normByR_depth"}{"moy"}-$_)**2; }
								$Regions[$r]{"normByR_depth"}{"std"} = ($sqtotal / (scalar@Autosomes-1))**0.5;
							}
						}
					} else {
						if(@{ $Regions[$r]{"all_normByS_depths_fem"} }) {
							foreach (@cnvVal) {
								$Regions[$r]{"normByR_depth_fem"}{$_} = norm($_,\@{ $Regions[$r]{"all_normByS_depths_fem"} });
							}
							if ($norm eq "std") {
								my$sqtotal = 0;
								foreach (@{ $Regions[$r]{"all_normByS_depths_fem"} })
									{ $sqtotal += ($Regions[$r]{"normByR_depth_fem"}{"moy"}-$_)**2; }
								$Regions[$r]{"normByR_depth_fem"}{"std"} = ($sqtotal / (scalar@{ $Regions[$r]{"all_normByS_depths_fem"} }-1))**0.5;
							}
						}
						if(@{ $Regions[$r]{"all_normByS_depths_males"} }) {
							foreach (@cnvVal) {
								$Regions[$r]{"normByR_depth_males"}{$_} = norm($_,\@{ $Regions[$r]{"all_normByS_depths_males"} });
							}
							if ($norm eq "std") {
								my$sqtotal = 0;
								foreach (@{ $Regions[$r]{"all_normByS_depths_males"} })
									{ $sqtotal += ($Regions[$r]{"normByR_depth_males"}{"moy"}-$_)**2; }
								$Regions[$r]{"normByR_depth_males"}{"std"} = ($sqtotal / (scalar@{ $Regions[$r]{"all_normByS_depths_males"} }-1))**0.5;
							}
						}
					}
				}

			
				# Un second parcours permet l'appel de CNV (en comparant la profondeur pour un patient par la profondeur moyenne de la region)
				foreach my$file (@Files) {

					if ($norm eq "std") { $prof_Moyenne_Inter = 0; }
					else  { $prof_Moyenne_Inter = 1; }

					if(! $Patients{$file}{"ecarte"}) {

						if ($normByGender eq "all") {

							if($Patients{$file}{"Sexe"} eq "F") {
								if ($Regions[$r]{"normByR_depth_fem"}{$norm}) {
									if ($Regions[$r]{"Chrom"} !~ /^(chrY|Y)$/) {
										if ($norm eq "std") {
											$Regions[$r]{$file}{"depth_ratio"}{"moy"} = $Regions[$r]{$file}{"normByS_depth"}/$Regions[$r]{"normByR_depth_fem"}{"moy"};
											$Regions[$r]{$file}{"depth_ratio"}{"std"} = ($Regions[$r]{$file}{"normByS_depth"}-$Regions[$r]{"normByR_depth_fem"}{"moy"})/$Regions[$r]{"normByR_depth_fem"}{"std"};
										} else {
											$Regions[$r]{$file}{"depth_ratio"}{$norm} = $Regions[$r]{$file}{"normByS_depth"}/$Regions[$r]{"normByR_depth_fem"}{$norm};
										}
										$prof_Moyenne_Inter = $Regions[$r]{$file}{"depth_ratio"}{$norm};
										print SORTIE "$prof_Moyenne_Inter\t";
									} else {
										# En faisant la distinction entre chromosomes X et Y, les femmes ne sont plus considérées pour l'appel de CNV des régions du dernier
										print SORTIE "NA\t";
									}
								} else {
									$Regions[$r]{"Appel"} = "No_Data";
									print SORTIE "NA\t";
								}

							} elsif ($Patients{$file}{"Sexe"} eq "H") {
								if ($Regions[$r]{"normByR_depth_males"}{$norm}) {
									if ($norm eq "std") {
										$Regions[$r]{$file}{"depth_ratio"}{"moy"} = $Regions[$r]{$file}{"normByS_depth"}/$Regions[$r]{"normByR_depth_males"}{"moy"};
										$Regions[$r]{$file}{"depth_ratio"}{"std"} = ($Regions[$r]{$file}{"normByS_depth"}-$Regions[$r]{"normByR_depth_males"}{"moy"})/$Regions[$r]{"normByR_depth_males"}{"std"};
									} else {
										$Regions[$r]{$file}{"depth_ratio"}{$norm} = $Regions[$r]{$file}{"normByS_depth"}/$Regions[$r]{"normByR_depth_males"}{$norm};
									}
									$prof_Moyenne_Inter = $Regions[$r]{$file}{"depth_ratio"}{$norm};
									print SORTIE "$prof_Moyenne_Inter\t";
								} else {
									$Regions[$r]{"Appel"} = "No_Data";
									print SORTIE "NA\t";
								}
							}

						} else {
							# Pour les autosomes
							if ($Regions[$r]{"Chrom"} !~ /^(chr[XY]|[XY])$/) {
								if ($Regions[$r]{"normByR_depth"}{$norm}) {
									if ($norm eq "std") {
										$Regions[$r]{$file}{"depth_ratio"}{"moy"} = $Regions[$r]{$file}{"normByS_depth"}/$Regions[$r]{"normByR_depth"}{"moy"};
										$Regions[$r]{$file}{"depth_ratio"}{"std"} = ($Regions[$r]{$file}{"normByS_depth"}-$Regions[$r]{"normByR_depth"}{"moy"})/$Regions[$r]{"normByR_depth"}{"std"};
									} else {
										$Regions[$r]{$file}{"depth_ratio"}{$norm} = $Regions[$r]{$file}{"normByS_depth"}/$Regions[$r]{"normByR_depth"}{$norm};
									}
									$prof_Moyenne_Inter = $Regions[$r]{$file}{"depth_ratio"}{$norm};
									print SORTIE "$prof_Moyenne_Inter\t";
								} else {
									$Regions[$r]{"Appel"} = "No_Data";
									print SORTIE "NA\t";
								}
							} else {
							# Pour les gonosomes
								if($Patients{$file}{"Sexe"} eq "F") {
									if ($Regions[$r]{"normByR_depth_fem"}{$norm}) {
										if ($Regions[$r]{"Chrom"} =~ /^(chrX|X)$/) {
											if ($norm eq "std") {
												$Regions[$r]{$file}{"depth_ratio"}{"moy"} = $Regions[$r]{$file}{"normByS_depth"}/$Regions[$r]{"normByR_depth_fem"}{"moy"};
												$Regions[$r]{$file}{"depth_ratio"}{"std"} = ($Regions[$r]{$file}{"normByS_depth"}-$Regions[$r]{"normByR_depth_fem"}{"moy"})/$Regions[$r]{"normByR_depth_fem"}{"std"};
											} else {
												$Regions[$r]{$file}{"depth_ratio"}{$norm} = $Regions[$r]{$file}{"normByS_depth"}/$Regions[$r]{"normByR_depth_fem"}{$norm};
											}
											$prof_Moyenne_Inter = $Regions[$r]{$file}{"depth_ratio"}{$norm};
											print SORTIE "$prof_Moyenne_Inter\t";
										} else {
											# En faisant la distinction entre chromosomes X et Y, les femmes ne sont plus considérées pour l'appel de CNV des régions du dernier
											print SORTIE "NA\t";
										}
									} else {
										$Regions[$r]{"Appel"} = "No_Data";
										print SORTIE "NA\t";
									}
								} elsif ($Patients{$file}{"Sexe"} eq "H") {
									if ($Regions[$r]{"normByR_depth_males"}{$norm}) {
										if ($norm eq "std") {
											$Regions[$r]{$file}{"depth_ratio"}{"moy"} = $Regions[$r]{$file}{"normByS_depth"}/$Regions[$r]{"normByR_depth_males"}{"moy"};
											$Regions[$r]{$file}{"depth_ratio"}{"std"} = ($Regions[$r]{$file}{"normByS_depth"}-$Regions[$r]{"normByR_depth_males"}{"moy"})/$Regions[$r]{"normByR_depth_males"}{"std"};
										} else {
											$Regions[$r]{$file}{"depth_ratio"}{$norm} = $Regions[$r]{$file}{"normByS_depth"}/$Regions[$r]{"normByR_depth_males"}{$norm};
										}
										$prof_Moyenne_Inter = $Regions[$r]{$file}{"depth_ratio"}{$norm};
										print SORTIE "$prof_Moyenne_Inter\t";
									} else {
										$Regions[$r]{"Appel"} = "No_Data";
										print SORTIE "NA\t";
									}
								}
							}
						}
			
						#print $Regions[$r]{$patient}{"depth_ratio"}."\t";
						if ($prof_Moyenne_Inter < $seuil_deletion) {
							$Regions[$r]{"nb_CNV"}{"DEL"}++;
							$nb_evts++;
							$CNV{$file}++;
							$Results{$file}{$r} = "DEL";
						} elsif ($prof_Moyenne_Inter > $seuil_duplication) {
							$Regions[$r]{"nb_CNV"}{"DUP"}++;
							$nb_evts++;
							$CNV{$file}++;
							$Results{$file}{$r} = "DUP";
						}

					} else {
						print SORTIE "NA\t";
					}
				}

				my$recurrent="";
				if(@{$Regions[$r]{"all_normByS_depths_fem"}} && @{$Regions[$r]{"all_normByS_depths_males"}} && (($nb_evts/(scalar@{$Regions[$r]{"all_normByS_depths_fem"}} + scalar@{$Regions[$r]{"all_normByS_depths_males"}})) > $seuil_region) && ($nb_parcours > 0)) { $recurrent = 1; }
				elsif(@{$Regions[$r]{"all_normByS_depths_fem"}} && !@{$Regions[$r]{"all_normByS_depths_males"}} && (($nb_evts/scalar@{$Regions[$r]{"all_normByS_depths_fem"}}) > $seuil_region) && ($nb_parcours > 0)) { $recurrent = 1; }
				elsif(!@{$Regions[$r]{"all_normByS_depths_fem"}} && @{$Regions[$r]{"all_normByS_depths_males"}} && (($nb_evts/scalar@{$Regions[$r]{"all_normByS_depths_males"}}) > $seuil_region) && ($nb_parcours > 0)) { $recurrent = 1; }
				if ($recurrent) {
					print SORTIE "CNV_Recurrent";
					print LOG "Region ecartee \: ". $Regions[$r]{"allLine"}."\n";
					$Regions[$r]{"Appel"} = "CNV_Recurrent";
					$regions_ecartees++;
					$continuer = 1;

				} else {		
					print SORTIE "OK";
				}
				print SORTIE "\n";


			} else {

				@{ $Regions[$r]{"all_normByS_depths"} }=();
				# Pour chaque patient, un premier parcours permet la pondération intra-patient (en divisant la profondeur de la region par la profondeur totale recalculée du patient)
				foreach my$file (@Files) {
					if(! $Patients{$file}{"ecarte"}) {
						$Regions[$r]{$file}{"normByS_depth"} = $Regions[$r]{$file}{"Mean"} / $Patients{$file}{"Ref_depth"};
						if ( ($Patients{$file}{"Sexe"} eq "H") && ($Regions[$r]{"Chrom"} =~ /^(chrX|X)$/) )
							{ push(@{ $Regions[$r]{"all_normByS_depths"} }, ($Regions[$r]{$file}{"normByS_depth"}*2)); }
						elsif ( ($Patients{$file}{"Sexe"} eq "F") && ($Regions[$r]{"Chrom"} =~ /^(chrY|Y)$/) )
							{ next ; }
						else { push(@{ $Regions[$r]{"all_normByS_depths"} }, $Regions[$r]{$file}{"normByS_depth"}); }
					}
				}
				# Nous calculons pour chaque region la mediane/moyenne de profondeur (ponderee pour chaque patient par rapport aux autres regions)
				if(@{ $Regions[$r]{"all_normByS_depths"} }) {
					foreach (@cnvVal) { 
						$Regions[$r]{"normByR_depth"}{$_} = norm($_,\@{ $Regions[$r]{"all_normByS_depths"} });
						}
					if ($norm eq "std") {
						my$sqtotal = 0;
						foreach (@{ $Regions[$r]{"all_normByS_depths"} })
							{ $sqtotal += ($Regions[$r]{"normByR_depth"}{"moy"}-$_)**2; }
						$Regions[$r]{"normByR_depth"}{"std"} = ($sqtotal / (scalar@{ $Regions[$r]{"all_normByS_depths"} }-1))**0.5;
					}
				}


				# Un second parcours permet l'appel de CNV (en comparant la profondeur pour un patient par la profondeur moyenne de la region)
				foreach my$file (@Files) {

					if ($norm eq "std") { $prof_Moyenne_Inter = 0; }
					else  { $prof_Moyenne_Inter = 1; }

					if(! $Patients{$file}{"ecarte"}) {

						if ($Regions[$r]{"normByR_depth"}{$norm}) {
							if($Patients{$file}{"Sexe"} eq "F") {
								if ($Regions[$r]{"Chrom"} !~ /^(chrY|Y)$/) {
									if ($norm eq "std") {
										$Regions[$r]{$file}{"depth_ratio"}{"moy"} = $Regions[$r]{$file}{"normByS_depth"}/$Regions[$r]{"normByR_depth"}{"moy"};
										$Regions[$r]{$file}{"depth_ratio"}{"std"} = ($Regions[$r]{$file}{"normByS_depth"}-$Regions[$r]{"normByR_depth"}{"moy"})/$Regions[$r]{"normByR_depth"}{"std"};
									} else {
										$Regions[$r]{$file}{"depth_ratio"}{$norm} = $Regions[$r]{$file}{"normByS_depth"}/$Regions[$r]{"normByR_depth"}{$norm};
									}
									$prof_Moyenne_Inter = $Regions[$r]{$file}{"depth_ratio"}{$norm};
									print SORTIE "$prof_Moyenne_Inter\t";
								} else {
									#les femmes ne sont pas considérées pour l'appel de CNV des régions du Y
									print SORTIE "NA\t";
								}
							} else {
								#*2 for chrX
								if ($Regions[$r]{"Chrom"} =~ /^(chrX|X)$/) {
									if ($norm eq "std") {
										$Regions[$r]{$file}{"depth_ratio"}{"moy"} = $Regions[$r]{$file}{"normByS_depth"}*2/$Regions[$r]{"normByR_depth"}{"moy"};
										$Regions[$r]{$file}{"depth_ratio"}{"std"} = ($Regions[$r]{$file}{"normByS_depth"}*2-$Regions[$r]{"normByR_depth"}{"moy"})/$Regions[$r]{"normByR_depth"}{"std"};
									} else {
										$Regions[$r]{$file}{"depth_ratio"}{$norm} = $Regions[$r]{$file}{"normByS_depth"}*2/$Regions[$r]{"normByR_depth"}{$norm};
									}
								} else {
									if ($norm eq "std") {
										$Regions[$r]{$file}{"depth_ratio"}{"moy"} = $Regions[$r]{$file}{"normByS_depth"}/$Regions[$r]{"normByR_depth"}{"moy"};
										$Regions[$r]{$file}{"depth_ratio"}{"std"} = ($Regions[$r]{$file}{"normByS_depth"}-$Regions[$r]{"normByR_depth"}{"moy"})/$Regions[$r]{"normByR_depth"}{"std"};
									} else {
										$Regions[$r]{$file}{"depth_ratio"}{$norm} = $Regions[$r]{$file}{"normByS_depth"}/$Regions[$r]{"normByR_depth"}{$norm};
									}
								}
								$prof_Moyenne_Inter = $Regions[$r]{$file}{"depth_ratio"}{$norm};
								print SORTIE "$prof_Moyenne_Inter\t";
							}
							if ($prof_Moyenne_Inter < $seuil_deletion) {
								$Regions[$r]{"nb_CNV"}{"DEL"}++;
								$nb_evts++;
								$CNV{$file}++;
								$Results{$file}{$r} = "DEL";
							} elsif ($prof_Moyenne_Inter > $seuil_duplication) {
								$Regions[$r]{"nb_CNV"}{"DUP"}++;
								$nb_evts++;
								$CNV{$file}++;
								$Results{$file}{$r} = "DUP";
							}
						} else {
							$Regions[$r]{"Appel"} = "No_Data";
							#$prof_Moyenne_Inter = 1;
							print SORTIE "NA\t";
						}
					} else {
						print SORTIE "NA\t";
					}

				}

				if(@{ $Regions[$r]{"all_normByS_depths"} } && (($nb_evts/scalar@{ $Regions[$r]{"all_normByS_depths"} }) > $seuil_region) && ($nb_parcours > 0)) {
					print SORTIE "CNV_Recurrent";
					print LOG "Region ecartee \: ". $Regions[$r]{"allLine"}."\n";
					$Regions[$r]{"Appel"} = "CNV_Recurrent";
					$regions_ecartees++;
					$continuer = 1;

				} else {		
					print SORTIE "OK";
				}
				print SORTIE "\n";
			}

		# SI LA REGION EST "MOCHE"
		} else {
			foreach my$file (@Files) {
				print SORTIE "NA\t";
			}
			print SORTIE $Regions[$r]{"Appel"}."\n";
		}

	}

	print LOG "Nombre de regions ecartees lors de cette iteration \: ".$regions_ecartees."\n\n\n";	

	my $patients_ecartes = 0;

	foreach my$file (@Files) {

		if (! $Patients{$file}{"ecarte"}) {
			my $prcent_evt;
			if (defined($CNV{$file})) {

				$prcent_evt = $CNV{$file}/$nb_regions_conservees;

				if ($prcent_evt >= $seuil_patient) {
					print LOG "Patient ecarte \: ".$sampleName{$file}."\n";
					$Patients{$file}{"ecarte"} = 1;
					$continuer = 1;
					$patients_ecartes++;			
				} else {
					$Patients{$file}{"ecarte"} = 0;
				}


			} else {
				print $sampleName{$file}." \: aucun CNV identifie\n";
				$Patients{$file}{"ecarte"} = 0;
			}


		}

	}

	print LOG "Nombre de patients ecartes lors de cette iteration \: ".$patients_ecartes."\n";

	#print "\n";
	close(SORTIE);
	close(LOG);
	$nb_parcours++;
	
}


for my$file (@Files) {
	mkdir "$outdir/$sampleName{$file}";
	}


##print CNV foreach sample

## $Results{$file}{$r} = "DEL";
## @{ $RegionOrder{$Chrom} } = [ regions sorted by pos ]
my(%uniqReg,%RegInArray,%regionOrder,%regionIndice);
for my$r (0..$#Regions) {
	unless (exists $uniqReg{ $Regions[$r]{"Chrom"}."-".$Regions[$r]{"Start"}."-".$Regions[$r]{"End"} }) {
		push (@{ $RegInArray{ $Regions[$r]{"Chrom"} } }, $r);
		$uniqReg{ $Regions[$r]{"Chrom"}."-".$Regions[$r]{"Start"}."-".$Regions[$r]{"End"} } = 1;
		}
	my@tab = split(/\t/, $Regions[$r]{"allLine"});
	if ($tab[3]) {
		my@tab2 = split(/:|,/, $tab[3]);
		$Regions[$r]{"label"} = substr($tab2[0], 0, 25);
		$Regions[$r]{"Gene"} = $tab2[0];
		}
	else { $Regions[$r]{"label"} = "$tab[1]-$tab[2]"; }
	}
foreach my$Chrom (keys%RegInArray) {
	##sort by region ends (in case several regions with same start) then by starts
	@{ $regionOrder{$Chrom} } = sort{$Regions[$a]{"End"}<=>$Regions[$b]{"End"}}@{ $RegInArray{$Chrom} }; 
	@{ $regionOrder{$Chrom} } = sort{$Regions[$a]{"Start"}<=>$Regions[$b]{"Start"}}@{ $RegInArray{$Chrom} }; 
	for (my$i=0;$i<scalar@{ $regionOrder{$Chrom} };$i++) { $regionIndice{ $regionOrder{$Chrom}[$i] } = $i; }
	}

## @{ $orderedCNV{$patient}{$Chrom} } = [r1, r2,...]
my%orderedCNV;
foreach my$file (keys%Results) {
	my(%uniqCNV,%CNVinArray);
	foreach my$r (sort{$a<=>$b}keys%{ $Results{$file} }) {
		if (!exists $uniqCNV{$Regions[$r]{"Chrom"}."-".$Regions[$r]{"Start"}."-".$Regions[$r]{"End"}}) { 
			push (@{ $CNVinArray{ $Regions[$r]{"Chrom"} } }, $r);
			$uniqCNV{ $Regions[$r]{"Chrom"}."-".$Regions[$r]{"Start"}."-".$Regions[$r]{"End"} } = 1;
			}
		}
	foreach my$Chrom (keys%CNVinArray) {
		##sort by region ends (in case several regions with same start) then by starts
		@{ $orderedCNV{$file}{$Chrom} } = sort{$Regions[$a]{"End"}<=>$Regions[$b]{"End"}}@{ $CNVinArray{$Chrom} }; 
		@{ $orderedCNV{$file}{$Chrom} } = sort{$Regions[$a]{"Start"}<=>$Regions[$b]{"Start"}}@{ $CNVinArray{$Chrom} }; 
		}
	}

unless ($minCNV) { $minCNV = 1; }
my(%Result2,%Result3,%Result4);
foreach my$file (keys%Results) {
	open (CNV1,">$outdir/$sampleName{$file}/CNV_$sampleName{$file}.summary.txt") or die "Pb lors de l'ecriture du fichier sortie $!\n";
	print CNV1 "#Region\tCNV\tclean_intervals_Nbr\tclean_ratio_to_$norm\tdirty_intervals_Nbr\tdirty_ratio_to_$norm\toverlapping_Samples\n";
	open (CNV2,">$outdir/$sampleName{$file}/CNV_$sampleName{$file}.allIntervals.txt") or die "Pb lors de l'ecriture du fichier sortie $!\n";
	print CNV2 "#Chrom\tStart\tEnd\tLength\tInfo\tinterval_Order\tCNV_type\tratio_to_$norm\toccurences";
	if (@cnvFields) {
		foreach (@cnvFields) { 
			if ($_ eq "norm") {
				if ($norm eq "std") { print CNV2 "\tmoy\tstd"; }
				else  { print CNV2 "\t$norm"; }
				}
			else { print CNV2 "\t$_"; }
			}
		}
	print CNV2 "\n";
	foreach my$Chromosome (@ChromOrder) {
		my$Chrom = $Chromosome ; $Chrom =~ s/chr//i;
		if (exists $orderedCNV{$file}{$Chrom}) {
			my$r=0;		##index in @{ $orderedCNV{$file}{$Chrom} }
			while ($r < scalar@{ $orderedCNV{$file}{$Chrom} } ) {
				##merge consecutive CNVs
				my@refSub = mergeConsecutiveCNV("",$maxNonCNV,$orderedCNV{$file}{$Chrom}[$r],$regionIndice{$orderedCNV{$file}{$Chrom}[$r]},\@{ $regionOrder{$Chrom} },\%{ $Results{$file} },\@Regions);
				my$i=$refSub[0]; my$cnvOK=$refSub[1]; my$nonCNVtot=$refSub[2];
				my@nextReg = @{ $refSub[3] };
				if ($maxNonCNVrate) {
					while (($nonCNVtot/($cnvOK+1)) > $maxNonCNVrate) {
						@refSub = mergeConsecutiveCNV(($cnvOK-1),$maxNonCNV,$orderedCNV{$file}{$Chrom}[$r],$regionIndice{$orderedCNV{$file}{$Chrom}[$r]},\@{ $regionOrder{$Chrom} },\%{ $Results{$file} },\@Regions);
						$i=$refSub[0]; my$cnvOK=$refSub[1]; my$nonCNVtot=$refSub[2];
						@nextReg = @{ $refSub[3] };
						}
					}
				##overlapping samples
				my$overlapCNT=0; my%overlapSMPL=();
				for (my$j=0;$j<=$cnvOK;$j++) {
					foreach my$other (keys%Patients) {
						if ( (exists $Results{$other}{$nextReg[$j]}) && ($Results{$other}{$nextReg[$j]} eq $Results{$file}{$nextReg[0]}) && (!exists $overlapSMPL{$other}) ) { 
							$overlapCNT++;
							$overlapSMPL{$other}=1;
							}
						}
					}
				##print patient's summary and allIntervals
				if ($cnvOK >= ($minCNV-1)) {
					if ($cnvOK == 0) {
						$Result2{$file}{$nextReg[0]} = $Results{$file}{$nextReg[0]};
						$Result3{$file}{$nextReg[0]} = $Results{$file}{$nextReg[0]};
						$Result4{$file}{$Chrom}{$Regions[$nextReg[0]]{"Start"}}{$Regions[$nextReg[0]]{"End"}} = $Results{$file}{$nextReg[0]};
						if (exists $regionIndice{$nextReg[0]}) {
							print CNV1 $Chromosome.":".$Regions[$nextReg[0]]{"Start"}."-".$Regions[$nextReg[0]]{"End"}."\t".$Results{$file}{$nextReg[0]}."\t1\t".sprintf("%.3f",$Regions[$nextReg[0]]{$file}{"depth_ratio"}{$norm})."\t.\t.\t$overlapCNT\n";
							print CNV2 $Chromosome."\t".$Regions[$nextReg[0]]{"Start"}."\t".$Regions[$nextReg[0]]{"End"}."\t".($Regions[$nextReg[0]]{"End"}-$Regions[$nextReg[0]]{"Start"}+1)." bp\t".$Regions[$nextReg[0]]{"label"}."\t".($regionIndice{$nextReg[0]}+1)."\t".$Results{$file}{$nextReg[0]}."\t".sprintf("%.3f",$Regions[$nextReg[0]]{$file}{"depth_ratio"}{$norm})."\t".$Regions[$nextReg[0]]{"nb_CNV"}{$Results{$file}{$nextReg[0]}};
							if (@cnvFields) {
								my$txt = printCNVfields($normByGender,\@cnvFields,\%{ $Regions[$nextReg[0]] },\%{ $Patients{$file} });
								print CNV2 "$txt";
								}
							print CNV2 "\n\n";
							}
						}
					else {
						my@cleanCNV = (); my@dirtyCNV = ();
						for (my$j=0;$j<=$cnvOK;$j++) {
							if (exists $Results{$file}{$nextReg[$j]}) {
								$Result2{$file}{$nextReg[$j]} = $Results{$file}{$nextReg[$j]};
								$Result4{$file}{$Chrom}{$Regions[$nextReg[$j]]{"Start"}}{$Regions[$nextReg[$j]]{"End"}} = $Results{$file}{$nextReg[$j]};
								push(@cleanCNV, $Regions[$nextReg[$j]]{$file}{"depth_ratio"}{$norm});
								push(@dirtyCNV, $Regions[$nextReg[$j]]{$file}{"depth_ratio"}{$norm});
								if (exists $regionIndice{$nextReg[$j]}) {
									print CNV2 $Chromosome."\t".$Regions[$nextReg[$j]]{"Start"}."\t".$Regions[$nextReg[$j]]{"End"}."\t".($Regions[$nextReg[$j]]{"End"}-$Regions[$nextReg[$j]]{"Start"}+1)." bp\t".$Regions[$nextReg[$j]]{"label"}."\t".($regionIndice{$nextReg[$j]}+1)."\t".$Results{$file}{$nextReg[$j]}."\t".sprintf("%.3f",$Regions[$nextReg[$j]]{$file}{"depth_ratio"}{$norm})."\t".$Regions[$nextReg[$j]]{"nb_CNV"}{$Results{$file}{$nextReg[$j]}};
									if (@cnvFields) {
										my$txt = printCNVfields($normByGender,\@cnvFields,\%{ $Regions[$nextReg[$j]] },\%{ $Patients{$file} });
										print CNV2 "$txt";
										}
									print CNV2 "\n";
									}
								}
							else {
								$Result2{$file}{$nextReg[$j]} = "NA";
								$Result4{$file}{$Chrom}{$Regions[$nextReg[$j]]{"Start"}}{$Regions[$nextReg[$j]]{"End"}} = "NA";
								if (exists $regionIndice{$nextReg[$j]}) {
									print CNV2 $Chromosome."\t".$Regions[$nextReg[$j]]{"Start"}."\t".$Regions[$nextReg[$j]]{"End"}."\t".($Regions[$nextReg[$j]]{"End"}-$Regions[$nextReg[$j]]{"Start"}+1)." bp\t".$Regions[$nextReg[$j]]{"label"}."\t".($regionIndice{$nextReg[$j]}+1)."\t";
									if (exists $Regions[$nextReg[$j]]{"Appel"}) { print CNV2 "NA\t"; }
									else { print CNV2 "no\t"; }
									if (exists $Regions[$nextReg[$j]]{$file}{"depth_ratio"}{$norm}) {
										print CNV2 sprintf("%.3f",$Regions[$nextReg[$j]]{$file}{"depth_ratio"}{$norm})."\t".$Regions[$nextReg[$j]]{"nb_CNV"}{$Results{$file}{$nextReg[0]}};
										push(@dirtyCNV, $Regions[$nextReg[$j]]{$file}{"depth_ratio"}{$norm});
										}
									else { print CNV2 "na\tna"; }
									if (@cnvFields) {
										my$txt = printCNVfields($normByGender,\@cnvFields,\%{ $Regions[$nextReg[$j]] },\%{ $Patients{$file} });
										print CNV2 "$txt";
										}
									print CNV2 "\n";
									}
								}
							}
						print CNV2 "\n";
						my$cleanAverage = 0;
						foreach (@cleanCNV) { $cleanAverage += $_; }
						$cleanAverage /= scalar@cleanCNV;
						my$dirtyAverage = 0;
						foreach (@dirtyCNV) { $dirtyAverage += $_; }
						$dirtyAverage /= scalar@dirtyCNV;
						print CNV1 $Chromosome.":".$Regions[$nextReg[0]]{"Start"}."-".$Regions[$nextReg[$cnvOK]]{"End"}."\t".$Results{$file}{$nextReg[0]}."\t".scalar@cleanCNV."\t".sprintf("%.3f",$cleanAverage)."\t".($cnvOK+1)."\t".sprintf("%.3f",$dirtyAverage)."\t$overlapCNT\n";
						$Result3{$file}{$nextReg[0]} = $Results{$file}{$nextReg[0]};
						}
					}
				if ($i >= 1) { $r += $i; }
				else { $r++; }
				}
			}
		}
	close CNV1;
	close CNV2;
	}


##print graph foreach Chrom
if ($CNVgraph) {
	for my$file (@Files) {
		graphChr1($nGraf,"$outdir/$sampleName{$file}",$norm,$seuil_deletion,$seuil_duplication,$file,$maxDepthGraph,\@Files,\%sampleName,\@Regions,\@ChromOrder,\%regionOrder,\%Patients,\%Result2);
		}
	}


##summary
open (OUT,">$outdir/CNV.summary.txt") or die "Pb lors de l'ecriture du fichier sortie $!\n";
print OUT "CNV analysis\n";
print OUT "\nparameters:
	intervals normalization: $norm
	intervals detected as deletion if normalized depth below $seuil_deletion of $norm
	intervals detected as duplication if normalized depth above $seuil_duplication of $norm\n";
if ($seuil_patient <1) { print OUT "\tsamples kept if NO more than ".(100*$seuil_patient)."% of intervals are CNVs\n"; }
if ($seuil_region <1) { print OUT "\tintervals kept if NO more than ".(100*$seuil_region)."% of samples are CNVs\n"; }
if ($minCov) { print OUT "\tintervals discarded if at least one sample covered less than $minCov\n"; }
if ($minDP) { print OUT "\tintervals discarded if mean/median depth less than $minDP\n"; }
if ($minCNV) { print OUT "\tCNV calling requires at least $minCNV consecutive called intervals\n"; }
if ($maxNonCNV) { print OUT "\tCNV calling tolerates NO more than $maxNonCNV non-CNV inside\n"; }
if ($normByGender eq "all") { print OUT "\tinterval depth normalization by gender for all chromosomes\n"; }
elsif ($normByGender eq "gono") { print OUT "\tinterval depth normalization by gender for sex chromosomes only\n"; }
if ($Ref eq "mean") { print OUT "\teach sample normalized by the mean of all region average depths\n"; }
else { print OUT "\teach sample normalized by the sum of sequenced bp\n"; }
if ($RefByGender) { print OUT "\t(during sample normalization, depths within chrX are doubled for males, and those within chrY are skipped)\n"; }
else { print OUT "\t(during sample normalization, depths from all chr are taken, whatever the sex)\n"; }

if ($depthTxt) { print OUT "bases counts:\n$depthTxt\n"; }
if ($sexTxt) { print OUT "\ngender determination:\n$sexTxt\n"; }

print OUT "\nintervals nber: ".scalar@Regions."\n";
print OUT "\nintervals discarded:\n";
my$N_noData=0; my$N_CNV_Recurrent=0; my$N_lowCov=0;
foreach my$r (0..$#Regions) {
	if(defined $Regions[$r]{"Appel"}) {
		if ($Regions[$r]{"Appel"} eq "No_Data") { $N_noData++; }
		elsif ($Regions[$r]{"Appel"} eq "CNV_Recurrent") { $N_CNV_Recurrent++; }
		elsif ($Regions[$r]{"Appel"} eq "Couverture_Faible")  { $N_lowCov++; }
		}
	}
print OUT "\tnot enough data: $N_noData\n";
print OUT "\tCNV_Recurrent: $N_CNV_Recurrent\n";
if ($minDP || $minCov) { print OUT "\tlow coverage :  $N_lowCov\n"; }

print OUT "\nsamples discarded:\n";
my$N_discard=0;
foreach my$file (@Files) {
		if ($Patients{$file}{"ecarte"}) { print OUT "\t$sampleName{$file}\n"; $N_discard++; }
		}
unless ($N_discard) { print OUT "\tnone\n"; }

print OUT "\nResults
patient\tCNV\tnber\n";
foreach my$file (@Files) {
	unless($Patients{$file}{"ecarte"}) {
		my$N_dup=0; my$N_del=0;
		foreach (keys%{ $Result3{$file} }) {
			if ($Result3{$file}{$_} eq "DUP") { $N_dup++; }
			elsif ($Result3{$file}{$_} eq "DEL") { $N_del++; }
			}
		print OUT "\n$sampleName{$file}: 
		DUP: $N_dup
		DEL: $N_del\n";
		}
	}
close OUT;


return(\%Patients,\%Result4);

}



####################
####################
####################


#CNV_tool::CNV_reAnalyse($fichier_cov,$outdir,$nGraf,$fichier_sexe,\%CNV_opt);
sub CNV_reAnalyse
{
my($fichier_cov,$outdir,$nGraf,$fichier_sexe,$h1)=@_;
my%CNV_opt = %$h1;
my$norm = $CNV_opt{"norm"};
my$normByGender = $CNV_opt{"normByGender"};
my$Ref = $CNV_opt{"RefDepth"};
my$RefByGender = $CNV_opt{"RefByGender"};
my$seuil_region = $CNV_opt{"seuil_region"};
my$seuil_patient = $CNV_opt{"seuil_patient"};
my$seuil_cov = $CNV_opt{"seuil_cov"};
my$seuil_deletion = $CNV_opt{"seuil_deletion"};
my$seuil_duplication = $CNV_opt{"seuil_duplication"};
my$minCNV = $CNV_opt{"min_following_CNV"};
my$minDP = $CNV_opt{"min_DP"} ;
my$maxNonCNV = $CNV_opt{"max_Non_CNV"};
my$maxNonCNVrate = $CNV_opt{"max_Non_CNV_rate"};
my$graphByChr = $CNV_opt{"chromGraph"};
my@cnvFields = @{ $CNV_opt{"fields"} };
my$ok = 0;
foreach my$i (0..$#cnvFields) {
	if ($cnvFields[$i] eq $norm) { 
		$ok = 1;
		if ($norm eq "std") { @cnvFields = (@cnvFields[0..($i-1)],"moy",@cnvFields[$i..$#cnvFields]); }
		last;
		}
	}
my@cnvVal = @cnvFields;		# = @cnvFields plus $norm (if not already present)
unless ($ok) { 
	if ($norm eq "std") { push(@cnvVal, ("moy","std")); }
	else { push(@cnvVal, $norm); }
	}

my$maxDepthGraph = 10;

my %Sexe;
if ($fichier_sexe) {
	open(PATIENTS, "<$fichier_sexe") or die "Fichier $fichier_sexe d'attribution des sexes inexistant ou impossible a lire\n";
	foreach my $ligne (<PATIENTS>) {
		$ligne =~ m/(\S+)\s+(\S+)/;
		$Sexe{$1} = $2;
	}
	close(PATIENTS);
	foreach my $patient (keys%Sexe) {
		if ($Sexe{$patient} =~ /^m|^h/i) { $Sexe{$patient} = "H"; }
		else { $Sexe{$patient} = "F"; }
	}
}


# Pour chaque region de la sortie DecoVa, on recupere les infos de couverture
my @idxP;	#idx of Patients
my %idxV;	#idx of Values foreach patient
my $patient_nbr;
my %Patients;
my $region_nbr = 0;
my $autosom_Regions = 0;
my $gonosom_Regions = 0;
my @Regions;
my @ChrOrder; my%allChr;
my $sexe_courant;

# PREMIER PARCOURS #
# PERMET DE CALCULER LA PROFONDEUR TOTALE POUR CHAQUE PATIENT, EN TENANT UNIQUEMENT COMPTE DES REGIONS CORRECTEMENT SEQUENCEES
my$l=0; #N° ligne
unless (-f "$fichier_cov") { die "Fichier $fichier_cov inexistant\n"; }
open(DECOVA, "<$fichier_cov") or die "Fichier $fichier_cov impossible a lire\n";
while (my $line = <DECOVA>) {
	$l++;
	chomp($line);
	my @INFOS = split("\t",$line);

	if ($l==1) {
		my$p=0;	#N° patient
		for (my$i=1;$i<scalar@INFOS;$i++) {
			if (($INFOS[$i-1]eq"")&&($INFOS[$i]ne"")) { 
				push(@idxP, $i);
				$Patients{$p}{"ID"} = $INFOS[$i];
				if ($fichier_sexe) {
					if (exists $Sexe{$INFOS[$i]})
						{ $Patients{$p}{"Sexe"} = $Sexe{$INFOS[$i]}; }
					else { print "patient $INFOS[$i] not found in $fichier_sexe\n"; }
					}
				$Patients{$p}{"tot_Autosomes_depth"} = 0;
				$Patients{$p}{"tot_chrX_depth"} = 0;
				$Patients{$p}{"tot_chrY_depth"} = 0;
				$Patients{$p}{"Ref_depth"} = 0;
				$p++;
			}
		$patient_nbr = scalar@idxP;
		}
	} elsif ($l==2) {
		my$p=0;
		for (my$i=$idxP[0];$i<scalar@INFOS;$i++) {
			if ($p < $#idxP) { 
				if ($i == $idxP[$p+1]) { $p++; }
			}
			if ($INFOS[$i] eq "mean") { $idxV{$p}{"mean"} = $i; } 
			elsif ($INFOS[$i] =~ /%>=.+X/) { $idxV{$p}{"cov"} = $i; }
			elsif ($INFOS[$i] eq "tot") { $idxV{$p}{"tot"} = $i; }
			if (($p > 0) && (!exists$idxV{($p-1)}{"tot"}) && ($Ref eq "tot")) { die "tot column not found in DeCovA output\n"; }
		}
	} elsif ( ($line =~ m/^\w+\t\d+\t\d+/) && ($line !~ m/^#/) ) {

		# On recupere les informations generales de la region
		$Regions[$region_nbr]{"Chrom"} = $INFOS[0];
		$Regions[$region_nbr]{"start"} = $INFOS[1];
		$Regions[$region_nbr]{"end"} = $INFOS[2];
		if ($idxP[0]>3) { $Regions[$region_nbr]{"Gene"} = $INFOS[3]; }
		else { $Regions[$region_nbr]{"Gene"} = "NA"; }
		for (my$i=0;$i<$idxP[0];$i++)
			{ $Regions[$region_nbr]{"allInfos"} .= $INFOS[$i]."\t"; }
		chop $Regions[$region_nbr]{"allInfos"};
		# to get chrom order from bed
		unless (exists $allChr{$INFOS[0]}) {
			push(@ChrOrder, $INFOS[0]);
			$allChr{$INFOS[0]} = 1;
		}

		# On parcourt la fin de la region (= infos sur les patients uniquement) une premiere fois pour savoir si la region est a conserver
		my $region_a_conserver = 1;
		if ($seuil_cov) {
			for (my$p=0;$p<scalar@idxP;$p++) {
				# Lecture de chaque colonne "%>=.X" : si aucun patient n'a 100% de bases séquencées à .X pour la région, elle est supprimée de l'analyse
				my$cov = $INFOS[$idxV{$p}{"cov"}];
				$cov =~ s/%$//; $cov /= 100;
				if ($cov > $seuil_cov) {
					$region_a_conserver = 0; last;
				}
			}
		}
		if ($minDP) {
			my@allDepth = ();
			for (my$p=0;$p<scalar@idxP;$p++) {
				push(@allDepth, $INFOS[$idxV{$p}{"mean"}]);
			}
			my$normDepth = norm($norm,\@allDepth);
			if ( $normDepth < $minDP)
				{ $region_a_conserver = 0; }
		}

		# Si la region est mal couverte, on l'imprime dans un fichier POUBELLE
		if ($region_a_conserver != 1) {
			$Regions[$region_nbr]{"Appel"} = "Couverture_Faible";

		} else {
		# Sinon, on parcourt la region une seconde fois pour recuperer les informations désirées pour chacun des patients (et donc la profondeur moyenne pour la region)
			for (my $patient = 0 ; $patient < $patient_nbr ; $patient++) {
				$Regions[$region_nbr]{$patient}{"raw_depth"} = $INFOS[$idxV{$patient}{"mean"}];
				if ($Ref eq "mean") {
					$Regions[$region_nbr]{$patient}{"for_Tot_depth"} = $INFOS[$idxV{$patient}{"mean"}];
				} else {
					$Regions[$region_nbr]{$patient}{"for_Tot_depth"} = $INFOS[$idxV{$patient}{"tot"}];
				}
				if ($INFOS[0] !~ m/^(chr[XY]|[XY])$/) {
					$autosom_Regions++;
					$Patients{$patient}{"tot_Autosomes_depth"} += $Regions[$region_nbr]{$patient}{"for_Tot_depth"};
				} else {
					if ($INFOS[0] =~ m/^(chrX|X)$/) {
						$gonosom_Regions++;
						$Patients{$patient}{"tot_chrX_depth"} += $Regions[$region_nbr]{$patient}{"for_Tot_depth"};
					} else {
						$Patients{$patient}{"tot_chrY_depth"} += $Regions[$region_nbr]{$patient}{"for_Tot_depth"};

					} 
				}
			}
		}
		$region_nbr++;
	}
}
close(DECOVA);

if ($seuil_region < (1/$patient_nbr)) { $seuil_region = (1/$patient_nbr); }
if ($seuil_patient < (1/$region_nbr)) { $seuil_patient = (1/$region_nbr); }

my$sexTxt="";
for (my $patient = 0 ; $patient < $patient_nbr ; $patient++) {
	$sexTxt .= $Patients{$patient}{"ID"}.":\n";
	if (! defined $Patients{$patient}{"Sexe"}) {
		if ($Patients{$patient}{"tot_chrX_depth"} && $gonosom_Regions && $Patients{$patient}{"tot_Autosomes_depth"} && $autosom_Regions) {
			$sexTxt .= "\tautoZ: ".($Patients{$patient}{"tot_Autosomes_depth"}/$autosom_Regions)."\n\tchrX: ".($Patients{$patient}{"tot_chrX_depth"}/$gonosom_Regions)."\n";
			if ( ($Patients{$patient}{"tot_chrX_depth"}/$gonosom_Regions) > (1.2*$Patients{$patient}{"tot_Autosomes_depth"}/$autosom_Regions) ) { $sexTxt .= "\t-> sexe ambigu\n"; $Patients{$patient}{"Sexe"} = "F"; } 
			elsif ( (($Patients{$patient}{"tot_chrX_depth"}/$gonosom_Regions) <= (1.2*$Patients{$patient}{"tot_Autosomes_depth"}/$autosom_Regions)) && (($Patients{$patient}{"tot_chrX_depth"}/$gonosom_Regions) >= (0.8*$Patients{$patient}{"tot_Autosomes_depth"}/$autosom_Regions)) ) { $Patients{$patient}{"Sexe"} = "F"; } 
			elsif ( (($Patients{$patient}{"tot_chrX_depth"}/$gonosom_Regions) < (0.8*$Patients{$patient}{"tot_Autosomes_depth"}/$autosom_Regions)) && (($Patients{$patient}{"tot_chrX_depth"}/$gonosom_Regions) > (1.2*0.5*$Patients{$patient}{"tot_Autosomes_depth"}/$autosom_Regions)) ) { $sexTxt .= "\t-> sexe ambigu\n"; $Patients{$patient}{"Sexe"} = "F"; } 
			elsif ( (($Patients{$patient}{"tot_chrX_depth"}/$gonosom_Regions) <= (1.2*0.5*$Patients{$patient}{"tot_Autosomes_depth"}/$autosom_Regions)) && (($Patients{$patient}{"tot_chrX_depth"}/$gonosom_Regions) >= (0.8*0.5*$Patients{$patient}{"tot_Autosomes_depth"}/$autosom_Regions)) ) { $Patients{$patient}{"Sexe"} = "H"; }
			else { $sexTxt .= "\t-> sexe ambigu\n"; $Patients{$patient}{"Sexe"} = "H"; }
			
		} else { $Patients{$patient}{"Sexe"} = "F"; }
	}
	$sexTxt .= "\t-> $Patients{$patient}{Sexe}\n";

	if ($RefByGender) {
		if ($Patients{$patient}{"Sexe"} eq "H") {
			$Patients{$patient}{"Ref_depth"} = ($Patients{$patient}{"tot_chrX_depth"} * 2) + $Patients{$patient}{"tot_Autosomes_depth"}; 
		} else {
			$Patients{$patient}{"Ref_depth"} = $Patients{$patient}{"tot_chrX_depth"} + $Patients{$patient}{"tot_Autosomes_depth"};
		}
	} else {
		$Patients{$patient}{"Ref_depth"} = $Patients{$patient}{"tot_chrX_depth"} + $Patients{$patient}{"tot_chrY_depth"} + $Patients{$patient}{"tot_Autosomes_depth"};
	}
	unless ($Patients{$patient}{"Ref_depth"}) { $Patients{$patient}{"ecarte"} = 1; }
}
print "gender:\n$sexTxt\n";

my$meanRef;
for (my $patient = 0 ; $patient < $patient_nbr ; $patient++) {
	if ($Patients{$patient}{"Ref_depth"})
		{ $meanRef += $Patients{$patient}{"Ref_depth"}; }
	else { $Patients{$patient}{"ecarte"} = 1; }
	}
$meanRef /= $patient_nbr;
for (my $patient = 0 ; $patient < $patient_nbr ; $patient++) {
	$Patients{$patient}{"Ref_depth"} /= $meanRef; 
	print "normalization factor for $Patients{$patient}{ID} : ".$Patients{$patient}{"Ref_depth"}."\n";
	}



##iterations
my $nb_parcours = 0;
my $continuer = 1;
my %Results;
while ($continuer == 1 || $nb_parcours <= 1) {

	my $iteration = $nb_parcours + 1;
	print "Iteration numero \: ".$iteration."\n";

	undef %Results;
	my %CNV;
	my $nb_regions_conservees = 0;
	my $regions_ecartees = 0;
	$continuer = 0;
	my $sortie = "$outdir/CNV_Iteration_".$iteration."\.tab";
	open(SORTIE,">", $sortie) or die("Pb lors de l'ecriture du fichier sortie $!\n");
	my $log = "$outdir/Logs_Iteration_".$iteration."\.txt";
	open(LOG,">", $log) or die("Pb lors de l'ecriture du fichier log $!\n");

	print SORTIE "Chrom"."\t";
	print SORTIE "start"."\t";
	print SORTIE "end"."\t";
	print SORTIE "Gene"."\t";
	print SORTIE "Numero_Region"."\t";

	for (my $patient = 0 ; $patient < $patient_nbr ; $patient++) {
		printf SORTIE $Patients{$patient}{"ID"}."\t";
	}

	print SORTIE "Statut_Region\n";


	# PREMIER PARCOURS BIS #
	# RECALCUL DE LA PROFONDEUR TOTALE POUR LA PONDERATION INTRA
	# CE (RE)CALCUL A LIEU SI DES REGIONS ETAIENT "MOCHES" APRES L'APPEL DE CNV
	if ($nb_parcours > 0) {

		# REINITIALISATION DE LA PROFONDEUR TOTALE PAR PATIENT
		for (my $patient = 0 ; $patient < $patient_nbr ; $patient++) {
			$Patients{$patient}{"tot_Autosomes_depth"} = 0;
			$Patients{$patient}{"tot_sexChr_depth"} = 0;
		}

		for (my $r = 0 ; $r < $region_nbr ; $r++) {
			if(!(defined($Regions[$r]{"Appel"}))) {
				for (my $patient = 0 ; $patient < $patient_nbr ; $patient++) {
					if ($Regions[$r]{"Chrom"} !~ m/^chr[XY]$|[XY]$/) {
						$Patients{$patient}{"tot_Autosomes_depth"} += $Regions[$r]{$patient}{"for_Tot_depth"};
					} elsif ($Regions[$r]{"Chrom"} =~ m/^chrX$|^X$/) {
						$Patients{$patient}{"tot_sexChr_depth"} += $Regions[$r]{$patient}{"for_Tot_depth"};
					} else {
						if ((!$RefByGender) && ($Patients{$patient}{"Sexe"} eq "H")) {
							$Patients{$patient}{"tot_sexChr_depth"} += $Regions[$region_nbr]{$patient}{"for_Tot_depth"};
						}
					}
				}
			}
		}


		for (my $p = 0 ; $p < $patient_nbr ; $p++) {
			if ($RefByGender && ($Patients{$p}{"Sexe"} eq "H")) {
				$Patients{$p}{"Ref_depth"} = ($Patients{$p}{"tot_sexChr_depth"} * 2) + $Patients{$p}{"tot_Autosomes_depth"};
			} else {
				$Patients{$p}{"Ref_depth"} = $Patients{$p}{"tot_sexChr_depth"} + $Patients{$p}{"tot_Autosomes_depth"};
			}
		}

		for (my $p = 0 ; $p < $patient_nbr ; $p++)
			{ $meanRef += $Patients{$p}{"Ref_depth"}; }
		$meanRef /= $patient_nbr;
		for (my $p = 0 ; $p < $patient_nbr ; $p++)
			{ $Patients{$p}{"Ref_depth"} /= $meanRef; }

	}
	print LOG "sample normalization factors :\n";
	for (my $p = 0 ; $p < $patient_nbr ; $p++) {
		print LOG "\t$Patients{$p}{ID} : $Patients{$p}{Ref_depth}\n";
	}
	print LOG "\n";


	# SECOND PARCOURS DES REGIONS #
	# PERMET DE PONDERER LA PROFONDEUR PAR LES AUTRES REGIONS (INTRA) ET ENTRE LES PATIENTS (INTER)
	for (my $r = 0 ; $r < $region_nbr ; $r++) {

		print SORTIE $Regions[$r]{"Chrom"}."\t";
		print SORTIE $Regions[$r]{"start"}."\t";
		print SORTIE $Regions[$r]{"end"}."\t";
		print SORTIE $Regions[$r]{"Gene"}."\t";
		print SORTIE $r."\t";

		# SI LA REGION N'EST PAS "MOCHE"
		if (!defined($Regions[$r]{"Appel"})) {

			$nb_regions_conservees++;

			$Regions[$r]{"nb_CNV"}{"DEL"} = 0;
			$Regions[$r]{"nb_CNV"}{"DUP"} = 0;
			my $nb_evts = 0;
			my $prof_Moyenne_Inter = "";


			if ($normByGender) {

				# Pour chaque patient, un premier parcours permet la pondération intra-patient (en divisant la profondeur de la region par la profondeur totale du patient)
				@{ $Regions[$r]{"all_normByS_depths_fem"} } = ();
				@{ $Regions[$r]{"all_normByS_depths_males"} } = ();

				for (my $patient = 0 ; $patient < $patient_nbr ; $patient++) {

					if($Patients{$patient}{"Sexe"} eq "F" && !($Patients{$patient}{"ecarte"})) {

						# Le controle permet de diviser uniquement par le nombre de patients restant apres l filtre, et non par l'ensemble des patients
						# Si la region concerne un gonosome et le patient est une femme, la référence est celle qui correspond aux autosomes + au chromosome X
						$Regions[$r]{$patient}{"normByS_depth"} = $Regions[$r]{$patient}{"raw_depth"}/$Patients{$patient}{"Ref_depth"};
						# On sauvegarde la profondeur ponderee intra patient, pour determiner ensuite la moyenne de ces profondeurs ponderees
						push( @{ $Regions[$r]{"all_normByS_depths_fem"} }, $Regions[$r]{$patient}{"normByS_depth"} );


					} elsif ($Patients{$patient}{"Sexe"} eq "H" && !($Patients{$patient}{"ecarte"})) {

						# Le controle permet de divisier uniquement par le nombre de patients restant apres l filtre, et non par l'ensemble des patients
						# Si la region concerne un gonosome et le patient est un homme, la référence est celle qui correspond uniquement au chromosome X
						$Regions[$r]{$patient}{"normByS_depth"} = $Regions[$r]{$patient}{"raw_depth"}/$Patients{$patient}{"Ref_depth"};
						# On sauvegarde la profondeur ponderee intra patient, pour determiner ensuite la moyenne de ces profondeurs ponderees
						push( @{ $Regions[$r]{"all_normByS_depths_males"} }, $Regions[$r]{$patient}{"normByS_depth"} );

					}

				}

				# Nous calculons pour chaque region et chaque sexe la moyenne de profondeur (ponderee pour chaque patient par rapport aux autres regions)
				if ($normByGender eq "all") {
					if(@{ $Regions[$r]{"all_normByS_depths_fem"} }) {
						foreach (@cnvVal) { $Regions[$r]{"normByR_depth_fem"}{$_} = norm($_,\@{ $Regions[$r]{"all_normByS_depths_fem"} }); }
						if ($norm eq "std") {
							my$sqtotal = 0;
							foreach (@{ $Regions[$r]{"all_normByS_depths_fem"} })
								{ $sqtotal += ($Regions[$r]{"normByR_depth_fem"}{"moy"}-$_)**2; }
							$Regions[$r]{"normByR_depth_fem"}{"std"} = ($sqtotal / (scalar@{ $Regions[$r]{"all_normByS_depths_fem"} }-1))**0.5;
						}
					}
					if(@{ $Regions[$r]{"all_normByS_depths_males"} }) {
						foreach (@cnvVal) { $Regions[$r]{"normByR_depth_males"}{$_} = norm({$_},\@{ $Regions[$r]{"all_normByS_depths_males"} }); }
						if ($norm eq "std") {
							my$sqtotal = 0;
							foreach (@{ $Regions[$r]{"all_normByS_depths_males"} })
								{ $sqtotal += ($Regions[$r]{"normByR_depth_males"}{"moy"}-$_)**2; }
							$Regions[$r]{"normByR_depth_males"}{"std"} = ($sqtotal / (scalar@{ $Regions[$r]{"all_normByS_depths_males"} }-1))**0.5;
						}
					}
				} else {
					if ($Regions[$r]{"Chrom"} !~ m/^chr[XY]$|^[XY]$/) {
						my@Autosomes=();
						if(@{ $Regions[$r]{"all_normByS_depths_fem"} }) {
							push(@Autosomes, @{ $Regions[$r]{"all_normByS_depths_fem"} });
						}
						if(@{ $Regions[$r]{"all_normByS_depths_males"} }) {
							push(@Autosomes, @{ $Regions[$r]{"all_normByS_depths_males"} });
						}
						if (@Autosomes) {
							foreach (@cnvVal) { $Regions[$r]{"normByR_depth"}{$_} = norm($_,\@Autosomes); }
							if ($norm eq "std") {
								my$sqtotal = 0;
								foreach (@Autosomes)
									{ $sqtotal += ($Regions[$r]{"normByR_depth"}{"moy"}-$_)**2; }
								$Regions[$r]{"normByR_depth"}{"std"} = ($sqtotal / (scalar@Autosomes-1))**0.5;
							}
						}
					} else {
						if(@{ $Regions[$r]{"all_normByS_depths_fem"} }) {
							foreach (@cnvVal) { $Regions[$r]{"normByR_depth_fem"}{$_} = norm($_,\@{ $Regions[$r]{"all_normByS_depths_fem"} }); }
							if ($norm eq "std") {
								my$sqtotal = 0;
								foreach (@{ $Regions[$r]{"all_normByS_depths_fem"} })
									{ $sqtotal += ($Regions[$r]{"normByR_depth_fem"}{"moy"}-$_)**2; }
								$Regions[$r]{"normByR_depth_fem"}{"std"} = ($sqtotal / (scalar@{ $Regions[$r]{"all_normByS_depths_fem"} }-1))**0.5;
							}
						}
						if(@{ $Regions[$r]{"all_normByS_depths_males"} }) {
							foreach (@cnvVal) { $Regions[$r]{"normByR_depth_males"}{$_} = norm($_,\@{ $Regions[$r]{"all_normByS_depths_males"} }); }
							if ($norm eq "std") {
								my$sqtotal = 0;
								foreach (@{ $Regions[$r]{"all_normByS_depths_males"} })
									{ $sqtotal += ($Regions[$r]{"normByR_depth_males"}{"moy"}-$_)**2; }
								$Regions[$r]{"normByR_depth_males"}{"std"} = ($sqtotal / (scalar@{ $Regions[$r]{"all_normByS_depths_males"} }-1))**0.5;
							}
						}
					}
				}
				
				# Un second parcours permet l'appel de CNV (en comparant la profondeur pour un patient par la profondeur moyenne de la region)
				for (my $patient = 0 ; $patient < $patient_nbr ; $patient++) {

					if ($norm eq "std") { $prof_Moyenne_Inter = 0; }
					else  { $prof_Moyenne_Inter = 1; }

					unless ($Patients{$patient}{"ecarte"}) {

						if ($normByGender eq "all") {

							if($Patients{$patient}{"Sexe"} eq "F") {
								if ($Regions[$r]{"normByR_depth_fem"}{$norm}) {
									if ($Regions[$r]{"Chrom"} !~ m/^(chrY|Y)$/) {
										if ($norm eq "std") {
											$Regions[$r]{$patient}{"depth_ratio"}{"moy"} = $Regions[$r]{$patient}{"normByS_depth"}/$Regions[$r]{"normByR_depth_fem"}{"moy"};
											$Regions[$r]{$patient}{"depth_ratio"}{"std"} = ($Regions[$r]{$patient}{"normByS_depth"}-$Regions[$r]{"normByR_depth_fem"}{"moy"})/$Regions[$r]{"normByR_depth_fem"}{"std"};
										} else {
											$Regions[$r]{$patient}{"depth_ratio"}{$norm} = $Regions[$r]{$patient}{"normByS_depth"}/$Regions[$r]{"normByR_depth_fem"}{$norm};
										}
										$prof_Moyenne_Inter = $Regions[$r]{$patient}{"depth_ratio"}{$norm};
										print SORTIE "$prof_Moyenne_Inter\t";
									} else {
										# En faisant la distinction entre chromosomes X et Y, les femmes ne sont plus considérées pour l'appel de CNV des régions du dernier
										print SORTIE "NA\t";
									}
								} else {
									$Regions[$r]{"Appel"} = "No_Data";
									print SORTIE "NA\t";
								}

							} elsif ($Patients{$patient}{"Sexe"} eq "H") {
								if ($Regions[$r]{"normByR_depth_males"}{$norm}) {
									if ($norm eq "std") {
										$Regions[$r]{$patient}{"depth_ratio"}{"moy"} = $Regions[$r]{$patient}{"normByS_depth"}/$Regions[$r]{"normByR_depth_males"}{"moy"};
										$Regions[$r]{$patient}{"depth_ratio"}{"std"} = ($Regions[$r]{$patient}{"normByS_depth"}-$Regions[$r]{"normByR_depth_males"}{"moy"})/$Regions[$r]{"normByR_depth_males"}{"std"};
									} else {
										$Regions[$r]{$patient}{"depth_ratio"}{$norm} = $Regions[$r]{$patient}{"normByS_depth"}/$Regions[$r]{"normByR_depth_males"}{$norm};
									}
									$prof_Moyenne_Inter = $Regions[$r]{$patient}{"depth_ratio"}{$norm};
									print SORTIE "$prof_Moyenne_Inter\t";
								} else {
									$Regions[$r]{"Appel"} = "No_Data";
									print SORTIE "NA\t";
								}
							}
						} else {
							# Pour les autosomes
							if ($Regions[$r]{"Chrom"} !~ m/^chr[XY]$|^[XY]$/) {
								if ($Regions[$r]{"normByR_depth"}{$norm}) {
									if ($norm eq "std") {
										$Regions[$r]{$patient}{"depth_ratio"}{"moy"} = $Regions[$r]{$patient}{"normByS_depth"} / $Regions[$r]{"normByR_depth"}{"moy"};
										$Regions[$r]{$patient}{"depth_ratio"}{"std"} = ($Regions[$r]{$patient}{"normByS_depth"} - $Regions[$r]{"normByR_depth"}{"moy"}) / $Regions[$r]{"normByR_depth"}{"std"};
									} else {
										$Regions[$r]{$patient}{"depth_ratio"}{$norm} = $Regions[$r]{$patient}{"normByS_depth"} / $Regions[$r]{"normByR_depth"}{$norm};
									}
									$prof_Moyenne_Inter = $Regions[$r]{$patient}{"depth_ratio"}{$norm};
									print SORTIE "$prof_Moyenne_Inter\t";
								} else {
									$Regions[$r]{"Appel"} = "No_Data";
									print SORTIE "NA\t";
								}
							} else {
								if($Patients{$patient}{"Sexe"} eq "F") {
									if ($Regions[$r]{"normByR_depth_fem"}{$norm}) {	
										if ($Regions[$r]{"Chrom"} =~ m/^(chrX|X)$/) {
											if ($norm eq "std") {
												$Regions[$r]{$patient}{"depth_ratio"}{"moy"} = $Regions[$r]{$patient}{"normByS_depth"}/$Regions[$r]{"normByR_depth_fem"}{"moy"};
												$Regions[$r]{$patient}{"depth_ratio"}{"std"} = ($Regions[$r]{$patient}{"normByS_depth"}-$Regions[$r]{"normByR_depth_fem"}{"moy"})/$Regions[$r]{"normByR_depth_fem"}{"std"};
											} else {
												$Regions[$r]{$patient}{"depth_ratio"}{$norm} = $Regions[$r]{$patient}{"normByS_depth"}/$Regions[$r]{"normByR_depth_fem"}{$norm};
											}
											$prof_Moyenne_Inter = $Regions[$r]{$patient}{"depth_ratio"}{$norm};
											print SORTIE "$prof_Moyenne_Inter\t";
										} else {
											# En faisant la distinction entre chromosomes X et Y, les femmes ne sont plus considérées pour l'appel de CNV des régions du dernier
											print SORTIE "NA\t";
										}
									} else {
										$Regions[$r]{"Appel"} = "No_Data";
										print SORTIE "NA\t";
									}
								} elsif ($Patients{$patient}{"Sexe"} eq "H") {
									if ($Regions[$r]{"normByR_depth_males"}{$norm}) {
										if ($norm eq "std") {
											$Regions[$r]{$patient}{"depth_ratio"}{"moy"} = $Regions[$r]{$patient}{"normByS_depth"}/$Regions[$r]{"normByR_depth_males"}{"moy"};
											$Regions[$r]{$patient}{"depth_ratio"}{"std"} = ($Regions[$r]{$patient}{"normByS_depth"}-$Regions[$r]{"normByR_depth_males"}{"moy"})/$Regions[$r]{"normByR_depth_males"}{"std"};
										} else {
											$Regions[$r]{$patient}{"depth_ratio"}{$norm} = $Regions[$r]{$patient}{"normByS_depth"}/$Regions[$r]{"normByR_depth_males"}{$norm};
										}
										$prof_Moyenne_Inter = $Regions[$r]{$patient}{"depth_ratio"}{$norm};
										print SORTIE "$prof_Moyenne_Inter\t";
									} else {
										$Regions[$r]{"Appel"} = "No_Data";
										print SORTIE "NA\t";
									}
								}
							}
						}
				
						if ($prof_Moyenne_Inter < $seuil_deletion) {
							$Regions[$r]{"nb_CNV"}{"DEL"}++;
							$nb_evts++;
							$CNV{$patient}++;
							$Results{$patient}{$r} = "DEL";
						} elsif ($prof_Moyenne_Inter > $seuil_duplication) {
							$Regions[$r]{"nb_CNV"}{"DUP"}++;
							$nb_evts++;
							$CNV{$patient}++;
							$Results{$patient}{$r} = "DUP";
						}

					} else {
						print SORTIE "NA\t";
					}

				}

				my$recurrent="";
				if(@{ $Regions[$r]{"all_normByS_depths_fem"} } && @{ $Regions[$r]{"all_normByS_depths_males"} } && (($nb_evts/(scalar@{ $Regions[$r]{"all_normByS_depths_fem"} } + scalar@{ $Regions[$r]{"all_normByS_depths_males"} })) > $seuil_region) && ($nb_parcours > 0)) { $recurrent = 1; }
				elsif(@{ $Regions[$r]{"all_normByS_depths_fem"} } && !@{ $Regions[$r]{"all_normByS_depths_males"} } && (($nb_evts/scalar@{ $Regions[$r]{"all_normByS_depths_fem"} }) > $seuil_region) && ($nb_parcours > 0)) { $recurrent = 1; }
				elsif(!@{ $Regions[$r]{"all_normByS_depths_fem"} } && @{ $Regions[$r]{"all_normByS_depths_males"} } && (($nb_evts/scalar@{ $Regions[$r]{"all_normByS_depths_males"} }) > $seuil_region) && ($nb_parcours > 0)) { $recurrent = 1; }
				if ($recurrent) {
					print SORTIE "CNV_Recurrent";
					print LOG "Region ecartee \: ". $Regions[$r]{"Chrom"}."\t".$Regions[$r]{"start"}."\t".$Regions[$r]{"end"}."\t".$Regions[$r]{"Gene"}."\n";
					$Regions[$r]{"Appel"} = "CNV_Recurrent";
					$regions_ecartees++;
					$continuer = 1;
				} else {		
					print SORTIE "OK";
				}

				print SORTIE "\n";



			} else {

				@{ $Regions[$r]{"all_normByS_depths"} }=();
				# Pour chaque patient, un premier parcours permet la pondération intra-patient (en divisant la profondeur de la region par la profondeur totale recalculée du patient)
				for (my $patient = 0 ; $patient < $patient_nbr ; $patient++) {
					unless ($Patients{$patient}{"ecarte"}) {
						$Regions[$r]{$patient}{"normByS_depth"} = $Regions[$r]{$patient}{"raw_depth"}/$Patients{$patient}{"Ref_depth"};
						if ( ($Patients{$patient}{"Sexe"} eq "H") && ($Regions[$r]{"Chrom"} =~ m/^(chrX|X)$/) )
							{ push(@{ $Regions[$r]{"all_normByS_depths"} }, ($Regions[$r]{$patient}{"normByS_depth"}*2)); }
						elsif ( ($Patients{$patient}{"Sexe"} eq "F") && ($Regions[$r]{"Chrom"} =~ m/^(chrY|Y)$/) )
							{ next ; }
						else { push(@{ $Regions[$r]{"all_normByS_depths"} }, $Regions[$r]{$patient}{"normByS_depth"}); }
					}
				}
				# Nous calculons pour chaque region la mediane/moyenne de profondeur (ponderee pour chaque patient par rapport aux autres regions)
				if(@{ $Regions[$r]{"all_normByS_depths"} }) {
					foreach (@cnvVal) {
						$Regions[$r]{"normByR_depth"}{$_} = norm($_,\@{ $Regions[$r]{"all_normByS_depths"} });
						}
					if ($norm eq "std") {
						my$sqtotal = 0;
						foreach (@{ $Regions[$r]{"all_normByS_depths"} })
							{ $sqtotal += ($Regions[$r]{"normByR_depth"}{"moy"}-$_)**2; }
						$Regions[$r]{"normByR_depth"}{"std"} = ($sqtotal / (scalar@{ $Regions[$r]{"all_normByS_depths"} }-1))**0.5;
					}
				}

				# Un second parcours permet l'appel de CNV (en comparant la profondeur pour un patient par la profondeur moyenne de la region)
				for (my $patient = 0 ; $patient < $patient_nbr ; $patient++) {

					if ($norm eq "std") { $prof_Moyenne_Inter = 0; }
					else  { $prof_Moyenne_Inter = 1; }

					unless ($Patients{$patient}{"ecarte"}) {

						if ($Regions[$r]{"normByR_depth"}{$norm}) {
							if($Patients{$patient}{"Sexe"} eq "F") {
								if ($Regions[$r]{"Chrom"} !~ m/^(chrY|Y)$/) {
									if ($norm eq "std") {
										$Regions[$r]{$patient}{"depth_ratio"}{"moy"} = $Regions[$r]{$patient}{"normByS_depth"}/$Regions[$r]{"normByR_depth"}{"moy"};
										$Regions[$r]{$patient}{"depth_ratio"}{"std"} = ($Regions[$r]{$patient}{"normByS_depth"} - $Regions[$r]{"normByR_depth"}{"moy"}) / $Regions[$r]{"normByR_depth"}{"std"};
									} else {
										$Regions[$r]{$patient}{"depth_ratio"}{$norm} = $Regions[$r]{$patient}{"normByS_depth"}/$Regions[$r]{"normByR_depth"}{$norm};
									}
									$prof_Moyenne_Inter = $Regions[$r]{$patient}{"depth_ratio"}{$norm};
									print SORTIE "$prof_Moyenne_Inter\t";
								} else {
									#les femmes ne sont pas considérées pour l'appel de CNV des régions du Y
									print SORTIE "NA\t";
								}
							} else {
								#*2 for chrX
								if ($Regions[$r]{"Chrom"} =~ m/^(chrX|X)$/) {
									if ($norm eq "std") {
										$Regions[$r]{$patient}{"depth_ratio"}{"moy"} = $Regions[$r]{$patient}{"normByS_depth"}*2 / $Regions[$r]{"normByR_depth"}{"moy"};
										$Regions[$r]{$patient}{"depth_ratio"}{"std"} = ($Regions[$r]{$patient}{"normByS_depth"}*2 - $Regions[$r]{"normByR_depth"}{"moy"}) / $Regions[$r]{"normByR_depth"}{"std"};
									} else {
										$Regions[$r]{$patient}{"depth_ratio"}{$norm} = $Regions[$r]{$patient}{"normByS_depth"}*2 / $Regions[$r]{"normByR_depth"}{$norm};
									}
								} else {
									if ($norm eq "std") {
										$Regions[$r]{$patient}{"depth_ratio"}{"moy"} = $Regions[$r]{$patient}{"normByS_depth"} / $Regions[$r]{"normByR_depth"}{"moy"};
										$Regions[$r]{$patient}{"depth_ratio"}{"std"} = ($Regions[$r]{$patient}{"normByS_depth"} - $Regions[$r]{"normByR_depth"}{"moy"}) / $Regions[$r]{"normByR_depth"}{"std"};
									} else {
										$Regions[$r]{$patient}{"depth_ratio"}{$norm} = $Regions[$r]{$patient}{"normByS_depth"} / $Regions[$r]{"normByR_depth"}{$norm};
									}
								}
								$prof_Moyenne_Inter = $Regions[$r]{$patient}{"depth_ratio"}{$norm};
								print SORTIE "$prof_Moyenne_Inter\t";
							}
							if ($prof_Moyenne_Inter < $seuil_deletion) {
								$Regions[$r]{"nb_CNV"}{"DEL"}++;
								$nb_evts++;
								$CNV{$patient}++;
								$Results{$patient}{$r} = "DEL";
							} elsif ($prof_Moyenne_Inter > $seuil_duplication) {
								$Regions[$r]{"nb_CNV"}{"DUP"}++;
								$nb_evts++;
								$CNV{$patient}++;
								$Results{$patient}{$r} = "DUP";
							}
						} else {
							$Regions[$r]{"Appel"} = "No_Data";
							print SORTIE "NA\t";
						}
					} else {
						print SORTIE "NA\t";
					}
				}

				if(@{ $Regions[$r]{"all_normByS_depths"} } && (($nb_evts/scalar@{ $Regions[$r]{"all_normByS_depths"} }) > $seuil_region) && ($nb_parcours > 0)) {
					print SORTIE "CNV_Recurrent";
					print LOG "Region ecartee \: ". $Regions[$r]{"Chrom"}."\t".$Regions[$r]{"start"}."\t".$Regions[$r]{"end"}."\t".$Regions[$r]{"Gene"}."\n";
					$Regions[$r]{"Appel"} = "CNV_Recurrent";
					$regions_ecartees++;
					$continuer = 1;

				} else {		
					print SORTIE "OK";
				}
				print SORTIE "\n";
			}


		# SI LA REGION EST "MOCHE"
		} else {

			for (my $patient = 0 ; $patient < $patient_nbr ; $patient++) {
				print SORTIE "NA\t";
			}
			print SORTIE $Regions[$r]{"Appel"}."\n";

		}

	}

	print LOG "Nombre de regions ecartees lors de cette iteration \: ".$regions_ecartees."\n\n\n";	

	my $patients_ecartes = 0;

	for (my $patient = 0 ; $patient < $patient_nbr ; $patient++) {

		unless ($Patients{$patient}{"ecarte"}) {
			my $prcent_evt;
			if (defined($CNV{$patient})) {

				$prcent_evt = $CNV{$patient}/$nb_regions_conservees;

				if ($prcent_evt >= $seuil_patient) {
					print LOG "Patient ecarte \: ".$Patients{$patient}{"ID"}."\n";
					$Patients{$patient}{"ecarte"} = 1;
					$continuer = 1;
					$patients_ecartes++;			
				} else {
					$Patients{$patient}{"ecarte"} = 0;
				}


			} else {
				print $Patients{$patient}{"ID"}." \: aucun CNV identifie\n";
				$Patients{$patient}{"ecarte"} = 0;
			}


		}

	}

	print LOG "Nombre de patients ecartes lors de cette iteration \: ".$patients_ecartes."\n";

	#print "\n";
	close(SORTIE);
	close(LOG);
	$nb_parcours++;
	
}


## @{ $regionOrder{$Chrom} } = [ regions sorted by pos ]
##put the first region with some coordinates in array (in case several regions with same coordinates exists)
my(%uniqReg,%RegInArray,%regionOrder,%regionIndice);
for my$r (0..$#Regions) {
	if (!exists $uniqReg{ $Regions[$r]{"Chrom"}."-".$Regions[$r]{"start"}."-".$Regions[$r]{"end"} }) {
		push (@{ $RegInArray{ $Regions[$r]{"Chrom"} } }, $r);
		$uniqReg{ $Regions[$r]{"Chrom"}."-".$Regions[$r]{"start"}."-".$Regions[$r]{"end"} } = $r;
		}
	if ($Regions[$r]{"Gene"} eq "NA")
		{ $Regions[$r]{"label"} = $Regions[$r]{"start"}."-".$Regions[$r]{"end"}; }
	else { 
		my@tab = split(/:|,/,$Regions[$r]{"Gene"});
		$Regions[$r]{"label"} = substr($tab[0], 0, 25);
		$Regions[$r]{"geneID"} = $tab[0];
		}
	}
foreach my$Chrom (keys%RegInArray) {
	##sort by ascending region ends (in case several regions with same start) then by ascending starts
	@{ $regionOrder{$Chrom} } = sort{$Regions[$a]{"end"}<=>$Regions[$b]{"end"}}@{ $RegInArray{$Chrom} }; 
	@{ $regionOrder{$Chrom} } = sort{$Regions[$a]{"start"}<=>$Regions[$b]{"start"}}@{ $RegInArray{$Chrom} };
	for (my$i=0;$i<scalar@{ $regionOrder{$Chrom} };$i++) { $regionIndice{ $regionOrder{$Chrom}[$i] } = $i; }
	}


##print CNV foreach sample
#$results{$patient}{$region} = "DUP";
#@{ $orderedCNV{$patient}{$Chrom} } = [r1, r2,...]
my%orderedCNV;
foreach my$patient (keys%Results) {
	my(%uniqCNV,%CNVinArray);
	foreach my$region (sort{$a<=>$b}keys%{ $Results{$patient} }) {
		if (!exists $uniqCNV{ $Regions[$region]{"Chrom"}."-".$Regions[$region]{"start"}."-".$Regions[$region]{"end"} }) { 
			push (@{ $CNVinArray{ $Regions[$region]{"Chrom"} } }, $region);
			$uniqCNV{ $Regions[$region]{"Chrom"}."-".$Regions[$region]{"start"}."-".$Regions[$region]{"end"} } = 1;
			}
		}
	foreach my$Chrom (keys%CNVinArray) {
		##sort by region ends (in case several regions with same start) then by starts
		@{ $orderedCNV{$patient}{$Chrom} } = sort{$Regions[$a]{"end"}<=>$Regions[$b]{"end"}}@{ $CNVinArray{$Chrom} }; 
		@{ $orderedCNV{$patient}{$Chrom} } = sort{$Regions[$a]{"start"}<=>$Regions[$b]{"start"}}@{ $CNVinArray{$Chrom} }; 
		}
	}

unless($minCNV) { $minCNV = 1; }
my%Result2;	##selected from %results
my%Result3;	##selected from %results, starting with first interval , ending with last interval
foreach my$patient (keys%Results) {
	print $Patients{$patient}{"ID"}."\n";
	open (CNV1,">$outdir/CNV_$Patients{$patient}{ID}.summary.txt") or die "Pb lors de l'ecriture du fichier sortie $!\n";
	print CNV1 "#Region\tCNV\tN_clean_intervals\tclean_ratio_to_$norm\tN_dirty_intervals\tdirty_ratio_to_$norm\toverlapping_Samples\n";
	open (CNV2,">$outdir/CNV_$Patients{$patient}{ID}.allIntervals.txt") or die "Pb lors de l'ecriture du fichier sortie $!\n";
	print CNV2 "#Chrom\tStart\tEnd\tLength\tInfo\tinterval_Order\tCNV_type\tratio_to_$norm\toccurences";
	if (@cnvFields) {
		foreach (@cnvFields) { 
			if ($_ eq "norm") {
				if ($norm eq "std") { print CNV2 "\tmoy\tstd"; }
				else  { print CNV2 "\t$norm"; }
				}
			else { print CNV2 "\t$_"; }
			}
		}
	print CNV2 "\n";
	foreach my$Chrom (@ChrOrder) {
		if (exists $orderedCNV{$patient}{$Chrom}) {
			my$r=0;		##index in @{ $orderedCNV{$patient}{$Chrom} }
			while ($r < scalar@{ $orderedCNV{$patient}{$Chrom} } ) {
				##merge consecutive CNVs
				my@refSub = mergeConsecutiveCNV("",$maxNonCNV,$orderedCNV{$patient}{$Chrom}[$r],$regionIndice{$orderedCNV{$patient}{$Chrom}[$r]},\@{ $regionOrder{$Chrom} },\%{ $Results{$patient} },\@Regions);
				my$i=$refSub[0]; my$cnvOK=$refSub[1]; my$nonCNVtot=$refSub[2];
				my@nextReg = @{ $refSub[3] };
				if ($maxNonCNVrate) {
					while (($nonCNVtot/($cnvOK+1)) > $maxNonCNVrate) {
						@refSub = mergeConsecutiveCNV(($cnvOK-1),$maxNonCNV,$orderedCNV{$patient}{$Chrom}[$r],$regionIndice{$orderedCNV{$patient}{$Chrom}[$r]},\@{ $regionOrder{$Chrom} },\%{ $Results{$patient} },\@Regions);
						$i=$refSub[0]; my$cnvOK=$refSub[1]; my$nonCNVtot=$refSub[2];
						@nextReg = @{ $refSub[3] };
						}
					}
				##overlapping samples
				my$overlapCNT=0; my%overlapSMPL=();
				for (my$j=0;$j<=$cnvOK;$j++) {
					foreach my$other (keys%Patients) {
						if ( (exists $Results{$other}{$nextReg[$j]}) && ($Results{$other}{$nextReg[$j]} eq $Results{$patient}{$nextReg[0]}) && (!exists $overlapSMPL{$other}) ) { 
							$overlapCNT++;
							$overlapSMPL{$other}=1;
							}
						}
					}
				##print patient's summary and allIntervals
				if ($cnvOK >= ($minCNV-1)) {
					if ($cnvOK == 0) {
						$Result2{$patient}{$nextReg[0]} = $Results{$patient}{$nextReg[0]};
						$Result3{$patient}{$nextReg[0]} = $Results{$patient}{$nextReg[0]};
						if (exists $regionIndice{$nextReg[0]}) {
							print CNV1 $Chrom.":".$Regions[$nextReg[0]]{"start"}."-".$Regions[$nextReg[0]]{"end"}."\t".$Results{$patient}{$nextReg[0]}."\t1\t".sprintf("%.3f",$Regions[$nextReg[0]]{$patient}{"depth_ratio"}{$norm})."\t.\t.\t$overlapCNT\n";
							print CNV2 $Chrom."\t".$Regions[$nextReg[0]]{"start"}."\t".$Regions[$nextReg[0]]{"end"}."\t".($Regions[$nextReg[0]]{"end"}-$Regions[$nextReg[0]]{"start"})." bp\t".$Regions[$nextReg[0]]{"label"}."\t".($regionIndice{$nextReg[0]}+1)."\t".$Results{$patient}{$nextReg[0]}."\t".sprintf("%.3f",$Regions[$nextReg[0]]{$patient}{"depth_ratio"}{$norm})."\t".$Regions[$nextReg[0]]{"nb_CNV"}{$Results{$patient}{$nextReg[0]}};
							if (@cnvFields) {
								my$txt = printCNVfields($normByGender,\@cnvFields,\%{ $Regions[$nextReg[0]] },\%{ $Patients{$patient} });
								print CNV2 "$txt";
								}
							print CNV2 "\n\n";
							}
						}
					else {
						my@cleanCNV = (); my@dirtyCNV = ();
						for (my$j=0;$j<=$cnvOK;$j++) {
							if (exists $Results{$patient}{$nextReg[$j]}) {
								$Result2{$patient}{$nextReg[$j]} = $Results{$patient}{$nextReg[$j]};
								push(@cleanCNV, $Regions[$nextReg[$j]]{$patient}{"depth_ratio"}{$norm});
								push(@dirtyCNV, $Regions[$nextReg[$j]]{$patient}{"depth_ratio"}{$norm});
								if (exists $regionIndice{$nextReg[$j]}) {
									print CNV2 $Chrom."\t".$Regions[$nextReg[$j]]{"start"}."\t".$Regions[$nextReg[$j]]{"end"}."\t".($Regions[$nextReg[$j]]{"end"}-$Regions[$nextReg[$j]]{"start"})." bp\t".$Regions[$nextReg[$j]]{"label"}."\t".($regionIndice{$nextReg[$j]}+1)."\t".$Results{$patient}{$nextReg[$j]}."\t".sprintf("%.3f",$Regions[$nextReg[$j]]{$patient}{"depth_ratio"}{$norm})."\t".$Regions[$nextReg[$j]]{"nb_CNV"}{$Results{$patient}{$nextReg[$j]}};
									if (@cnvFields) {
										my$txt = printCNVfields($normByGender,\@cnvFields,\%{ $Regions[$nextReg[$j]] },\%{ $Patients{$patient} });
										print CNV2 "$txt";
										}
									print CNV2 "\n";
									}
								}
							else	{ 
								$Result2{$patient}{$nextReg[$j]} = "NA";
								if (exists $regionIndice{$nextReg[$j]}) {
									print CNV2 $Chrom."\t".$Regions[$nextReg[$j]]{"start"}."\t".$Regions[$nextReg[$j]]{"end"}."\t".($Regions[$nextReg[$j]]{"end"}-$Regions[$nextReg[$j]]{"start"})." bp\t".$Regions[$nextReg[$j]]{"label"}."\t".($regionIndice{$nextReg[$j]}+1)."\t";
									if (exists $Regions[$nextReg[$j]]{"Appel"}) { print CNV2 "NA\t"; }
									else { print CNV2 "no\t"; }
									if (exists $Regions[$nextReg[$j]]{$patient}{"depth_ratio"}{$norm}) {
										print CNV2 sprintf("%.3f",$Regions[$nextReg[$j]]{$patient}{"depth_ratio"}{$norm})."\t".$Regions[$nextReg[$j]]{"nb_CNV"}{$Results{$patient}{$nextReg[0]}};
										push(@dirtyCNV, $Regions[$nextReg[$j]]{$patient}{"depth_ratio"}{$norm});
										}
									else { print CNV2 "na\tna"; }
									if (@cnvFields) {
										my$txt = printCNVfields($normByGender,\@cnvFields,\%{ $Regions[$nextReg[$j]] },\%{ $Patients{$patient} });
										print CNV2 "$txt";
										}
									print CNV2 "\n";
									}
								}
							}
						print CNV2 "\n";
						my$cleanAverage = 0;
						foreach (@cleanCNV) { $cleanAverage += $_; }
						$cleanAverage /= scalar@cleanCNV;
						my$dirtyAverage = 0;
						foreach (@dirtyCNV) { $dirtyAverage += $_; }
						$dirtyAverage /= scalar@dirtyCNV;
						print CNV1 $Chrom.":".$Regions[$nextReg[0]]{"start"}."-".$Regions[$nextReg[$cnvOK]]{"end"}."\t".$Results{$patient}{$nextReg[0]}."\t".scalar@cleanCNV."\t".sprintf("%.3f",$cleanAverage)."\t".($cnvOK+1)."\t".sprintf("%.3f",$dirtyAverage)."\t$overlapCNT\n";
						$Result3{$patient}{$nextReg[0]} = $Results{$patient}{$nextReg[0]};
						}
					}
				if ($i > 1) { $r += $i; }
				else { $r++; }
				}
			}
		}
	close CNV1;
	close CNV2;
	}


##print graph foreach Chrom
if ($graphByChr) {
	for my$patient (keys%Patients) {
		graphChr2($outdir,$norm,$seuil_deletion,$seuil_duplication,$patient,$maxDepthGraph,\%Patients,\@Regions,\@ChrOrder,\%regionOrder,\%Result2);
		}
	}

##summary
open (OUT,">$outdir/CNV.summary.txt") or die "Pb lors de l'ecriture du fichier sortie $!\n";
print OUT "CNV analysis\n";
print OUT "\nparameters:
	intervals normalization: $norm
	intervals detected as deletion if normalized depth below $seuil_deletion of $norm
	intervals detected as duplication if normalized depth above $seuil_duplication of $norm\n";
if ($seuil_patient <1) { print OUT "\tsamples kept if NO more than ".(100*$seuil_patient)."% of intervals are CNVs\n"; }
if ($seuil_region <1) { print OUT "\tintervals kept if NO more than ".(100*$seuil_region)."% of samples are CNVs\n"; }
if ($seuil_cov) { print OUT "\tintervals discarded if at least one sample covered less than $seuil_cov\n"; }
if ($minDP) { print OUT "\tintervals discarded if mean/median depth less than $minDP\n"; }
if ($minCNV) { print OUT "\tCNV calling requires at least $minCNV consecutive called intervals\n"; }
if ($maxNonCNV) { print OUT "\tCNV calling tolerates NO more than $maxNonCNV non-CNV inside\n"; }
if ($normByGender eq "all") { print OUT "\tinterval depth normalization by gender for all chromosomes\n"; }
elsif ($normByGender eq "gono") { print OUT "\tinterval depth normalization by gender for sex chromosomes only\n"; }
if ($Ref eq "mean") { print OUT "\teach sample normalized by the mean of all region average depths\n"; }
else { print OUT "\teach sample normalized by the sum of sequenced bp\n"; }
if ($RefByGender) { print OUT "\t(during sample normalization, depths within chrX are doubled for males, and those within chrY are skipped)\n"; }
else { print OUT "\t(during sample normalization, depths from all chr are taken, whatever the sex)\n"; }

if ($sexTxt) {
	print OUT "\ngender determination:\n$sexTxt\n";
	}

print OUT "\nintervals : ".scalar@Regions."\n";
print OUT "\nintervals discarded:\n";
my$N_noData=0; my$N_CNV_Recurrent=0; my$N_lowCov=0;
foreach my$r (0..$#Regions) {
	if(defined $Regions[$r]{"Appel"}) {
		if ($Regions[$r]{"Appel"} eq "No_Data") { $N_noData++; }
		elsif ($Regions[$r]{"Appel"} eq "CNV_Recurrent") { $N_CNV_Recurrent++; }
		elsif ($Regions[$r]{"Appel"} eq "Couverture_Faible")  { $N_lowCov++; }
		}
	}
print OUT "\tnot enough data: $N_noData\n";
print OUT "\tCNV_Recurrent: $N_CNV_Recurrent\n";
if ($minDP || $seuil_cov) { print OUT "\tlow coverage :  $N_lowCov\n"; }

print OUT "\nsamples discarded:\n";
my$N_discard=0;
foreach my$patient (sort(keys%Patients)) {
		if ($Patients{$patient}{"ecarte"}) { print OUT "\t$Patients{$patient}{ID}\n"; $N_discard++; }
		}
unless ($N_discard) { print OUT "\tnone\n"; }

print OUT "\nResults
patient\tCNV\tnber\n";
foreach my$patient (sort(keys%Patients)) {
	unless($Patients{$patient}{"ecarte"}) {
		my$N_dup=0; my$N_del=0;
		foreach (keys%{ $Result3{$patient} }) {
			if ($Result3{$patient}{$_} eq "DUP") { $N_dup++; }
			elsif ($Result3{$patient}{$_} eq "DEL") { $N_del++; }
			}
		print OUT "\n$Patients{$patient}{ID}: 
		DUP: $N_dup
		DEL: $N_del\n";
		}
	}
close OUT;


}


####################

sub norm {
my($norm,$h1)=@_;
my@allDepth=@$h1;
my$normDepth;
if ($norm eq "med" || $norm eq "min" || $norm eq "max") {
	@allDepth = sort{$a<=>$b}@allDepth;
	if ($norm eq "min") { $normDepth = $allDepth[0]; }
	elsif ($norm eq "max") { $normDepth = $allDepth[-1]; }
	else {
		#odd?
		if(scalar@allDepth%2) 
			{ $normDepth = $allDepth[int(scalar@allDepth/2)]; }
		#even
		else { $normDepth = ( $allDepth[int(scalar@allDepth/2)-1] + $allDepth[int(scalar@allDepth/2)] )/2; }
		}
	}

else {
	foreach (@allDepth)
		{ $normDepth += $_; } 
	$normDepth /= scalar@allDepth;
	}

return($normDepth);
}


####################

sub printCNVfields {
my($normByGender,$h1,$h2,$h3) = @_;
my@cnvFields = @$h1;
my%Regions = %$h2;
my%Patient = %$h3;
my$txt = "";
foreach my$val (@cnvFields) {
	if ($Regions{"normByR_depth"}{$val}) {
		if (!$normByGender || ($normByGender && ($normByGender eq "gono" && $Regions{"Chrom"} !~ m/^chr[XY]$|^[XY]$/))) {
			$txt .= "\t".sprintf("%.1f",$Regions{"normByR_depth"}{$val});
			}
		else {
			if ($Patient{"Sexe"} eq "F") {
				$txt .= "\t".sprintf("%.1f",$Regions{"normByR_depth_fem"}{$val});
				}
			else {
				$txt .= "\t".sprintf("%.1f",$Regions{"normByR_depth_males"}{$val});
				}
			}
		}
	else { $txt .= "\tna"; }
	}
return ($txt);
}


####################

sub mergeConsecutiveCNV {
my($maxI,$maxNonCNV,$orderedCNV,$regionIndice,$h1,$h2,$h3)=@_;
my@regionOrder = @$h1;
my%Results = %$h2;
my@Regions = @$h3;
my$ok=1; my$i=0; my$cnvOK=0; my$nonCNV=0; my$nonCNVtot=0;
my@nextReg=($orderedCNV);
while ($ok) {
	$i++;
	if (exists $regionOrder[($regionIndice + $i)]) {
		push(@nextReg, $regionOrder[($regionIndice + $i)]);
		if (exists $Results{$nextReg[$i]}) {
			## same CNV type ?
			if($Results{$nextReg[$i]} eq $Results{$nextReg[0]}) {
				$cnvOK = $i;
				#$nonCNVtot += $nonCNV;
				$nonCNV=0;		##to reset
				}
			else { $ok=0; }
			}
		else {
			if (!exists $Regions[$nextReg[$i]]{"Appel"}) { 
				if ($maxNonCNV) {
					$nonCNV++;
					if ($nonCNV > $maxNonCNV) { $ok=0; }
					}
				else { $ok=0; }
				}
			}
		}
	else { $ok=0; }
	if ($maxI && $i>$maxI) { $ok=0; }
	}
return($i,$cnvOK,$nonCNVtot,\@nextReg);
}


####################

sub graphChr1 {

my($nGraf,$outdir,$norm,$seuil_deletion,$seuil_duplication,$file,$maxDepthGraph,$h1,$h2,$h3,$h4,$h5,$h6,$h7)= @_;
my@Files = @$h1;
my%sampleName = %$h2;
my@Regions = @$h3;
my@ChromOrder = @$h4;
my%regionOrder = %$h5;
my%Patients = %$h6;
my%Results = %$h7;

print "cmdR, for sample $sampleName{$file}\n";

my$normGraf = $norm;
if ($norm eq "std") { $normGraf = "moy"; }

my$Nbr_Chr= scalar(keys%regionOrder);
if ($nGraf eq "max") { $nGraf = $Nbr_Chr; }

my$maxX=0;
foreach my$Chrom (keys%regionOrder) {	
	if (scalar@{ $regionOrder{$Chrom} } > $maxX) { $maxX = scalar@{ $regionOrder{$Chrom} }; }
	}

my$cmdR = "";
my$c=1; #chr iteration
my$n=1; #chr iteration, stepped back to 0 each time a graph is done
my$N=1; #graph iteration
my$maxGeneLab=200;
my$maxGeneSep=500;
#my$maxLabLen=25;

foreach my$Chrom (@ChromOrder) { 
	
	$Chrom =~ s/^chr//i;
	
	if (exists $regionOrder{$Chrom}) {

		my$maxYsup=$seuil_duplication; my$maxYinf=$seuil_deletion;	
		foreach my$region (@{ $regionOrder{$Chrom} }) {
			foreach my$f (@Files) {
				if (exists$Regions[$region]{$f}{"depth_ratio"}{$normGraf}) {
					if ($Regions[$region]{$f}{"depth_ratio"}{$normGraf} > $maxYsup)
						{ $maxYsup = $Regions[$region]{$f}{"depth_ratio"}{$normGraf}; }
					if ($normGraf eq "std") {
						if ($Regions[$region]{$f}{"depth_ratio"}{$normGraf} < $maxYinf)
							{ $maxYinf = $Regions[$region]{$f}{"depth_ratio"}{$normGraf}; }
						}
					}
				}
			}
		if ($maxDepthGraph && $maxYsup > $maxDepthGraph) { $maxYsup = $maxDepthGraph; }
		if ($normGraf eq "std") {
			if ($maxDepthGraph && $maxYinf < (-$maxDepthGraph)) { $maxYinf = (-$maxDepthGraph); }
			}

		if ($normGraf eq "std") { $cmdR .= "par(fig=c(0,1,".(1-(($n-0.05)/$nGraf)).",".(1-(($n-0.95)/$nGraf))."), new=TRUE)
	plot (c(0,0), xlim=c(0,$maxX), ylim=c($maxYinf,$maxYsup), type =\"n\", main=\"chrom: $Chrom\", xlab=\"\", ylab=\"depth_ratio_to_$normGraf\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n"; }
		else { $cmdR .= "par(fig=c(0,1,".(1-(($n-0.05)/$nGraf)).",".(1-(($n-0.95)/$nGraf))."), new=TRUE)
	plot (c(0,0), xlim=c(0,$maxX), ylim=c(0,$maxYsup), type =\"n\", main=\"chrom: $Chrom\", xlab=\"\", ylab=\"depth_ratio_to_$normGraf\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n"; }

		my$Nbr_Reg = scalar@{ $regionOrder{$Chrom} };

		##gene separations
		my$currentGene=""; my$tmpTxt=""; my$Nbr_gene=0;
		for (my$r=0;$r<$Nbr_Reg;$r++) {
			if ($Regions[$regionOrder{$Chrom}[$r]]{"Gene"}) {
				if ($Regions[$regionOrder{$Chrom}[$r]]{"Gene"} ne "NA" && $Regions[$regionOrder{$Chrom}[$r]]{"Gene"} ne $currentGene)  {
					$tmpTxt .= "abline(v=".($r+0.5).", col=\"blue\", lty = \"dotted\", lwd=1)\n";
					$currentGene = $Regions[$regionOrder{$Chrom}[$r]]{"Gene"};
					$Nbr_gene++;
					}
				}
			}
		if ($Nbr_gene < $maxGeneSep) { $cmdR .= $tmpTxt; }

		my$Nbr_CNV=0;
		for (my$r=0;$r<$Nbr_Reg;$r++) {
			if (exists$Results{$file}{$regionOrder{$Chrom}[$r]})
				{ $Nbr_CNV++; }
			}

		##region labels
		my@printReg=();
		if ($maxX < $maxGeneLab) {
			##in grey if invalid
			for (my$r=0;$r<$Nbr_Reg;$r++) { 
				if (defined $Regions[$regionOrder{$Chrom}[$r]]{"Appel"}) { push(@printReg,$r) ; }
				}
			if (@printReg) {
				$cmdR .= "axis(1, at=c(";
				foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
				chop $cmdR;
				$cmdR .="), labels=c(";
				foreach my$r (@printReg) { $cmdR .= "\"".$Regions[$regionOrder{$Chrom}[$r]]{"label"}."\","; }
				chop $cmdR;
				$cmdR .= "), col.axis=\"darkgrey\", las=2";
				if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
				else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_Reg))." )\n"; }
				}
			##in black if valid
			@printReg=();
			for (my$r=0;$r<$Nbr_Reg;$r++) { 
				if(!defined $Regions[$regionOrder{$Chrom}[$r]]{"Appel"}) { push(@printReg,$r) ; }
				}
			if (@printReg) {
				$cmdR .= "axis(1, at=c(";
				foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
				chop $cmdR;
				$cmdR .="), labels=c(";
				foreach my$r (@printReg) { $cmdR .= "\"".$Regions[$regionOrder{$Chrom}[$r]]{"label"}."\","; }
				chop $cmdR;
				$cmdR .= "), col.axis=\"black\", las=2";
				if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
				else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_Reg))." )\n"; }
				}
			}
		else {
			##in grey if invalid; only ticks
			for (my$r=0;$r<$Nbr_Reg;$r++) { 
				if (defined $Regions[$regionOrder{$Chrom}[$r]]{"Appel"}) { push(@printReg,$r) ; }
				}
			if (@printReg && scalar@printReg<$maxGeneSep) {
				$cmdR .= "axis(1, at=c(";
				foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
				chop $cmdR;
				$cmdR .= "), labels = FALSE, col.ticks=\"darkgrey\")\n";
				}
			##in black if valid
			@printReg=();
			if ($Nbr_gene < $maxGeneLab) { #&& ($Nbr_gene+$Nbr_CNV)>=$maxGeneLab) {
				$currentGene="";
				for (my$r=0;$r<$Nbr_Reg;$r++) {
					if ($Regions[$regionOrder{$Chrom}[$r]]{"Gene"}) {
						if ($Regions[$regionOrder{$Chrom}[$r]]{"Gene"} ne "NA" && $Regions[$regionOrder{$Chrom}[$r]]{"Gene"} ne $currentGene)  { 
							push(@printReg,$r); 
							$currentGene = $Regions[$regionOrder{$Chrom}[$r]]{"Gene"};
							}
						}
					}
				}
			#elsif ($Nbr_gene<$maxGeneLab && ($Nbr_gene+$Nbr_CNV)<$maxGeneLab) {
			#	for (my$r=0;$r<$Nbr_Reg;$r++) {
			#		if ( ($Regions[$RegionOrder{$Chrom}[$r]]{"Gene"} && $Regions[$RegionOrder{$Chrom}[$r]]{"Gene"} ne "NA" && $Regions[$RegionOrder{$Chrom}[$r]]{"Gene"} ne $currentGene) || (exists$Results{$file}{$RegionOrder{$Chrom}[$r]}) ) {
			#			push(@printReg,$r);
			#			$currentGene = $Regions[$RegionOrder{$Chrom}[$r]]{"Gene"};
			#			}
			#		}
			#	}
			if (@printReg) {
				$cmdR .= "axis(1, at=c(";
				foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
				chop $cmdR;
				$cmdR .="), labels=c(";
				foreach my$r (@printReg) { $cmdR .= "\"".$Regions[$regionOrder{$Chrom}[$r]]{"Gene"}."\","; }
				chop $cmdR;
				$cmdR .= "), col.axis=\"black\", las=2";
				if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
				else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_gene))." )\n"; }
				}
			##in red if CNV; only ticks
			@printReg=();
			for (my$r=0;$r<$Nbr_Reg;$r++) {
				if (exists $Results{$file}{$regionOrder{$Chrom}[$r]}) { push(@printReg,$r) ; }
				}
			if (@printReg && scalar@printReg<$maxGeneSep) {
				$cmdR .= "axis(1, at=c(";
				foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
				chop $cmdR;
				$cmdR .= "), labels = FALSE, col.ticks=\"red\")\n";
				}
			
			
			}

		##all not target sample lines (black):
		for (my$f=0;$f<scalar@Files;$f++) {
			unless ($Files[$f] eq $file) {
				my$r1=0;
				while ($r1<$Nbr_Reg) {
					my$r2=$r1;
					while ($r2<$Nbr_Reg) {
						if (exists $Regions[$regionOrder{$Chrom}[$r2]]{$Files[$f]}{"depth_ratio"}{$normGraf}) { $r2++;}
						else { last; }
						}
					if (($r2-1) > $r1) {
						$cmdR .= "lines( c(";
						for (my$r=$r1;$r<$r2;$r++)
							{ $cmdR .= ($r+1).","; }
						chop $cmdR;
						$cmdR .= "), c(";
						for (my$r=$r1;$r<$r2;$r++)
							{ $cmdR .= $Regions[$regionOrder{$Chrom}[$r]]{$Files[$f]}{"depth_ratio"}{$normGraf}.","; }
						chop $cmdR;
						$cmdR .= "), type =\"l\", lwd=1, col=\"black\")\n";
						}
					elsif (($r2-1) == $r1) {
						$cmdR .= "lines( c(".($r1+1)."), c(".$Regions[$regionOrder{$Chrom}[$r1]]{$Files[$f]}{"depth_ratio"}{$normGraf}."), type =\"p\", lwd=1, col=\"black\")\n";
						}
					$r1 = ($r2+1);
					}
				}
			}

		##threshold lines:
		if ($norm eq "std" && $normGraf eq "moy") {
			foreach my$gender ("normByR_depth","normByR_depth_fem","normByR_depth_males") {
				my$r1=0;
				while ($r1<$Nbr_Reg) {
					my$r2=$r1;
					while ($r2<$Nbr_Reg) {
						if ($Regions[$regionOrder{$Chrom}[$r2]]{$gender}{"moy"}) { $r2++;}
						else { last; }
						}
					if (($r2-1) > $r1) {
						foreach my$fold ($seuil_duplication,$seuil_deletion) {
							$cmdR .= "lines( c(";
							for (my$r=$r1;$r<$r2;$r++)
								{ $cmdR .= ($r+1).","; }
							chop $cmdR;
							$cmdR .= "), c(";
							for (my$r=$r1;$r<$r2;$r++) {
								if ( $Regions[$regionOrder{$Chrom}[$r]]{$gender}{"moy"} )
									{ $cmdR .= (1+$fold*($Regions[$regionOrder{$Chrom}[$r]]{$gender}{"std"} / $Regions[$regionOrder{$Chrom}[$r]]{$gender}{"moy"})).","; }
								else 	{ $cmdR .= "1,"; }
								}
							chop $cmdR;
							$cmdR .= "), type =\"l\", lwd=1, col=\"darkgrey\")\n";
							}
						}
					elsif (($r2-1) == $r1) {
						$cmdR .= "lines( c(".($r1+1)."), c(".(1+($Regions[$regionOrder{$Chrom}[$r1]]{$gender}{"std"} / $Regions[$regionOrder{$Chrom}[$r1]]{$gender}{"moy"}))."), type =\"p\", lwd=1, col=\"darkgrey\")\n";
						$cmdR .= "lines( c(".($r1+1)."), c(".(1-($Regions[$regionOrder{$Chrom}[$r1]]{$gender}{"std"} / $Regions[$regionOrder{$Chrom}[$r1]]{$gender}{"moy"}))."), type =\"p\", lwd=1, col=\"darkgrey\")\n";
						}
					$r1 = ($r2+1);
					}
				}
			}
		else {
			$cmdR .= "abline(h=$seuil_deletion, col=\"darkgrey\", lty = \"dashed\", lwd=1)\n";
			$cmdR .= "abline(h=$seuil_duplication, col=\"darkgrey\", lty = \"dashed\", lwd=1)\n";
			}

		##target sample line (green)
		my$r1=0;
		while ($r1<$Nbr_Reg) {
			my$r2=$r1;
			while ($r2<$Nbr_Reg) {
				if (exists$Regions[$regionOrder{$Chrom}[$r2]]{$file}{"depth_ratio"}{$normGraf}) { $r2++;}
				else { last; }
				}
			if (($r2-1) > $r1) {
				$cmdR .= "lines( c(";
				for (my$r=$r1;$r<$r2;$r++)
					{ $cmdR .= ($r+1).","; }
				chop $cmdR;
				$cmdR .= "), c(";
				for (my$r=$r1;$r<$r2;$r++)
					{ $cmdR .= $Regions[$regionOrder{$Chrom}[$r]]{$file}{"depth_ratio"}{$normGraf}.","; }
				chop $cmdR;
				$cmdR .= "), type =\"l\", lwd=1, col=\"green\")\n";
				}
			elsif (($r2-1) == $r1) {
				$cmdR .= "lines( c(".($r1+1)."), c(".$Regions[$regionOrder{$Chrom}[$r1]]{$file}{"depth_ratio"}{$normGraf}."), type =\"p\", lwd=1, col=\"green\")\n";
				}
			$r1 = ($r2+1);
			}
		##points for CNVs
		my$points .= "points( c(";
		for (my$r=0;$r<$Nbr_Reg;$r++) {
			if (exists$Results{$file}{$regionOrder{$Chrom}[$r]} && exists$Regions[$regionOrder{$Chrom}[$r]]{$file}{"depth_ratio"}{$normGraf})
				{ $points .= ($r+1).","; }
			}
		chop $points;
		$points .= "), c(";
		foreach my$region (@{ $regionOrder{$Chrom} }) {
			if (exists$Results{$file}{$region} && exists$Regions[$region]{$file}{"depth_ratio"}{$normGraf})
				{ $points .= $Regions[$region]{$file}{"depth_ratio"}{$normGraf}.","; }
			}
		chop $points;
		$points .= "), type =\"p\", pch = 16, lwd=2, col=\"red\")\n";
		if ($Nbr_CNV) { $cmdR .= $points; }

		#if ($normGraf eq "std") { $cmdR .= "abline(h=0, col=\"darkgrey\", lty = \"dashed\", lwd=1)\n"; }
		#else { $cmdR .= "abline(h=1, col=\"darkgrey\", lty = \"dashed\", lwd=1)\n"; }
		$cmdR .= "abline(h=1, col=\"darkgrey\", lty = \"dashed\", lwd=1)\n";


		if ($c==$Nbr_Chr || $n==$nGraf) {
			open (CMDR, ">$outdir/$sampleName{$file}\_temp.R") || die;
			print CMDR "#!/usr/bin/env Rscript\n\n" ;
			if ($nGraf==$Nbr_Chr) { print CMDR "pdf(\"".$outdir."/CNV_$sampleName{$file}.pdf\", width=11.69, height=".($nGraf*3).")\npar(mfrow=c($nGraf,1))\n"; }
			else {
				if ($N>1) { print CMDR "pdf(\"".$outdir."/CNV_$sampleName{$file}\_$N.pdf\", width=11.69, height=".($nGraf*3).")\npar(mfrow=c($nGraf,1))\n"; }
				else { print CMDR "pdf(\"".$outdir."/CNV_$sampleName{$file}\_$N.pdf\", width=11.69, height=".($n*3).")\npar(mfrow=c($n,1))\n"; }
				}
			print CMDR "$cmdR";
			print CMDR "title(main=\"sample: $sampleName{$file}";
			if ($Patients{$file}{"ecarte"}) { print CMDR " (invalid)\", col.main=\"red\""; }
			else { print CMDR "\""; }
			print CMDR ", outer=TRUE, line=-2, cex.main=2)\n";
			print CMDR "dev.off();\n";
			close CMDR;
			system "Rscript $outdir/$sampleName{$file}\_temp.R";
			unlink "$outdir/$sampleName{$file}\_temp.R";
			$cmdR="";
			$n=0;
			$N++;
			}

		$c++;$n++;

		}	
	}
}


####################
sub graphChr2 {

my($outdir,$norm,$seuil_deletion,$seuil_duplication,$patient,$maxDepthGraph,$h1,$h2,$h3,$h4,$h5)= @_;
my%Patients = %$h1;
my@Regions = @$h2;
my@ChrOrder = @$h3;
my%regionOrder = %$h4;
my%Results = %$h5;	#$results{$patient}{$region}

my$normGraf = $norm;
if ($norm eq "std") { $normGraf = "moy"; }

my$maxX=0;
foreach my$Chrom (keys%regionOrder) {	
	if (scalar@{ $regionOrder{$Chrom} } > $maxX) { $maxX = scalar@{ $regionOrder{$Chrom} }; }
	}

my$Nbr_Chr= scalar(keys%regionOrder);
print "cmdR,  for $Patients{$patient}{ID}\n";
open (CMDR, ">$outdir/$Patients{$patient}{ID}\_temp.R") || die;
print CMDR "#!/usr/bin/env Rscript\n\n" ;
print CMDR "pdf(\"".$outdir."/CNV_$Patients{$patient}{ID}.pdf\", width=11.69, height=".($Nbr_Chr*3).")\n
par(mfrow=c($Nbr_Chr,1))\n";

my$cmdR = "";
my$c=0; #chr iteration
my$maxGeneLab=200;
my$maxGeneSep=500;
#my$maxLabLen=25;

for my$Chrom (@ChrOrder) {

	my$maxYsup=$seuil_duplication; my$maxYinf=$seuil_deletion;	
	foreach my$region (@{ $regionOrder{$Chrom} }) {
		for (my$p=0;$p<scalar(keys%Patients);$p++) {
			if (exists $Regions[$region]{$p}{"depth_ratio"}{$normGraf}) {
				if ($Regions[$region]{$p}{"depth_ratio"}{$normGraf} > $maxYsup)
					{ $maxYsup = $Regions[$region]{$p}{"depth_ratio"}{$normGraf}; }
				if ($normGraf eq "std") {
					if ($Regions[$region]{$p}{"depth_ratio"}{$normGraf} < $maxYinf)
						{ $maxYinf = $Regions[$region]{$p}{"depth_ratio"}{$normGraf}; }
					}
				}
			}
		}
	if ($maxDepthGraph && $maxYsup > $maxDepthGraph) { $maxYsup = $maxDepthGraph; }
	if ($normGraf eq "std") {
		if ($maxDepthGraph && $maxYinf < (-$maxDepthGraph)) { $maxYinf = (-$maxDepthGraph); }
		}

	my$ChrName = $Chrom; $ChrName =~ s/^chr//;

	if ($normGraf eq "std") {
		$cmdR .= "par(fig=c(0,1,".(1-(($c+0.95)/$Nbr_Chr)).",".(1-(($c+0.05)/$Nbr_Chr))."), new=TRUE)
plot (c(0,0), xlim=c(0,$maxX), ylim=c($maxYinf,$maxYsup), type =\"n\", main=\"chrom: $ChrName\", xlab=\"\", ylab=\"depth_ratio_to_$normGraf\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n";
		}
	else {
		$cmdR .= "par(fig=c(0,1,".(1-(($c+0.95)/$Nbr_Chr)).",".(1-(($c+0.05)/$Nbr_Chr))."), new=TRUE)
plot (c(0,0), xlim=c(0,$maxX), ylim=c(0,$maxYsup), type =\"n\", main=\"chrom: $ChrName\", xlab=\"\", ylab=\"depth_ratio_to_$normGraf\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n";
		}

	my$Nbr_Reg = scalar@{ $regionOrder{$Chrom} };

	#gene vertical separations
	my$currentGene=""; my$tmpTxt=""; my$Nbr_gene=0;
	for (my$r=0;$r<$Nbr_Reg;$r++) {
		if ($Regions[$regionOrder{$Chrom}[$r]]{"Gene"} ne "NA" && $Regions[$regionOrder{$Chrom}[$r]]{"geneID"} ne $currentGene) {
			$tmpTxt .= "abline(v=".($r+0.5).", col=\"blue\", lty = \"dotted\", lwd=1)\n";
			$currentGene = $Regions[$regionOrder{$Chrom}[$r]]{"geneID"};
			$Nbr_gene++;
			}
		}
	if ($Nbr_gene < $maxGeneSep) { $cmdR .= $tmpTxt; }

	my$Nbr_CNV=0;
	for (my$r=0;$r<$Nbr_Reg;$r++) {
		if (exists$Results{$patient}{$regionOrder{$Chrom}[$r]})
			{ $Nbr_CNV++; }
		}

	#x labels
	my@printReg=();
	if ($maxX < $maxGeneLab) {
		##in grey if invalid
		for (my$r=0;$r<$Nbr_Reg;$r++) { 
			if (defined $Regions[$regionOrder{$Chrom}[$r]]{"Appel"}) { push(@printReg,$r) ; }
			}
		if (@printReg) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
			chop $cmdR;
			$cmdR .="), labels=c(";
			foreach my$r (@printReg) { $cmdR .= "\"".$Regions[$regionOrder{$Chrom}[$r]]{"label"}."\","; }
			chop $cmdR;
			$cmdR .= "), col.axis=\"darkgrey\", las=2";
			if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
			else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_Reg))." )\n"; }
			}
		##in black if valid
		@printReg=();
		for (my$r=0;$r<$Nbr_Reg;$r++) { 
			if(!defined $Regions[$regionOrder{$Chrom}[$r]]{"Appel"}) { push(@printReg,$r) ; }
			}
		if (@printReg) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
			chop $cmdR;
			$cmdR .="), labels=c(";
			foreach my$r (@printReg) { $cmdR .= "\"".$Regions[$regionOrder{$Chrom}[$r]]{"label"}."\","; }
			chop $cmdR;
			$cmdR .= "), col=\"black\", las=2";
			if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
			else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_Reg))." )\n"; }
			}
		}
	else {
		##in grey if invalid; only ticks
		for (my$r=0;$r<$Nbr_Reg;$r++) { 
			if (defined $Regions[$regionOrder{$Chrom}[$r]]{"Appel"}) { push(@printReg,$r) ; }
			}
		if (@printReg && scalar@printReg<$maxGeneSep) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
			chop $cmdR;
			$cmdR .= "), labels = FALSE, col.ticks=\"darkgrey\")\n";
			}
		##in black if valid
		@printReg=();
		if ($Nbr_gene < $maxGeneLab) {	# && ($Nbr_gene+$Nbr_CNV)>=$maxGeneLab) {
			$currentGene="";
			for (my$r=0;$r<$Nbr_Reg;$r++) {
				if ($Regions[$regionOrder{$Chrom}[$r]]{"Gene"} ne "NA" && $Regions[$regionOrder{$Chrom}[$r]]{"geneID"} ne $currentGene) {
					push(@printReg,$r);
					$currentGene = $Regions[$regionOrder{$Chrom}[$r]]{"geneID"};
					}
				}
			}
		#elsif ($Nbr_gene<$maxGeneLab && ($Nbr_gene+$Nbr_CNV)<$maxGeneLab) {
		#	for (my$r=0;$r<$Nbr_Reg;$r++) {
		#		if ( ($Regions[$regionOrder{$Chrom}[$r]]{"Gene"} ne "NA" && $Regions[$regionOrder{$Chrom}[$r]]{"geneID"} ne $currentGene) || (exists$results{$patient}{$regionOrder{$Chrom}[$r]}) ) {
		#			push(@printReg,$r);
		#			$currentGene = $Regions[$regionOrder{$Chrom}[$r]]{"geneID"};
		#			}
		#		}
		#	}
		if (@printReg) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
			chop $cmdR;
			$cmdR .="), labels=c(";
			foreach my$r (@printReg) { $cmdR .= "\"".$Regions[$regionOrder{$Chrom}[$r]]{"label"}."\","; }
			chop $cmdR;
			$cmdR .= "), col.axis=\"black\", las=2";
			if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
			else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_gene))." )\n"; }
			}
		##in red if CNV; only ticks
		@printReg=();
		for (my$r=0;$r<$Nbr_Reg;$r++) {
			if (exists $Results{$patient}{$regionOrder{$Chrom}[$r]}) { push(@printReg,$r) ; }
			}
		if (@printReg && scalar@printReg<$maxGeneSep) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
			chop $cmdR;
			$cmdR .= "), labels = FALSE, col.ticks=\"red\")\n";
			}
		}

	#all not target sample lines (black):
	for (my$p=0;$p<scalar(keys%Patients);$p++) {
		unless ($p == $patient) {

			my$r1=0;
			while ($r1<$Nbr_Reg) {
				my$r2=$r1;
				while ($r2<$Nbr_Reg) {
					if (exists $Regions[$regionOrder{$Chrom}[$r2]]{$p}{"depth_ratio"}{$normGraf}) { $r2++;}
					else { last; }
					}
				if (($r2-1) > $r1) {
					$cmdR .= "lines( c(";
					for (my$r=$r1;$r<$r2;$r++)
						{ $cmdR .= ($r+1).","; }
					chop $cmdR;
					$cmdR .= "), c(";
					for (my$r=$r1;$r<$r2;$r++)
						{ $cmdR .= $Regions[$regionOrder{$Chrom}[$r]]{$p}{"depth_ratio"}{$normGraf}.","; }
					chop $cmdR;
					$cmdR .= "), type =\"l\", lwd=1, col=\"black\")\n";
					}
				elsif (($r2-1) == $r1) {
					$cmdR .= "lines( c(".($r1+1)."), c(".$Regions[$regionOrder{$Chrom}[$r1]]{$p}{"depth_ratio"}{$normGraf}."), type =\"p\", lwd=1, col=\"black\")\n";
					}
				$r1 = ($r2+1);
				}
			}
		}

	##threshold lines:
	if ($norm eq "std" && $normGraf eq "moy") {
		foreach my$gender ("normByR_depth","normByR_depth_fem","normByR_depth_males") {
			my$r1=0;
			while ($r1<$Nbr_Reg) {
				my$r2=$r1;
				while ($r2<$Nbr_Reg) {
					if ($Regions[$regionOrder{$Chrom}[$r2]]{$gender}{"moy"}) { $r2++; }
					else { last; }
					}
				if (($r2-1) > $r1) {
					#foreach my$fold (1,$seuil_duplication,(-1),$seuil_deletion) {
					foreach my$fold ($seuil_duplication,$seuil_deletion) {
						$cmdR .= "lines( c(";
						for (my$r=$r1;$r<$r2;$r++)
							{ $cmdR .= ($r+1).","; }
						chop $cmdR;
						$cmdR .= "), c(";
						for (my$r=$r1;$r<$r2;$r++) {
							if ( $Regions[$regionOrder{$Chrom}[$r]]{$gender}{"moy"} )
								{ $cmdR .= (1+$fold*($Regions[$regionOrder{$Chrom}[$r]]{$gender}{"std"} / $Regions[$regionOrder{$Chrom}[$r]]{$gender}{"moy"})).","; }
							else 	{ $cmdR .= "1,"; }
							}
						chop $cmdR;
						$cmdR .= "), type =\"l\", lwd=1, col=\"darkgrey\")\n";
						}
					}
				elsif (($r2-1) == $r1) {
					$cmdR .= "lines( c(".($r1+1)."), c(".(1+($Regions[$regionOrder{$Chrom}[$r1]]{$gender}{"std"} / $Regions[$regionOrder{$Chrom}[$r1]]{$gender}{"moy"}))."), type =\"p\", lwd=1, col=\"darkgrey\")\n";
					$cmdR .= "lines( c(".($r1+1)."), c(".(1-($Regions[$regionOrder{$Chrom}[$r1]]{$gender}{"std"} / $Regions[$regionOrder{$Chrom}[$r1]]{$gender}{"moy"}))."), type =\"p\", lwd=1, col=\"darkgrey\")\n";
					}
				$r1 = ($r2+1);
				}
			}
		}
	else {
		$cmdR .= "abline(h=$seuil_deletion, col=\"darkgrey\", lty = \"dashed\", lwd=1)\n";
		$cmdR .= "abline(h=$seuil_duplication, col=\"darkgrey\", lty = \"dashed\", lwd=1)\n";
		}

	#target sample line (green):
	my$r1=0;
	while ($r1<$Nbr_Reg) {
		my$r2=$r1;
		while ($r2<$Nbr_Reg) {
			if (exists$Regions[$regionOrder{$Chrom}[$r2]]{$patient}{"depth_ratio"}{$normGraf}) { $r2++;}
			else { last; }
			}
		if (($r2-1) > $r1) {
			$cmdR .= "lines( c(";
			for (my$r=$r1;$r<$r2;$r++)
				{ $cmdR .= ($r+1).","; }
			chop $cmdR;
			$cmdR .= "), c(";
			for (my$r=$r1;$r<$r2;$r++)
				{ $cmdR .= $Regions[$regionOrder{$Chrom}[$r]]{$patient}{"depth_ratio"}{$normGraf}.","; }
			chop $cmdR;
			$cmdR .= "), type =\"l\", lwd=1, col=\"green\")\n";
			}
		elsif (($r2-1) == $r1) {
			$cmdR .= "lines( c(".($r1+1)."), c(".$Regions[$regionOrder{$Chrom}[$r1]]{$patient}{"depth_ratio"}{$normGraf}."), type =\"p\", lwd=1, col=\"green\")\n";
			}
		$r1 = ($r2+1);
		}

	#points for CNVs
	my$points .= "points( c(";
	for (my$r=0;$r<$Nbr_Reg;$r++) {
		if (exists$Results{$patient}{$regionOrder{$Chrom}[$r]} && exists$Regions[$regionOrder{$Chrom}[$r]]{$patient}{"depth_ratio"}{$normGraf})
			{ $points .= ($r+1).","; }
		}
	chop $points;
	$points .= "), c(";
	foreach my$region (@{ $regionOrder{$Chrom} }) {
		if (exists$Results{$patient}{$region} && exists$Regions[$region]{$patient}{"depth_ratio"}{$normGraf})
			{ $points .= $Regions[$region]{$patient}{"depth_ratio"}{$normGraf}.","; }
		}
	chop $points;
	$points .= "), type =\"p\", pch = 16, lwd=2, col=\"red\")\n";
	if ($Nbr_CNV) { $cmdR .= $points; }

	#if ($normGraf eq "std") { $cmdR .= "abline(h=0, col=\"darkgrey\", lty = \"dashed\", lwd=1)\n"; }
	#else { $cmdR .= "abline(h=1, col=\"darkgrey\", lty = \"dashed\", lwd=1)\n"; }
	$cmdR .= "abline(h=1, col=\"darkgrey\", lty = \"dashed\", lwd=1)\n";


	$c++;
	}

#print "$cmdR\n";
print CMDR "$cmdR";
print CMDR "title(main=\"sample: $Patients{$patient}{ID}";
if ($Patients{$patient}{"ecarte"}) { print CMDR " (invalid)\", col.main=\"red\""; }
else { print CMDR "\""; }
print CMDR ", outer=TRUE, line=-2, cex.main=2)\n";
print CMDR "dev.off();\n";
close CMDR;
system "Rscript $outdir/$Patients{$patient}{ID}\_temp.R";
unlink "$outdir/$Patients{$patient}{ID}\_temp.R";


}

####################

1;

