package Bio::NGS::HCL::DeCovA::CNV_tool;

##version1.5

use strict;
use warnings;


#@hashSub = CNV_tool::CNV_detect(\@Files,\%sName2,"$outdir/CNV_analysis",$pThreshold,$nGraf,$fichier_sexe,\%CNV_opt,$mbq,$mmq,$refBedLines,\@ChromOrder,\%fai);
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
my$CNVgraph = $CNV_opt{"chromGraph"};

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
	$Patients{$file}{"Ref_Autosomes_Patient"} = 0;
	$Patients{$file}{"Ref_X_Patient"} = 0;
	$Patients{$file}{"Ref_Y_Patient"} = 0;
	$Patients{$file}{"Ref_Profondeur_Patient"} = 0;
	}

my$n_Regions = 0;
my$autosom_Regions = 0;
my$gonosom_Regions = 0;

#$bedLines[$i]{"allLine"} = $line;
#$bedLines[$i]{"Chr"} = $tab[0]; $bedLines[$i]{"Start"} = ($tab[1]+1); $bedLines[$i]{"End"} = $tab[2];
#$bedLines[$i]{$file}{"Cov"}
#$bedLines[$i]{$file}{"Mean"}
for my$r (0..$#Regions) {
	my $region_a_conserver = 0;
	if ($minCov) {
		foreach my$file (@Files) {
			if ( $Regions[$r]{$file}{"Cov"} >= $minCov)
				{ $region_a_conserver = 1; last; }
			}
		}
	else { $region_a_conserver = 1; }
	# si region à conserver, calcul Profondeur_Autosomes/gonosomes
	if ($region_a_conserver) {
		foreach my$file (@Files) {
			if ($Regions[$r]{"Chrom"} !~ /^(chr[XY]|[XY])$/) {
				$autosom_Regions++;
				if ($Ref eq "mean") { $Patients{$file}{"Ref_Autosomes_Patient"} += $Regions[$r]{$file}{"Mean"}; }
				}
			else {
				if ($Regions[$r]{"Chrom"} =~ /^(chrX|X)$/) {
					$gonosom_Regions++;
					if ($Ref eq "mean") { $Patients{$file}{"Ref_X_Patient"} += $Regions[$r]{$file}{"Mean"}; }
					}
				else {
					if ($Ref eq "mean") { $Patients{$file}{"Ref_Y_Patient"} += $Regions[$r]{$file}{"Mean"}; }
					}
				}
			}
		}
	# Si region mal couverte, on l'imprime dans un fichier POUBELLE
	else
		{ $Regions[$r]{"Appel"} = "Couverture_Faible"; }
	}


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
		$Patients{$file}{"Ref_Autosomes_Patient"} = `$cmd`;
		chomp $Patients{$file}{"Ref_Autosomes_Patient"};
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
		$Patients{$file}{"Ref_X_Patient"} = `$cmd`;
		chomp $Patients{$file}{"Ref_X_Patient"};
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
		$Patients{$file}{"Ref_Y_Patient"} = `$cmd`;
		chomp $Patients{$file}{"Ref_Y_Patient"};
		

		print "total sequenced bases for $sampleName{$file} : \n\tautoZ:".$Patients{$file}{"Ref_Autosomes_Patient"}."\n\tchrX:".$Patients{$file}{"Ref_X_Patient"}."\n\tchrY:".$Patients{$file}{"Ref_Y_Patient"}."\n";
		}

	##gender if no sex file
	print $sampleName{$file}.":\n";
	unless (exists $Patients{$file}{"Sexe"}) {
		if ($Patients{$file}{"Ref_X_Patient"} && $gonosom_Regions && $Patients{$file}{"Ref_Autosomes_Patient"} && $autosom_Regions) {
			print "\tautoZ: ".($Patients{$file}{"Ref_Autosomes_Patient"}/$autosom_Regions)."\n\tchrX: ".($Patients{$file}{"Ref_X_Patient"}/$gonosom_Regions)."\n";
			if ( ($Patients{$file}{"Ref_X_Patient"}/$gonosom_Regions) > (1.2*$Patients{$file}{"Ref_Autosomes_Patient"}/$autosom_Regions) ) { print "\t-> sexe ambigu\n"; $Patients{$file}{"Sexe"} = "F"; }
			elsif ( (($Patients{$file}{"Ref_X_Patient"}/$gonosom_Regions) <= (1.2*$Patients{$file}{"Ref_Autosomes_Patient"}/$autosom_Regions)) && (($Patients{$file}{"Ref_X_Patient"}/$gonosom_Regions) >= (0.8*$Patients{$file}{"Ref_Autosomes_Patient"}/$autosom_Regions)) ) { $Patients{$file}{"Sexe"} = "F"; }
			elsif ( (($Patients{$file}{"Ref_X_Patient"}/$gonosom_Regions) < (0.8*$Patients{$file}{"Ref_Autosomes_Patient"}/$autosom_Regions)) && (($Patients{$file}{"Ref_X_Patient"}/$gonosom_Regions) > (1.2*0.5*$Patients{$file}{"Ref_Autosomes_Patient"}/$autosom_Regions)) ) { print "\t-> sexe ambigu\n"; $Patients{$file}{"Sexe"} = "F"; }
			elsif ( (($Patients{$file}{"Ref_X_Patient"}/$gonosom_Regions) <= (1.2*0.5*$Patients{$file}{"Ref_Autosomes_Patient"}/$autosom_Regions)) && (($Patients{$file}{"Ref_X_Patient"}/$gonosom_Regions) >= (0.8*0.5*$Patients{$file}{"Ref_Autosomes_Patient"}/$autosom_Regions)) ) { $Patients{$file}{"Sexe"} = "H"; }
			else { print "\t-> sexe ambigu\n"; $Patients{$file}{"Sexe"} = "H"; }
			}
		else {
			print "\tno data available to determine the gender\n";
			$Patients{$file}{"Sexe"} = "F"; 
			}
		}
	print "\t-> $Patients{$file}{Sexe}\n";

	##Ref_Profondeur_Patient
	if ($RefByGender) {
		if ($Patients{$file}{"Sexe"} eq "H") {
			$Patients{$file}{"Ref_Profondeur_Patient"} = ($Patients{$file}{"Ref_X_Patient"} * 2) + $Patients{$file}{"Ref_Autosomes_Patient"}; } 
		else {
			$Patients{$file}{"Ref_Profondeur_Patient"} = $Patients{$file}{"Ref_X_Patient"} + $Patients{$file}{"Ref_Autosomes_Patient"}; }
		}
	else {
		$Patients{$file}{"Ref_Profondeur_Patient"} = $Patients{$file}{"Ref_X_Patient"} + $Patients{$file}{"Ref_Y_Patient"} + $Patients{$file}{"Ref_Autosomes_Patient"};
		}

	}

my$meanRef;
foreach my$file (@Files) {
	if ($Patients{$file}{"Ref_Profondeur_Patient"})
		{ $meanRef += $Patients{$file}{"Ref_Profondeur_Patient"}; }
	else { $Patients{$file}{"ecarte"} = 1; }
	}
$meanRef /= scalar@Files;
foreach my$file (@Files) {
	$Patients{$file}{"Ref_Profondeur_Patient"} /= $meanRef; 
	print "Ratio_Profondeur_Patient for $sampleName{$file} : ".$Patients{$file}{"Ref_Profondeur_Patient"}."\n";
	}


my $nb_parcours = 0;
my $continuer = 1;
my %Results;
##iterations
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
				$Patients{$file}{"Ref_Autosomes_Patient"} = 0;
				$Patients{$file}{"Ref_Gonosomes_Patient"} = 0;
				}
			for (my $r = 0 ; $r < $#Regions ; $r++) {
				if(!defined $Regions[$r]{"Appel"}) {
					foreach my$file (@Files) {
						if ($Regions[$r]{"Chrom"} !~ /^(chr[XY]|[XY])$/) {
							$Patients{$file}{"Ref_Autosomes_Patient"} += $Regions[$r]{$file}{"Mean"};
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
					$Patients{$file}{"Ref_Profondeur_Patient"} = ($Patients{$file}{"Ref_Gonosomes_Patient"} * 2) + $Patients{$file}{"Ref_Autosomes_Patient"}; } 
				else {
					$Patients{$file}{"Ref_Profondeur_Patient"} = $Patients{$file}{"Ref_Gonosomes_Patient"} + $Patients{$file}{"Ref_Autosomes_Patient"}; }
				}
			foreach my$file (@Files)
				{ $meanRef += $Patients{$file}{"Ref_Profondeur_Patient"}; }
			$meanRef /= scalar@Files;
			foreach my$file (@Files) {
				$Patients{$file}{"Ref_Profondeur_Patient"} /= $meanRef; 
				print "Ratio_Profondeur_Patient for $sampleName{$file} : ".$Patients{$file}{"Ref_Profondeur_Patient"}."\n";
				}
			}
		}

	# SECOND PARCOURS DES REGIONS #
	# PERMET DE PONDERER LA PROFONDEUR PAR LES AUTRES REGIONS (INTRA) ET ENTRE LES PATIENTS (INTER)
	for my$r (0..$#Regions) {

		print SORTIE $Regions[$r]{"allLine"}."\t";
		print SORTIE $r."\t";

		# SI LA REGION N'EST PAS "MOCHE"
		if (!(defined($Regions[$r]{"Appel"}))) {

			$nb_regions_conservees++;
			my $nb_deletions = 0;
			my $nb_duplications = 0;
			my $nb_evts = 0;
			my$prof_Moyenne_Inter = "";

			if ($normByGender) {

				@{ $Regions[$r]{"Prof_Region_Ponderee_Femmes"} } = ();
				@{ $Regions[$r]{"Prof_Region_Ponderee_Hommes"} } = ();
				# Pour chaque patient, un premier parcours permet la pondération intra-patient (en divisant la profondeur de la region par la profondeur totale du patient)
				foreach my$file (@Files) {

					if ($Patients{$file}{"Sexe"} eq "F" && !($Patients{$file}{"ecarte"})) {

						# Le controle permet de diviser uniquement par le nombre de patients restant apres l filtre, et non par l'ensemble des patients
						# Si la region concerne un gonosome et le patient est une femme, la référence est celle qui correspond aux autosomes + au chromosome X
						$Regions[$r]{$file}{"Prof_Region_Ponderee"} = $Regions[$r]{$file}{"Mean"} / $Patients{$file}{"Ref_Profondeur_Patient"};
						# On sauvegarde la profondeur ponderee intra patient, pour determiner ensuite la moyenne de ces profondeurs ponderees
						push( @{ $Regions[$r]{"Prof_Region_Ponderee_Femmes"} }, $Regions[$r]{$file}{"Prof_Region_Ponderee"} );


					} elsif ($Patients{$file}{"Sexe"} eq "H" && !($Patients{$file}{"ecarte"})) {

						# Le controle permet de divisier uniquement par le nombre de patients restant apres l filtre, et non par l'ensemble des patients
						# Si la region concerne un gonosome et le patient est un homme, la référence est celle qui correspond uniquement au chromosome X
						$Regions[$r]{$file}{"Prof_Region_Ponderee"} = $Regions[$r]{$file}{"Mean"} / $Patients{$file}{"Ref_Profondeur_Patient"};
						# On sauvegarde la profondeur ponderee intra patient, pour determiner ensuite la moyenne de ces profondeurs ponderees
						push( @{ $Regions[$r]{"Prof_Region_Ponderee_Hommes"} }, $Regions[$r]{$file}{"Prof_Region_Ponderee"} );

					}

				}

				# Nous calculons pour chaque region et chaque sexe la moyenne de profondeur (ponderee pour chaque patient par rapport aux autres regions)
				if ($normByGender eq "all") {
					if(@{ $Regions[$r]{"Prof_Region_Ponderee_Femmes"} }) {
						 $Regions[$r]{"Norm_Prof_Ponderees_Femmes"} = norm($norm,\@{ $Regions[$r]{"Prof_Region_Ponderee_Femmes"} });
						if ($norm eq "std") {
							my$sqtotal = 0;
							foreach (@{ $Regions[$r]{"Prof_Region_Ponderee_Femmes"} })
								{ $sqtotal += ($Regions[$r]{"Norm_Prof_Ponderees_Femmes"}-$_)**2; }
							$Regions[$r]{"Std_Prof_Ponderees_Femmes"} = ($sqtotal / (scalar@{ $Regions[$r]{"Prof_Region_Ponderee_Femmes"} }-1))**0.5;
						}
					}
					if(@{ $Regions[$r]{"Prof_Region_Ponderee_Hommes"} }) {
						$Regions[$r]{"Norm_Prof_Ponderees_Hommes"} = norm($norm,\@{ $Regions[$r]{"Prof_Region_Ponderee_Hommes"} });
						if ($norm eq "std") {
							my$sqtotal = 0;
							foreach (@{ $Regions[$r]{"Prof_Region_Ponderee_Hommes"} })
								{ $sqtotal += ($Regions[$r]{"Norm_Prof_Ponderees_Hommes"}-$_)**2; }
							$Regions[$r]{"Std_Prof_Ponderees_Hommes"} = ($sqtotal / (scalar@{ $Regions[$r]{"Prof_Region_Ponderee_Hommes"} }-1))**0.5;
						}
					}
				} else {
					if ($Regions[$r]{"Chrom"} !~ /^(chr[XY]|[XY])$/) {
						my@Autosomes=();
						if(@{ $Regions[$r]{"Prof_Region_Ponderee_Femmes"} }) {
							push(@Autosomes, @{ $Regions[$r]{"Prof_Region_Ponderee_Femmes"} });
						}
						if(@{ $Regions[$r]{"Prof_Region_Ponderee_Hommes"} }) {
							push(@Autosomes, @{ $Regions[$r]{"Prof_Region_Ponderee_Hommes"} });
						}
						if (@Autosomes) {
							$Regions[$r]{"Norm_Prof_Ponderees_Autosomes"} = norm($norm,\@Autosomes);
							if ($norm eq "std") {
								my$sqtotal = 0;
								foreach (@Autosomes)
									{ $sqtotal += ($Regions[$r]{"Norm_Prof_Ponderees_Autosomes"}-$_)**2; }
								$Regions[$r]{"Std_Prof_Ponderees_Autosomes"} = ($sqtotal / (scalar@Autosomes-1))**0.5;
							}
						}
					} else {
						if(@{ $Regions[$r]{"Prof_Region_Ponderee_Femmes"} }) {
							 $Regions[$r]{"Norm_Prof_Ponderees_Femmes"} = norm($norm,\@{ $Regions[$r]{"Prof_Region_Ponderee_Femmes"} });
							if ($norm eq "std") {
								my$sqtotal = 0;
								foreach (@{ $Regions[$r]{"Prof_Region_Ponderee_Femmes"} })
									{ $sqtotal += ($Regions[$r]{"Norm_Prof_Ponderees_Femmes"}-$_)**2; }
								$Regions[$r]{"Std_Prof_Ponderees_Femmes"} = ($sqtotal / (scalar@{ $Regions[$r]{"Prof_Region_Ponderee_Femmes"} }-1))**0.5;
							}
						}
						if(@{ $Regions[$r]{"Prof_Region_Ponderee_Hommes"} }) {
							$Regions[$r]{"Norm_Prof_Ponderees_Hommes"} = norm($norm,\@{ $Regions[$r]{"Prof_Region_Ponderee_Hommes"} });
							if ($norm eq "std") {
								my$sqtotal = 0;
								foreach (@{ $Regions[$r]{"Prof_Region_Ponderee_Hommes"} })
									{ $sqtotal += ($Regions[$r]{"Norm_Prof_Ponderees_Hommes"}-$_)**2; }
								$Regions[$r]{"Std_Prof_Ponderees_Hommes"} = ($sqtotal / (scalar@{ $Regions[$r]{"Prof_Region_Ponderee_Hommes"} }-1))**0.5;
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
								if ($Regions[$r]{"Norm_Prof_Ponderees_Femmes"}) {
									if ($Regions[$r]{"Chrom"} !~ /^(chrY|Y)$/) {
										if ($norm eq "std") {
											$Regions[$r]{$file}{"Moyenne_Inter"} = ($Regions[$r]{$file}{"Prof_Region_Ponderee"}-$Regions[$r]{"Norm_Prof_Ponderees_Femmes"})/$Regions[$r]{"Std_Prof_Ponderees_Femmes"};
										} else {
											$Regions[$r]{$file}{"Moyenne_Inter"} = $Regions[$r]{$file}{"Prof_Region_Ponderee"}/$Regions[$r]{"Norm_Prof_Ponderees_Femmes"};
										}
										$prof_Moyenne_Inter = $Regions[$r]{$file}{"Moyenne_Inter"};
										print SORTIE $Regions[$r]{$file}{"Moyenne_Inter"}."\t";
									} else {
										# En faisant la distinction entre chromosomes X et Y, les femmes ne sont plus considérées pour l'appel de CNV des régions du dernier
										print SORTIE "NA\t";
									}
								} else {
									$Regions[$r]{"Appel"} = "Absence_Donnees";
									print SORTIE "NA\t";
								}

							} elsif ($Patients{$file}{"Sexe"} eq "H") {
								if ($Regions[$r]{"Norm_Prof_Ponderees_Hommes"}) {
									if ($norm eq "std") {
										$Regions[$r]{$file}{"Moyenne_Inter"} = ($Regions[$r]{$file}{"Prof_Region_Ponderee"}-$Regions[$r]{"Norm_Prof_Ponderees_Hommes"})/$Regions[$r]{"Std_Prof_Ponderees_Hommes"};
									} else {
										$Regions[$r]{$file}{"Moyenne_Inter"} = $Regions[$r]{$file}{"Prof_Region_Ponderee"}/$Regions[$r]{"Norm_Prof_Ponderees_Hommes"};
									}
									$prof_Moyenne_Inter = $Regions[$r]{$file}{"Moyenne_Inter"};
									print SORTIE $Regions[$r]{$file}{"Moyenne_Inter"}."\t";
								} else {
									$Regions[$r]{"Appel"} = "Absence_Donnees";
									print SORTIE "NA\t";
								}
							}

						} else {
							# Pour les autosomes
							if ($Regions[$r]{"Chrom"} !~ /^(chr[XY]|[XY])$/) {
								if ($Regions[$r]{"Norm_Prof_Ponderees_Autosomes"}) {
									if ($norm eq "std") {
										$Regions[$r]{$file}{"Moyenne_Inter"} = ($Regions[$r]{$file}{"Prof_Region_Ponderee"}-$Regions[$r]{"Norm_Prof_Ponderees_Autosomes"})/$Regions[$r]{"Std_Prof_Ponderees_Autosomes"};
									} else {
										$Regions[$r]{$file}{"Moyenne_Inter"} = $Regions[$r]{$file}{"Prof_Region_Ponderee"}/$Regions[$r]{"Norm_Prof_Ponderees_Autosomes"};
									}
									$prof_Moyenne_Inter = $Regions[$r]{$file}{"Moyenne_Inter"};
									print SORTIE $Regions[$r]{$file}{"Moyenne_Inter"}."\t";
								} else {
									$Regions[$r]{"Appel"} = "Absence_Donnees";
									print SORTIE "NA\t";
								}
							} else {
							# Pour les gonosomes
								if($Patients{$file}{"Sexe"} eq "F") {
									if ($Regions[$r]{"Norm_Prof_Ponderees_Femmes"}) {	
										if ($Regions[$r]{"Chrom"} =~ /^(chrX|X)$/) {
											if ($norm eq "std") {
												$Regions[$r]{$file}{"Moyenne_Inter"} = ($Regions[$r]{$file}{"Prof_Region_Ponderee"}-$Regions[$r]{"Norm_Prof_Ponderees_Femmes"})/$Regions[$r]{"Std_Prof_Ponderees_Femmes"};
											} else {
												$Regions[$r]{$file}{"Moyenne_Inter"} = $Regions[$r]{$file}{"Prof_Region_Ponderee"}/$Regions[$r]{"Norm_Prof_Ponderees_Femmes"};
											}
											$prof_Moyenne_Inter = $Regions[$r]{$file}{"Moyenne_Inter"};
											print SORTIE $Regions[$r]{$file}{"Moyenne_Inter"}."\t";
										} else {
											# En faisant la distinction entre chromosomes X et Y, les femmes ne sont plus considérées pour l'appel de CNV des régions du dernier
											print SORTIE "NA\t";
										}
									} else {
										$Regions[$r]{"Appel"} = "Absence_Donnees";
										print SORTIE "NA\t";
									}
								} elsif ($Patients{$file}{"Sexe"} eq "H") {
									if ($Regions[$r]{"Norm_Prof_Ponderees_Hommes"}) {
										if ($norm eq "std") {
											$Regions[$r]{$file}{"Moyenne_Inter"} = ($Regions[$r]{$file}{"Prof_Region_Ponderee"}-$Regions[$r]{"Norm_Prof_Ponderees_Hommes"})/$Regions[$r]{"Std_Prof_Ponderees_Hommes"};
										} else {
											$Regions[$r]{$file}{"Moyenne_Inter"} = $Regions[$r]{$file}{"Prof_Region_Ponderee"}/$Regions[$r]{"Norm_Prof_Ponderees_Hommes"};
										}
										$prof_Moyenne_Inter = $Regions[$r]{$file}{"Moyenne_Inter"};
										print SORTIE $Regions[$r]{$file}{"Moyenne_Inter"}."\t";
									} else {
										$Regions[$r]{"Appel"} = "Absence_Donnees";
										print SORTIE "NA\t";
									}
								}
							}
						}
			
						#print $Regions[$r]{$patient}{"Moyenne_Inter"}."\t";
						if ($prof_Moyenne_Inter < $seuil_deletion) {
							$nb_deletions++;
							$nb_evts++;
							$CNV{$file}++;
							$Results{$file}{$r} = "DEL\t".sprintf("%.3f",$Regions[$r]{$file}{"Moyenne_Inter"});
						} elsif ($prof_Moyenne_Inter > $seuil_duplication) {
							$nb_duplications++;
							$nb_evts++;
							$CNV{$file}++;
							$Results{$file}{$r} = "DUP\t".sprintf("%.3f",$Regions[$r]{$file}{"Moyenne_Inter"});
						}

					} else {
						print SORTIE "NA\t";
					}

				}

				my$recurrent="";
				if(@{$Regions[$r]{"Prof_Region_Ponderee_Femmes"}} && @{$Regions[$r]{"Prof_Region_Ponderee_Hommes"}} && (($nb_evts/(scalar@{$Regions[$r]{"Prof_Region_Ponderee_Femmes"}} + scalar@{$Regions[$r]{"Prof_Region_Ponderee_Hommes"}})) > $seuil_region) && ($nb_parcours > 0)) { $recurrent = 1; }
				elsif(@{$Regions[$r]{"Prof_Region_Ponderee_Femmes"}} && !@{$Regions[$r]{"Prof_Region_Ponderee_Hommes"}} && (($nb_evts/scalar@{$Regions[$r]{"Prof_Region_Ponderee_Femmes"}}) > $seuil_region) && ($nb_parcours > 0)) { $recurrent = 1; }
				elsif(!@{$Regions[$r]{"Prof_Region_Ponderee_Femmes"}} && @{$Regions[$r]{"Prof_Region_Ponderee_Hommes"}} && (($nb_evts/scalar@{$Regions[$r]{"Prof_Region_Ponderee_Hommes"}}) > $seuil_region) && ($nb_parcours > 0)) { $recurrent = 1; }
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

				@{ $Regions[$r]{"Prof_Region_Ponderee"} }=();
				# Pour chaque patient, un premier parcours permet la pondération intra-patient (en divisant la profondeur de la region par la profondeur totale recalculée du patient)
				foreach my$file (@Files) {
					if(! $Patients{$file}{"ecarte"}) {
						$Regions[$r]{$file}{"Prof_Region_Ponderee"} = $Regions[$r]{$file}{"Mean"} / $Patients{$file}{"Ref_Profondeur_Patient"};
						if ( ($Patients{$file}{"Sexe"} eq "H") && ($Regions[$r]{"Chrom"} =~ /^(chrX|X)$/) )
							{ push(@{ $Regions[$r]{"Prof_Region_Ponderee"} }, ($Regions[$r]{$file}{"Prof_Region_Ponderee"}*2)); }
						elsif ( ($Patients{$file}{"Sexe"} eq "F") && ($Regions[$r]{"Chrom"} =~ /^(chrY|Y)$/) )
							{ next ; }
						else { push(@{ $Regions[$r]{"Prof_Region_Ponderee"} }, $Regions[$r]{$file}{"Prof_Region_Ponderee"}); }
						}
					}
				# Nous calculons pour chaque region la mediane/moyenne de profondeur (ponderee pour chaque patient par rapport aux autres regions)
				if(@{ $Regions[$r]{"Prof_Region_Ponderee"} }) {
					$Regions[$r]{"Norm_Prof_Ponderees"} = norm($norm,\@{ $Regions[$r]{"Prof_Region_Ponderee"} });
					if ($norm eq "std") {
						my$sqtotal = 0;
						foreach (@{ $Regions[$r]{"Prof_Region_Ponderee"} })
							{ $sqtotal += ($Regions[$r]{"Norm_Prof_Ponderees"}-$_)**2; }
						$Regions[$r]{"Std_Prof_Ponderees"} = ($sqtotal / (scalar@{ $Regions[$r]{"Prof_Region_Ponderee"} }-1))**0.5;
					}
				}

				# Un second parcours permet l'appel de CNV (en comparant la profondeur pour un patient par la profondeur moyenne de la region)
				foreach my$file (@Files) {

					if ($norm eq "std") { $prof_Moyenne_Inter = 0; }
					else  { $prof_Moyenne_Inter = 1; }

					if(! $Patients{$file}{"ecarte"}) {

						if ($Regions[$r]{"Norm_Prof_Ponderees"}) {
							if($Patients{$file}{"Sexe"} eq "F") {
								if ($Regions[$r]{"Chrom"} !~ /^(chrY|Y)$/) {
									if ($norm eq "std") {
										$Regions[$r]{$file}{"Moyenne_Inter"} = ($Regions[$r]{$file}{"Prof_Region_Ponderee"}-$Regions[$r]{"Norm_Prof_Ponderees"})/$Regions[$r]{"Std_Prof_Ponderees"};
									} else {
									$Regions[$r]{$file}{"Moyenne_Inter"} = $Regions[$r]{$file}{"Prof_Region_Ponderee"}/$Regions[$r]{"Norm_Prof_Ponderees"};
									}
									$prof_Moyenne_Inter = $Regions[$r]{$file}{"Moyenne_Inter"};
									print SORTIE $Regions[$r]{$file}{"Moyenne_Inter"}."\t";
								} else {
									#les femmes ne sont pas considérées pour l'appel de CNV des régions du Y
									print SORTIE "NA\t";
								}
							} else {
								#*2 for chrX
								if ($Regions[$r]{"Chrom"} =~ /^(chrX|X)$/) {
									if ($norm eq "std") {
										$Regions[$r]{$file}{"Moyenne_Inter"} = ($Regions[$r]{$file}{"Prof_Region_Ponderee"}*2-$Regions[$r]{"Norm_Prof_Ponderees"})/$Regions[$r]{"Std_Prof_Ponderees"};
									} else {
									$Regions[$r]{$file}{"Moyenne_Inter"} = $Regions[$r]{$file}{"Prof_Region_Ponderee"}*2/$Regions[$r]{"Norm_Prof_Ponderees"};
									}
								} else {
									$Regions[$r]{$file}{"Moyenne_Inter"} = $Regions[$r]{$file}{"Prof_Region_Ponderee"}/$Regions[$r]{"Norm_Prof_Ponderees"};
									if ($norm eq "std") {
										$Regions[$r]{$file}{"Moyenne_Inter"} = ($Regions[$r]{$file}{"Prof_Region_Ponderee"}-$Regions[$r]{"Norm_Prof_Ponderees"})/$Regions[$r]{"Std_Prof_Ponderees"};
									} else {
									$Regions[$r]{$file}{"Moyenne_Inter"} = $Regions[$r]{$file}{"Prof_Region_Ponderee"}/$Regions[$r]{"Norm_Prof_Ponderees"};
									}
								}
								$prof_Moyenne_Inter = $Regions[$r]{$file}{"Moyenne_Inter"};
								print SORTIE $Regions[$r]{$file}{"Moyenne_Inter"}."\t";
							}
							if ($prof_Moyenne_Inter < $seuil_deletion) {
								$nb_deletions++;
								$nb_evts++;
								$CNV{$file}++;
								$Results{$file}{$r} = "DEL\t".sprintf("%.3f",$Regions[$r]{$file}{"Moyenne_Inter"});
							} elsif ($prof_Moyenne_Inter > $seuil_duplication) {
								$nb_duplications++;
								$nb_evts++;
								$CNV{$file}++;
								$Results{$file}{$r} = "DUP\t".sprintf("%.3f",$Regions[$r]{$file}{"Moyenne_Inter"});
							}
						} else {
							$Regions[$r]{"Appel"} = "Absence_Donnees";
							#$prof_Moyenne_Inter = 1;
							print SORTIE "NA\t";
						}
					} else {
						print SORTIE "NA\t";
					}
				}

				if(@{ $Regions[$r]{"Prof_Region_Ponderee"} } && (($nb_evts/scalar@{ $Regions[$r]{"Prof_Region_Ponderee"} }) > $seuil_region) && ($nb_parcours > 0)) {
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

## @{ $regionOrder{$Chrom} } = [ regions sorted by pos ]
my(%uniqReg,%RegInArray,%RegionOrder,%regionIndice);
for my$r (0..$#Regions) {
	unless (exists $uniqReg{ $Regions[$r]{"Chrom"}."-".$Regions[$r]{"Start"}."-".$Regions[$r]{"End"} }) {
		push (@{ $RegInArray{ $Regions[$r]{"Chrom"} } }, $r);
		$uniqReg{ $Regions[$r]{"Chrom"}."-".$Regions[$r]{"Start"}."-".$Regions[$r]{"End"} } = 1;
		}
	my@tab = split(/\t/, $Regions[$r]{"allLine"});
	if ($tab[3]) {
		my@tab2 = split(/:|,/, $tab[3]);
		my$label = "";
		if (length($tab2[0]) > 20) {
			my@tab3 =  split(//, $tab2[0]);
			for my $i (0..20) { $label .= $tab3[$i]; }
			}
		else { $label = $tab2[0]; }
		$Regions[$r]{"label"} = "$label";
		$Regions[$r]{"Gene"} = $tab2[0];
		}
	else { $Regions[$r]{"label"} = "$tab[1]-$tab[2]"; }
	}
foreach my$Chrom (keys%RegInArray) {
	##sort by region ends (in case several regions with same start) then by starts
	@{ $RegionOrder{$Chrom} } = sort{$Regions[$a]{"End"}<=>$Regions[$b]{"End"}}@{ $RegInArray{$Chrom} }; 
	@{ $RegionOrder{$Chrom} } = sort{$Regions[$a]{"Start"}<=>$Regions[$b]{"Start"}}@{ $RegInArray{$Chrom} }; 
	for (my$i=0;$i<scalar@{ $RegionOrder{$Chrom} };$i++) { $regionIndice{ $RegionOrder{$Chrom}[$i] } = $i; }
	}

##print CNV foreach sample
#$results{$patient}{$region} = "DUP (Moyenne_Inter) (Borne_3P-Borne_5P bp)";
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
my(%Result2,%Result3);
foreach my$file (keys%Results) {
	open (CNV1,">$outdir/$sampleName{$file}/CNV_$sampleName{$file}.summary.txt") or die "Pb lors de l'ecriture du fichier sortie $!\n";
	print CNV1 "#Region\tCNV\tN_intervals\n";
	open (CNV2,">$outdir/$sampleName{$file}/CNV_$sampleName{$file}.allIntervals.txt") or die "Pb lors de l'ecriture du fichier sortie $!\n";
	print CNV2 "#Chrom\tStart\tEnd\tInfo\tinterval_Order\tCNV_type\tratio_to_$norm\tinterval_Length\n";
	foreach my$Chrom (@ChromOrder) {
		if (exists $orderedCNV{$file}{$Chrom}) {
			my$r=0;		##index in @{ $orderedCNV{$file}{$Chrom} }
			while ($r < scalar@{ $orderedCNV{$file}{$Chrom} } ) {
				my$ok=1; my$i=0; my$cnvOK=0;
				my@nextReg=($orderedCNV{$file}{$Chrom}[$r]);
				my@tab = split(/\s/,$Results{$file}{$orderedCNV{$file}{$Chrom}[$r]});
				while ($ok) {
					$i++;
					if (exists $RegionOrder{$Chrom}[($regionIndice{$orderedCNV{$file}{$Chrom}[$r]} + $i)]) {
						push(@nextReg, $RegionOrder{$Chrom}[($regionIndice{$orderedCNV{$file}{$Chrom}[$r]} + $i)]);
						if (exists $Results{$file}{$nextReg[$i]}) {
							my@tab2 = split(/\s/, $Results{$file}{$nextReg[$i]});
							if($tab2[0] eq $tab[0]) { $cnvOK = $i; }	##same CNV type
							else { $ok=0; }
							}
						else {
							unless (exists $Regions[$nextReg[$i]]{"Appel"}) { $ok=0; }
							}
						}
					else { $ok=0; }
					}
				if ($cnvOK >= ($minCNV-1)) {
					if ($cnvOK == 0) {
						$Result2{$file}{$nextReg[0]} = $Results{$file}{$nextReg[0]};
						$Result3{$file}{$Chrom}{$Regions[$nextReg[0]]{"Start"}}{$Regions[$nextReg[0]]{"End"}} = $Results{$file}{$nextReg[0]};
						if (exists $regionIndice{$nextReg[0]}) {
							print CNV1 $Chrom.":".$Regions[$nextReg[0]]{"Start"}."-".$Regions[$nextReg[0]]{"End"}."\t".$tab[0]."\t1\n";
							print CNV2 $Chrom."\t".$Regions[$nextReg[0]]{"Start"}."\t".$Regions[$nextReg[0]]{"End"}."\t".$Regions[$nextReg[0]]{"label"}."\t".($regionIndice{$nextReg[0]}+1)."\t".$Results{$file}{$nextReg[0]}."\t".($Regions[$nextReg[0]]{"End"}-$Regions[$nextReg[0]]{"Start"}+1)." bp\n\n";
							}
						}
					else {
						print CNV1 $Chrom.":".$Regions[$nextReg[0]]{"Start"}."-".$Regions[$nextReg[$cnvOK]]{"End"}."\t".$tab[0]."\t".($cnvOK+1)."\n";
						for (my$j=0;$j<=$cnvOK;$j++) {
							if (exists $Results{$file}{$nextReg[$j]}) {
								$Result2{$file}{$nextReg[$j]} = $Results{$file}{$nextReg[$j]};
								$Result3{$file}{$Chrom}{$Regions[$nextReg[$j]]{"Start"}}{$Regions[$nextReg[$j]]{"End"}} = $Results{$file}{$nextReg[$j]};
								if (exists $regionIndice{$nextReg[$j]})
									{ print CNV2 $Chrom."\t".$Regions[$nextReg[$j]]{"Start"}."\t".$Regions[$nextReg[$j]]{"End"}."\t".$Regions[$nextReg[$j]]{"label"}."\t".($regionIndice{$nextReg[$j]}+1)."\t".$Results{$file}{$nextReg[$j]}."\t".($Regions[$nextReg[$j]]{"End"}-$Regions[$nextReg[$j]]{"Start"}+1)." bp\n"; }
								}
							else {
								$Result2{$file}{$nextReg[$j]} = "NA";
								$Result3{$file}{$Chrom}{$Regions[$nextReg[$j]]{"Start"}}{$Regions[$nextReg[$j]]{"End"}} = "NA";
								if (exists $regionIndice{$nextReg[$j]}) {
									print CNV2 $Chrom."\t".$Regions[$nextReg[$j]]{"Start"}."\t".$Regions[$nextReg[$j]]{"End"}."\t".$Regions[$nextReg[$j]]{"label"}."\t".($regionIndice{$nextReg[$j]}+1)."\tNA\t";
									if (exists $Regions[$nextReg[$j]]{$file}{"Moyenne_Inter"})
										{ print CNV2 sprintf("%.3f",$Regions[$nextReg[$j]]{$file}{"Moyenne_Inter"}); }
									else { print CNV2 "na"; }
									print CNV2 "\t".($Regions[$nextReg[$j]]{"End"}-$Regions[$nextReg[$j]]{"Start"})." bp"."\n";
									}
								}
							}
						print CNV2 "\n";
						}
					}
				if ($i > 1) { $r += ($i-1); }
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
		graphChr1($nGraf,"$outdir/$sampleName{$file}",$norm,$seuil_deletion,$seuil_duplication,$file,$maxDepthGraph,\@Files,\%sampleName,\@Regions,\@ChromOrder,\%RegionOrder,\%Patients,\%Result2);
		}
	}



return(\%Patients,\%Result3);

}


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
my$graphByChr = $CNV_opt{"chromGraph"};

my$maxDepthGraph = 10;

my %Sexe;
if ($fichier_sexe) {
	open(PATIENTS, "<$fichier_sexe") or die "Fichier $fichier_sexe d'attribution des sexes inexistant ou impossible a lire\n";
	foreach my $ligne (<PATIENTS>) {
		$ligne =~ m/(\S+)\s+(\S+)/;
		$Sexe{$1} = $2;
	}
	close(PATIENTS);
	foreach my $sample (keys%Sexe) {
		if ($Sexe{$sample} =~ /^m|^h/i) { $Sexe{$sample} = "H"; }
		else { $Sexe{$sample} = "F"; }
	}
}


# Pour chaque region de la sortie DecoVa, on recupere les infos de couverture
my @idxP;	#idx of Patients
my %idxV;	#idx of Values foreach patient
my $nombre_patients;
my %Patients;
my $nombre_regions = 0;
my $autosom_Regions = 0;
my $gonosom_Regions = 0;
my @Regions;
my @ChrOrder; my%allChr;
my $sexe_courant;

# PREMIER PARCOURS #
# PERMET DE CALCULER LA PROFONDEUR TOTALE POUR CHAQUE PATIENT, EN TENANT UNIQUEMENT COMPTE DES REGIONS CORRECTEMENT SEQUENCEES
my$l=0; #N° ligne
unless (-f "$fichier_cov") { die "Fichier $fichier_cov inexistant\n"; }
open(DECOVA, "<$fichier_cov") or die "Fichier $fichier_cov inexistant ou impossible a lire\n";
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
				$Patients{$p}{"Profondeur_Autosomes_Patient"} = 0;
				$Patients{$p}{"Profondeur_X_Patient"} = 0;
				$Patients{$p}{"Profondeur_Y_Patient"} = 0;
				$Patients{$p}{"Ref_Profondeur_Patient"} = 0;
				$p++;
			}
		$nombre_patients = scalar@idxP;
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
		$Regions[$nombre_regions]{"Chromosome"} = $INFOS[0];
		$Regions[$nombre_regions]{"Borne_5P"} = $INFOS[1];
		$Regions[$nombre_regions]{"Borne_3P"} = $INFOS[2];
		if ($idxP[0]>3) { $Regions[$nombre_regions]{"Gene"} = $INFOS[3]; }
		else { $Regions[$nombre_regions]{"Gene"} = "NA"; }
		for (my$i=0;$i<$idxP[0];$i++)
			{ $Regions[$nombre_regions]{"allInfos"} .= $INFOS[$i]."\t"; }
		chop $Regions[$nombre_regions]{"allInfos"};
		# to get chrom order from bed
		unless (exists $allChr{$INFOS[0]}) {
			push(@ChrOrder, $INFOS[0]);
			$allChr{$INFOS[0]} = 1;
		}

		# On parcourt la fin de la region (= infos sur les patients uniquement) une premiere fois pour savoir si la region est a conserver
		my $region_a_conserver = 0;
		if ($seuil_cov) {
			for (my$p=0;$p<scalar@idxP;$p++) {
				# Lecture de chaque colonne "%>=.X" : si aucun patient n'a 100% de bases séquencées à .X pour la région, elle est supprimée de l'analyse
				my$cov = $INFOS[$idxV{$p}{"cov"}];
				$cov =~ s/%$//; $cov /= 100;
				if ($cov > $seuil_cov) {
					$region_a_conserver = 1; last;
				}
			}
		} else { $region_a_conserver = 1; }

		# Si la region est mal couverte, on l'imprime dans un fichier POUBELLE
		if ($region_a_conserver != 1) {
			$Regions[$nombre_regions]{"Appel"} = "Couverture_Faible";

		} else {
		# Sinon, on parcourt la region une seconde fois pour recuperer les informations désirées pour chacun des patients (et donc la profondeur moyenne pour la region)
			for (my $patient = 0 ; $patient < $nombre_patients ; $patient++) {
				$Regions[$nombre_regions]{$patient}{"Prof_Moy_Region_Patient"} = $INFOS[$idxV{$patient}{"mean"}];
				if ($Ref eq "mean") {
					$Regions[$nombre_regions]{$patient}{"Prof_Tot_Region_Patient"} = $INFOS[$idxV{$patient}{"mean"}];
				} else {
					$Regions[$nombre_regions]{$patient}{"Prof_Tot_Region_Patient"} = $INFOS[$idxV{$patient}{"tot"}];
				}
				if ($INFOS[0] !~ m/^(chr[XY]|[XY])$/) {
					$autosom_Regions++;
					$Patients{$patient}{"Profondeur_Autosomes_Patient"} += $Regions[$nombre_regions]{$patient}{"Prof_Tot_Region_Patient"};
				} else {
					if ($INFOS[0] =~ m/^(chrX|X)$/) {
						$gonosom_Regions++;
						$Patients{$patient}{"Profondeur_X_Patient"} += $Regions[$nombre_regions]{$patient}{"Prof_Tot_Region_Patient"};
					} else {
						$Patients{$patient}{"Profondeur_Y_Patient"} += $Regions[$nombre_regions]{$patient}{"Prof_Tot_Region_Patient"};

					} 
				}
			}
		}
		$nombre_regions++;
	}
}
close(DECOVA);

if ($seuil_region < (1/$nombre_patients)) { $seuil_region = (1/$nombre_patients); }
if ($seuil_patient < (1/$nombre_regions)) { $seuil_patient = (1/$nombre_regions); }

for (my $patient = 0 ; $patient < $nombre_patients ; $patient++) {
	print $Patients{$patient}{"ID"}.":\n";
	if (! defined $Patients{$patient}{"Sexe"}) {
		if ($Patients{$patient}{"Profondeur_X_Patient"} && $gonosom_Regions && $Patients{$patient}{"Profondeur_Autosomes_Patient"} && $autosom_Regions) {
			print "\tautoZ: ".($Patients{$patient}{"Profondeur_Autosomes_Patient"}/$autosom_Regions)."\n\tchrX: ".($Patients{$patient}{"Profondeur_X_Patient"}/$gonosom_Regions)."\n";
			if ( ($Patients{$patient}{"Profondeur_X_Patient"}/$gonosom_Regions) > (1.2*$Patients{$patient}{"Profondeur_Autosomes_Patient"}/$autosom_Regions) ) { print "\t-> sexe ambigu\n"; $Patients{$patient}{"Sexe"} = "F"; } 
			elsif ( (($Patients{$patient}{"Profondeur_X_Patient"}/$gonosom_Regions) <= (1.2*$Patients{$patient}{"Profondeur_Autosomes_Patient"}/$autosom_Regions)) && (($Patients{$patient}{"Profondeur_X_Patient"}/$gonosom_Regions) >= (0.8*$Patients{$patient}{"Profondeur_Autosomes_Patient"}/$autosom_Regions)) ) { $Patients{$patient}{"Sexe"} = "F"; } 
			elsif ( (($Patients{$patient}{"Profondeur_X_Patient"}/$gonosom_Regions) < (0.8*$Patients{$patient}{"Profondeur_Autosomes_Patient"}/$autosom_Regions)) && (($Patients{$patient}{"Profondeur_X_Patient"}/$gonosom_Regions) > (1.2*0.5*$Patients{$patient}{"Profondeur_Autosomes_Patient"}/$autosom_Regions)) ) { print "\t-> sexe ambigu\n"; $Patients{$patient}{"Sexe"} = "F"; } 
			elsif ( (($Patients{$patient}{"Profondeur_X_Patient"}/$gonosom_Regions) <= (1.2*0.5*$Patients{$patient}{"Profondeur_Autosomes_Patient"}/$autosom_Regions)) && (($Patients{$patient}{"Profondeur_X_Patient"}/$gonosom_Regions) >= (0.8*0.5*$Patients{$patient}{"Profondeur_Autosomes_Patient"}/$autosom_Regions)) ) { $Patients{$patient}{"Sexe"} = "H"; }
			else { print "\t-> sexe ambigu\n"; $Patients{$patient}{"Sexe"} = "H"; }
			
		} else { $Patients{$patient}{"Sexe"} = "F"; }
	}
	print "\t-> $Patients{$patient}{Sexe}\n";

	if ($RefByGender) {
		if ($Patients{$patient}{"Sexe"} eq "H") {
			$Patients{$patient}{"Ref_Profondeur_Patient"} = ($Patients{$patient}{"Profondeur_X_Patient"} * 2) + $Patients{$patient}{"Profondeur_Autosomes_Patient"}; 
		} else {
			$Patients{$patient}{"Ref_Profondeur_Patient"} = $Patients{$patient}{"Profondeur_X_Patient"} + $Patients{$patient}{"Profondeur_Autosomes_Patient"};
		}
	} else {
		$Patients{$patient}{"Ref_Profondeur_Patient"} = $Patients{$patient}{"Profondeur_X_Patient"} + $Patients{$patient}{"Profondeur_Y_Patient"} + $Patients{$patient}{"Profondeur_Autosomes_Patient"};
	}
	unless ($Patients{$patient}{"Ref_Profondeur_Patient"}) { $Patients{$patient}{"ecarte"} = 1; }
}


my $nb_parcours = 0;
my $continuer = 1;
my %results;

##iterations
while ($continuer == 1 || $nb_parcours <= 1) {

	my $iteration = $nb_parcours + 1;
	print "Iteration numero \: ".$iteration."\n";

	undef %results;
	my %CNV;
	my $nb_regions_conservees = 0;
	my $regions_ecartees = 0;
	$continuer = 0;
	my $sortie = "$outdir/CNV_Iteration_".$iteration."\.tab";
	open(SORTIE,">", $sortie) or die("Pb lors de l'ecriture du fichier sortie $!\n");
	my $log = "$outdir/Logs_Iteration_".$iteration."\.txt";
	open(LOG,">", $log) or die("Pb lors de l'ecriture du fichier log $!\n");

	print SORTIE "Chromosome"."\t";
	print SORTIE "Borne_5P"."\t";
	print SORTIE "Borne_3P"."\t";
	print SORTIE "Gene"."\t";
	print SORTIE "Numero_Region"."\t";

	for (my $patient = 0 ; $patient < $nombre_patients ; $patient++) {
		printf SORTIE $Patients{$patient}{"ID"}."\t";

	}

	print SORTIE "Statut_Region\n";


	# PREMIER PARCOURS BIS #
	# RECALCUL DE LA PROFONDEUR TOTALE POUR LA PONDERATION INTRA
	# CE (RE)CALCUL A LIEU SI DES REGIONS ETAIENT "MOCHES" APRES L'APPEL DE CNV
	if ($nb_parcours > 0) {

		# REINITIALISATION DE LA PROFONDEUR TOTALE PAR PATIENT
		for (my $patient = 0 ; $patient < $nombre_patients ; $patient++) {
			$Patients{$patient}{"Profondeur_Autosomes_Patient"} = 0;
			$Patients{$patient}{"Profondeur_Gonosomes_Patient"} = 0;
		}

		for (my $region = 0 ; $region < $nombre_regions ; $region++) {
			if(!(defined($Regions[$region]{"Appel"}))) {

				for (my $patient = 0 ; $patient < $nombre_patients ; $patient++) {
					if ($Regions[$region]{"Chromosome"} !~ m/^chr[XY]$|[XY]$/) {
						$Patients{$patient}{"Profondeur_Autosomes_Patient"} += $Regions[$region]{$patient}{"Prof_Tot_Region_Patient"};
					} elsif ($Regions[$region]{"Chromosome"} =~ m/^chrX$|^X$/) {
						$Patients{$patient}{"Profondeur_Gonosomes_Patient"} += $Regions[$region]{$patient}{"Prof_Tot_Region_Patient"};
					} else {
						if ((!$RefByGender) && ($Patients{$patient}{"Sexe"} eq "H")) {
							$Patients{$patient}{"Profondeur_Gonosomes_Patient"} += $Regions[$nombre_regions]{$patient}{"Prof_Tot_Region_Patient"};
						}
					}
				}
			}
		}


		for (my $patient = 0 ; $patient < $nombre_patients ; $patient++) {
			if ($RefByGender && ($Patients{$patient}{"Sexe"} eq "H")) {
				$Patients{$patient}{"Ref_Profondeur_Patient"} = ($Patients{$patient}{"Profondeur_Gonosomes_Patient"} * 2) + $Patients{$patient}{"Profondeur_Autosomes_Patient"};
			} else {
				$Patients{$patient}{"Ref_Profondeur_Patient"} = $Patients{$patient}{"Profondeur_Gonosomes_Patient"} + $Patients{$patient}{"Profondeur_Autosomes_Patient"};
			}
		}
	}


	# SECOND PARCOURS DES REGIONS #
	# PERMET DE PONDERER LA PROFONDEUR PAR LES AUTRES REGIONS (INTRA) ET ENTRE LES PATIENTS (INTER)
	for (my $region = 0 ; $region < $nombre_regions ; $region++) {

		print SORTIE $Regions[$region]{"Chromosome"}."\t";
		print SORTIE $Regions[$region]{"Borne_5P"}."\t";
		print SORTIE $Regions[$region]{"Borne_3P"}."\t";
		print SORTIE $Regions[$region]{"Gene"}."\t";
		print SORTIE $region."\t";

		# SI LA REGION N'EST PAS "MOCHE"
		if (!defined($Regions[$region]{"Appel"})) {

			$nb_regions_conservees++;

			$Regions[$region]{"Profondeur_Totale_Region"} = 0;
			$Regions[$region]{"Somme_Profondeur_Region_Ponderee_Autosomes"} = 0;
			$Regions[$region]{"Somme_Profondeur_Region_Ponderee_Femmes"} = 0;
			$Regions[$region]{"Somme_Profondeur_Region_Ponderee_Hommes"} = 0;
			my $nb_deletions = 0;
			my $nb_duplications = 0;
			my $nb_evts = 0;
			my $prof_Moyenne_Inter = "";


			if ($normByGender) {

				# Pour chaque patient, un premier parcours permet la pondération intra-patient (en divisant la profondeur de la region par la profondeur totale du patient)
				@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Femmes"} } = ();
				@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Hommes"} } = ();

				for (my $patient = 0 ; $patient < $nombre_patients ; $patient++) {

					if($Patients{$patient}{"Sexe"} eq "F" && !($Patients{$patient}{"ecarte"})) {

						# Le controle permet de diviser uniquement par le nombre de patients restant apres l filtre, et non par l'ensemble des patients
						# Si la region concerne un gonosome et le patient est une femme, la référence est celle qui correspond aux autosomes + au chromosome X
						$Regions[$region]{$patient}{"Profondeur_Region_Ponderee"} = $Regions[$region]{$patient}{"Prof_Moy_Region_Patient"}/$Patients{$patient}{"Ref_Profondeur_Patient"};
						# On sauvegarde la profondeur ponderee intra patient, pour determiner ensuite la moyenne de ces profondeurs ponderees
						push( @{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Femmes"} }, $Regions[$region]{$patient}{"Profondeur_Region_Ponderee"} );


					} elsif ($Patients{$patient}{"Sexe"} eq "H" && !($Patients{$patient}{"ecarte"})) {

						# Le controle permet de divisier uniquement par le nombre de patients restant apres l filtre, et non par l'ensemble des patients
						# Si la region concerne un gonosome et le patient est un homme, la référence est celle qui correspond uniquement au chromosome X
						$Regions[$region]{$patient}{"Profondeur_Region_Ponderee"} = $Regions[$region]{$patient}{"Prof_Moy_Region_Patient"}/$Patients{$patient}{"Ref_Profondeur_Patient"};
						# On sauvegarde la profondeur ponderee intra patient, pour determiner ensuite la moyenne de ces profondeurs ponderees
						push( @{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Hommes"} }, $Regions[$region]{$patient}{"Profondeur_Region_Ponderee"} );

					}

				}

				# Nous calculons pour chaque region et chaque sexe la moyenne de profondeur (ponderee pour chaque patient par rapport aux autres regions)
				if ($normByGender eq "all") {
					if(@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Femmes"} }) {
						 $Regions[$region]{"Norm_Profondeurs_Ponderees_Femmes"} = norm($norm,\@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Femmes"} });
						if ($norm eq "std") {
							my$sqtotal = 0;
							foreach (@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Femmes"} })
								{ $sqtotal += ($Regions[$region]{"Norm_Profondeurs_Ponderees_Femmes"}-$_)**2; }
							$Regions[$region]{"Std_Profondeurs_Ponderees_Femmes"} = ($sqtotal / (scalar@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Femmes"} }-1))**0.5;
						}
					}
					if(@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Hommes"} }) {
						$Regions[$region]{"Norm_Profondeurs_Ponderees_Hommes"} = norm($norm,\@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Hommes"} });
						if ($norm eq "std") {
							my$sqtotal = 0;
							foreach (@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Hommes"} })
								{ $sqtotal += ($Regions[$region]{"Norm_Profondeurs_Ponderees_Hommes"}-$_)**2; }
							$Regions[$region]{"Std_Profondeurs_Ponderees_Hommes"} = ($sqtotal / (scalar@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Hommes"} }-1))**0.5;
						}
					}
				} else {
					if ($Regions[$region]{"Chromosome"} !~ m/^chr[XY]$|^[XY]$/) {
						my@Autosomes=();
						if(@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Femmes"} }) {
							push(@Autosomes, @{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Femmes"} });
						}
						if(@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Hommes"} }) {
							push(@Autosomes, @{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Hommes"} });
						}
						if (@Autosomes) {
							$Regions[$region]{"Norm_Profondeurs_Ponderees_Autosomes"} = norm($norm,\@Autosomes);
							if ($norm eq "std") {
								my$sqtotal = 0;
								foreach (@Autosomes)
									{ $sqtotal += ($Regions[$region]{"Norm_Profondeurs_Ponderees_Autosomes"}-$_)**2; }
								$Regions[$region]{"Std_Profondeurs_Ponderees_Autosomes"} = ($sqtotal / (scalar@Autosomes-1))**0.5;
							}
						}
					} else {
						if(@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Femmes"} }) {
							 $Regions[$region]{"Norm_Profondeurs_Ponderees_Femmes"} = norm($norm,\@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Femmes"} });
							if ($norm eq "std") {
								my$sqtotal = 0;
								foreach (@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Femmes"} })
									{ $sqtotal += ($Regions[$region]{"Norm_Profondeurs_Ponderees_Femmes"}-$_)**2; }
								$Regions[$region]{"Std_Profondeurs_Ponderees_Femmes"} = ($sqtotal / (scalar@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Femmes"} }-1))**0.5;
							}
						}
						if(@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Hommes"} }) {
							$Regions[$region]{"Norm_Profondeurs_Ponderees_Hommes"} = norm($norm,\@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Hommes"} });
							if ($norm eq "std") {
								my$sqtotal = 0;
								foreach (@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Hommes"} })
									{ $sqtotal += ($Regions[$region]{"Norm_Profondeurs_Ponderees_Hommes"}-$_)**2; }
								$Regions[$region]{"Std_Profondeurs_Ponderees_Hommes"} = ($sqtotal / (scalar@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Hommes"} }-1))**0.5;
							}
						}
					}
				}
				
				# Un second parcours permet l'appel de CNV (en comparant la profondeur pour un patient par la profondeur moyenne de la region)
				for (my $patient = 0 ; $patient < $nombre_patients ; $patient++) {

					if ($norm eq "std") { $prof_Moyenne_Inter = 0; }
					else  { $prof_Moyenne_Inter = 1; }

					unless ($Patients{$patient}{"ecarte"}) {

						if ($normByGender eq "all") {

							if($Patients{$patient}{"Sexe"} eq "F") {
								if ($Regions[$region]{"Norm_Profondeurs_Ponderees_Femmes"}) {
									if ($Regions[$region]{"Chromosome"} !~ m/^(chrY|Y)$/) {
										if ($norm eq "std") {
											$Regions[$region]{$patient}{"Moyenne_Inter"} = ($Regions[$region]{$patient}{"Profondeur_Region_Ponderee"}-$Regions[$region]{"Norm_Profondeurs_Ponderees_Femmes"})/$Regions[$region]{"Std_Profondeurs_Ponderees_Femmes"};
										} else {
											$Regions[$region]{$patient}{"Moyenne_Inter"} = $Regions[$region]{$patient}{"Profondeur_Region_Ponderee"}/$Regions[$region]{"Norm_Profondeurs_Ponderees_Femmes"};
										}
										$prof_Moyenne_Inter = $Regions[$region]{$patient}{"Moyenne_Inter"};
										print SORTIE $Regions[$region]{$patient}{"Moyenne_Inter"}."\t";
									} else {
										# En faisant la distinction entre chromosomes X et Y, les femmes ne sont plus considérées pour l'appel de CNV des régions du dernier
										print SORTIE "NA\t";
									}
								} else {
									$Regions[$region]{"Appel"} = "Absence_Donnees";
									print SORTIE "NA\t";
								}

							} elsif ($Patients{$patient}{"Sexe"} eq "H") {
								if ($Regions[$region]{"Norm_Profondeurs_Ponderees_Hommes"}) {
									if ($norm eq "std") {
										$Regions[$region]{$patient}{"Moyenne_Inter"} = ($Regions[$region]{$patient}{"Profondeur_Region_Ponderee"}-$Regions[$region]{"Norm_Profondeurs_Ponderees_Hommes"})/$Regions[$region]{"Std_Profondeurs_Ponderees_Hommes"};
									} else {
										$Regions[$region]{$patient}{"Moyenne_Inter"} = $Regions[$region]{$patient}{"Profondeur_Region_Ponderee"}/$Regions[$region]{"Norm_Profondeurs_Ponderees_Hommes"};
									}
									$prof_Moyenne_Inter = $Regions[$region]{$patient}{"Moyenne_Inter"};
									print SORTIE $Regions[$region]{$patient}{"Moyenne_Inter"}."\t";
								} else {
									$Regions[$region]{"Appel"} = "Absence_Donnees";
									print SORTIE "NA\t";
								}
							}

						} else {
							# Pour les autosomes
							if ($Regions[$region]{"Chromosome"} !~ m/^chr[XY]$|^[XY]$/) {
								if ($Regions[$region]{"Norm_Profondeurs_Ponderees_Autosomes"}) {
									if ($norm eq "std") {
										$Regions[$region]{$patient}{"Moyenne_Inter"} = ($Regions[$region]{$patient}{"Profondeur_Region_Ponderee"} - $Regions[$region]{"Norm_Profondeurs_Ponderees_Autosomes"}) / $Regions[$region]{"Std_Profondeurs_Ponderees_Autosomes"};
									} else {
										$Regions[$region]{$patient}{"Moyenne_Inter"} = $Regions[$region]{$patient}{"Profondeur_Region_Ponderee"} / $Regions[$region]{"Norm_Profondeurs_Ponderees_Autosomes"};
									}
									$prof_Moyenne_Inter = $Regions[$region]{$patient}{"Moyenne_Inter"};
									print SORTIE $Regions[$region]{$patient}{"Moyenne_Inter"}."\t";
								} else {
									$Regions[$region]{"Appel"} = "Absence_Donnees";
									print SORTIE "NA\t";
								}
							} else {
								if($Patients{$patient}{"Sexe"} eq "F") {
									if ($Regions[$region]{"Norm_Profondeurs_Ponderees_Femmes"}) {	
										if ($Regions[$region]{"Chromosome"} =~ m/^(chrX|X)$/) {
											if ($norm eq "std") {
												$Regions[$region]{$patient}{"Moyenne_Inter"} = ($Regions[$region]{$patient}{"Profondeur_Region_Ponderee"}-$Regions[$region]{"Norm_Profondeurs_Ponderees_Femmes"})/$Regions[$region]{"Std_Profondeurs_Ponderees_Femmes"};
											} else {
												$Regions[$region]{$patient}{"Moyenne_Inter"} = $Regions[$region]{$patient}{"Profondeur_Region_Ponderee"}/$Regions[$region]{"Norm_Profondeurs_Ponderees_Femmes"};
											}
											$prof_Moyenne_Inter = $Regions[$region]{$patient}{"Moyenne_Inter"};
											print SORTIE $Regions[$region]{$patient}{"Moyenne_Inter"}."\t";
										} else {
											# En faisant la distinction entre chromosomes X et Y, les femmes ne sont plus considérées pour l'appel de CNV des régions du dernier
											print SORTIE "NA\t";
										}
									} else {
										$Regions[$region]{"Appel"} = "Absence_Donnees";
										print SORTIE "NA\t";
									}
								} elsif ($Patients{$patient}{"Sexe"} eq "H") {
									if ($Regions[$region]{"Norm_Profondeurs_Ponderees_Hommes"}) {
										if ($norm eq "std") {
											$Regions[$region]{$patient}{"Moyenne_Inter"} = ($Regions[$region]{$patient}{"Profondeur_Region_Ponderee"}-$Regions[$region]{"Norm_Profondeurs_Ponderees_Hommes"})/$Regions[$region]{"Std_Profondeurs_Ponderees_Hommes"};
										} else {
											$Regions[$region]{$patient}{"Moyenne_Inter"} = $Regions[$region]{$patient}{"Profondeur_Region_Ponderee"}/$Regions[$region]{"Norm_Profondeurs_Ponderees_Hommes"};
										}
										$prof_Moyenne_Inter = $Regions[$region]{$patient}{"Moyenne_Inter"};
										print SORTIE $Regions[$region]{$patient}{"Moyenne_Inter"}."\t";
									} else {
										$Regions[$region]{"Appel"} = "Absence_Donnees";
										print SORTIE "NA\t";
									}
								}
							}
						}
				
						if ($prof_Moyenne_Inter < $seuil_deletion) {
							$nb_deletions++;
							$nb_evts++;
							$CNV{$patient}++;
							$results{$patient}{$region} = "DEL\t".sprintf("%.3f",$Regions[$region]{$patient}{"Moyenne_Inter"});
						} elsif ($prof_Moyenne_Inter > $seuil_duplication) {
							$nb_duplications++;
							$nb_evts++;
							$CNV{$patient}++;
							$results{$patient}{$region} = "DUP\t".sprintf("%.3f",$Regions[$region]{$patient}{"Moyenne_Inter"});
						}

					} else {
						print SORTIE "NA\t";
					}

				}

				my$recurrent="";
				if(@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Femmes"} } && @{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Hommes"} } && (($nb_evts/(scalar@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Femmes"} } + scalar@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Hommes"} })) > $seuil_region) && ($nb_parcours > 0)) { $recurrent = 1; }
				elsif(@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Femmes"} } && !@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Hommes"} } && (($nb_evts/scalar@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Femmes"} }) > $seuil_region) && ($nb_parcours > 0)) { $recurrent = 1; }
				elsif(!@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Femmes"} } && @{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Hommes"} } && (($nb_evts/scalar@{ $Regions[$region]{"ProfondeurS_Region_Ponderee_Hommes"} }) > $seuil_region) && ($nb_parcours > 0)) { $recurrent = 1; }
				if ($recurrent) {
					print SORTIE "CNV_Recurrent";
					print LOG "Region ecartee \: ". $Regions[$region]{"Chromosome"}."\t".$Regions[$region]{"Borne_5P"}."\t".$Regions[$region]{"Borne_3P"}."\t".$Regions[$region]{"Gene"}."\n";
					$Regions[$region]{"Appel"} = "CNV_Recurrent";
					$regions_ecartees++;
					$continuer = 1;
				} else {		
					print SORTIE "OK";
				}

				print SORTIE "\n";



			} else {

				@{ $Regions[$region]{"ProfondeurS_Region_Ponderee"} }=();
				# Pour chaque patient, un premier parcours permet la pondération intra-patient (en divisant la profondeur de la region par la profondeur totale recalculée du patient)
				for (my $patient = 0 ; $patient < $nombre_patients ; $patient++) {
					unless ($Patients{$patient}{"ecarte"}) {
						$Regions[$region]{$patient}{"Profondeur_Region_Ponderee"} = $Regions[$region]{$patient}{"Prof_Moy_Region_Patient"}/$Patients{$patient}{"Ref_Profondeur_Patient"};
						if ( ($Patients{$patient}{"Sexe"} eq "H") && ($Regions[$region]{"Chromosome"} =~ m/^(chrX|X)$/) )
							{ push(@{ $Regions[$region]{"ProfondeurS_Region_Ponderee"} }, ($Regions[$region]{$patient}{"Profondeur_Region_Ponderee"}*2)); }
						elsif ( ($Patients{$patient}{"Sexe"} eq "F") && ($Regions[$region]{"Chromosome"} =~ m/^(chrY|Y)$/) )
							{ next ; }
						else { push(@{ $Regions[$region]{"ProfondeurS_Region_Ponderee"} }, $Regions[$region]{$patient}{"Profondeur_Region_Ponderee"}); }
						}
					}
				# Nous calculons pour chaque region la mediane/moyenne de profondeur (ponderee pour chaque patient par rapport aux autres regions)
				if(@{ $Regions[$region]{"ProfondeurS_Region_Ponderee"} }) {
					 $Regions[$region]{"Norm_Profondeurs_Ponderees"} = norm($norm,\@{ $Regions[$region]{"ProfondeurS_Region_Ponderee"} });
					if ($norm eq "std") {
						my$sqtotal = 0;
						foreach (@{ $Regions[$region]{"ProfondeurS_Region_Ponderee"} })
							{ $sqtotal += ($Regions[$region]{"Norm_Profondeurs_Ponderees"}-$_)**2; }
						$Regions[$region]{"Std_Profondeurs_Ponderees"} = ($sqtotal / (scalar@{ $Regions[$region]{"ProfondeurS_Region_Ponderee"} }-1))**0.5;
					}
				}

				# Un second parcours permet l'appel de CNV (en comparant la profondeur pour un patient par la profondeur moyenne de la region)
				for (my $patient = 0 ; $patient < $nombre_patients ; $patient++) {

					if ($norm eq "std") { $prof_Moyenne_Inter = 0; }
					else  { $prof_Moyenne_Inter = 1; }

					unless ($Patients{$patient}{"ecarte"}) {

						if ($Regions[$region]{"Norm_Profondeurs_Ponderees"}) {
							if($Patients{$patient}{"Sexe"} eq "F") {
								if ($Regions[$region]{"Chromosome"} !~ m/^(chrY|Y)$/) {
									if ($norm eq "std") {
										$Regions[$region]{$patient}{"Moyenne_Inter"} = ($Regions[$region]{$patient}{"Profondeur_Region_Ponderee"}-$Regions[$region]{"Norm_Profondeurs_Ponderees"})/$Regions[$region]{"Std_Profondeurs_Ponderees"};
									} else {
									$Regions[$region]{$patient}{"Moyenne_Inter"} = $Regions[$region]{$patient}{"Profondeur_Region_Ponderee"}/$Regions[$region]{"Norm_Profondeurs_Ponderees"};
									}
									$prof_Moyenne_Inter = $Regions[$region]{$patient}{"Moyenne_Inter"};
									print SORTIE $Regions[$region]{$patient}{"Moyenne_Inter"}."\t";
								} else {
									#les femmes ne sont pas considérées pour l'appel de CNV des régions du Y
									print SORTIE "NA\t";
								}
							} else {
								#*2 for chrX
								if ($Regions[$region]{"Chromosome"} =~ m/^(chrX|X)$/) {
									if ($norm eq "std") {
										$Regions[$region]{$patient}{"Moyenne_Inter"} = ($Regions[$region]{$patient}{"Profondeur_Region_Ponderee"}*2-$Regions[$region]{"Norm_Profondeurs_Ponderees"})/$Regions[$region]{"Std_Profondeurs_Ponderees"};
									} else {
									$Regions[$region]{$patient}{"Moyenne_Inter"} = $Regions[$region]{$patient}{"Profondeur_Region_Ponderee"}*2/$Regions[$region]{"Norm_Profondeurs_Ponderees"};
									}
								} else {
									$Regions[$region]{$patient}{"Moyenne_Inter"} = $Regions[$region]{$patient}{"Profondeur_Region_Ponderee"}/$Regions[$region]{"Norm_Profondeurs_Ponderees"};
									if ($norm eq "std") {
										$Regions[$region]{$patient}{"Moyenne_Inter"} = ($Regions[$region]{$patient}{"Profondeur_Region_Ponderee"}-$Regions[$region]{"Norm_Profondeurs_Ponderees"})/$Regions[$region]{"Std_Profondeurs_Ponderees"};
									} else {
									$Regions[$region]{$patient}{"Moyenne_Inter"} = $Regions[$region]{$patient}{"Profondeur_Region_Ponderee"}/$Regions[$region]{"Norm_Profondeurs_Ponderees"};
									}
								}
								$prof_Moyenne_Inter = $Regions[$region]{$patient}{"Moyenne_Inter"};
								print SORTIE $Regions[$region]{$patient}{"Moyenne_Inter"}."\t";
							}
							if ($prof_Moyenne_Inter < $seuil_deletion) {
								$nb_deletions++;
								$nb_evts++;
								$CNV{$patient}++;
								$results{$patient}{$region} = "DEL\t".sprintf("%.3f",$Regions[$region]{$patient}{"Moyenne_Inter"});
							} elsif ($prof_Moyenne_Inter > $seuil_duplication) {
								$nb_duplications++;
								$nb_evts++;
								$CNV{$patient}++;
								$results{$patient}{$region} = "DUP\t".sprintf("%.3f",$Regions[$region]{$patient}{"Moyenne_Inter"});
							}
						} else {
							$Regions[$region]{"Appel"} = "Absence_Donnees";
							print SORTIE "NA\t";
						}
					} else {
						print SORTIE "NA\t";
					}
				}
	
				if(@{ $Regions[$region]{"ProfondeurS_Region_Ponderee"} } && (($nb_evts/scalar@{ $Regions[$region]{"ProfondeurS_Region_Ponderee"} }) > $seuil_region) && ($nb_parcours > 0)) {
					print SORTIE "CNV_Recurrent";
					print LOG "Region ecartee \: ". $Regions[$region]{"Chromosome"}."\t".$Regions[$region]{"Borne_5P"}."\t".$Regions[$region]{"Borne_3P"}."\t".$Regions[$region]{"Gene"}."\n";
					$Regions[$region]{"Appel"} = "CNV_Recurrent";
					$regions_ecartees++;
					$continuer = 1;

				} else {		
					print SORTIE "OK";
				}
				print SORTIE "\n";
			}


		# SI LA REGION EST "MOCHE"
		} else {

			for (my $patient = 0 ; $patient < $nombre_patients ; $patient++) {
				print SORTIE "NA\t";
			}
			print SORTIE $Regions[$region]{"Appel"}."\n";

		}

	}

	print LOG "Nombre de regions ecartees lors de cette iteration \: ".$regions_ecartees."\n\n\n";	

	my $patients_ecartes = 0;

	for (my $patient = 0 ; $patient < $nombre_patients ; $patient++) {

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
my(%uniqReg,%RegInArray,%regionOrder,%regionIndice);
for my$region (0..$#Regions) {
	if (!exists $uniqReg{ $Regions[$region]{"Chromosome"}."-".$Regions[$region]{"Borne_5P"}."-".$Regions[$region]{"Borne_3P"} }) {
		push (@{ $RegInArray{ $Regions[$region]{"Chromosome"} } }, $region);
		$uniqReg{ $Regions[$region]{"Chromosome"}."-".$Regions[$region]{"Borne_5P"}."-".$Regions[$region]{"Borne_3P"} } = 1;
		}
	if ($Regions[$region]{"Gene"} eq "NA")
		{ $Regions[$region]{"label"} = $Regions[$region]{"Borne_5P"}."-".$Regions[$region]{"Borne_3P"}; }
	else { 
		my@tab = split(/:|,/,$Regions[$region]{"Gene"});
		my$label = "";
		if (length($tab[0]) > 20) {
			my@tab2 =  split(//, $tab[0]);
			for my$i (0..20) { $label .= $tab2[$i]; }
			}
		else { $label = $tab[0]; }
		$Regions[$region]{"label"} = $label;
		$Regions[$region]{"geneID"} = $tab[0];
		}
	}
for my$Chrom (keys%RegInArray) {
	##sort by region ends (in case several regions with same start) then by starts
	@{ $regionOrder{$Chrom} } = sort{$Regions[$a]{"Borne_3P"}<=>$Regions[$b]{"Borne_3P"}}@{ $RegInArray{$Chrom} }; 
	@{ $regionOrder{$Chrom} } = sort{$Regions[$a]{"Borne_5P"}<=>$Regions[$b]{"Borne_5P"}}@{ $RegInArray{$Chrom} };
	for (my$i=0;$i<scalar@{ $regionOrder{$Chrom} };$i++) { $regionIndice{ $regionOrder{$Chrom}[$i] } = $i; }
	}


##print CNV foreach sample
#$results{$patient}{$region} = "DUP (Moyenne_Inter) (Borne_3P-Borne_5P bp)";
my%orderedCNV;
foreach my$patient (keys%results) {
	my(%uniqCNV,%CNVinArray);
	foreach my$region (sort{$a<=>$b}keys%{ $results{$patient} }) {
		if (!exists $uniqCNV{ $Regions[$region]{"Chromosome"}."-".$Regions[$region]{"Borne_5P"}."-".$Regions[$region]{"Borne_3P"} }) { 
			push (@{ $CNVinArray{ $Regions[$region]{"Chromosome"} } }, $region);
			$uniqCNV{ $Regions[$region]{"Chromosome"}."-".$Regions[$region]{"Borne_5P"}."-".$Regions[$region]{"Borne_3P"} } = 1;
			}
		}
	foreach my$Chrom (keys%CNVinArray) {
		##sort by region ends (in case several regions with same start) then by starts
		@{ $orderedCNV{$patient}{$Chrom} } = sort{$Regions[$a]{"Borne_3P"}<=>$Regions[$b]{"Borne_3P"}}@{ $CNVinArray{$Chrom} }; 
		@{ $orderedCNV{$patient}{$Chrom} } = sort{$Regions[$a]{"Borne_5P"}<=>$Regions[$b]{"Borne_5P"}}@{ $CNVinArray{$Chrom} }; 
		}
	}

unless($minCNV) { $minCNV = 1; }
my%Result2;
foreach my$patient (keys%results) {
	open (CNV1,">$outdir/CNV_$Patients{$patient}{ID}.summary.txt") or die "Pb lors de l'ecriture du fichier sortie $!\n";
	print CNV1 "#Region\tCNV\tN_intervals\n";
	open (CNV2,">$outdir/CNV_$Patients{$patient}{ID}.allIntervals.txt") or die "Pb lors de l'ecriture du fichier sortie $!\n";
	print CNV2 "#Chrom\tStart\tEnd\tInfo\tinterval_Order\tCNV_type\tratio_to_$norm\tinterval_Length\n";
	foreach my$Chrom (@ChrOrder) {
		if (exists $orderedCNV{$patient}{$Chrom}) {
			my$r=0;		##index in @{ $orderedCNV{$patient}{$Chrom} }
			while ($r < scalar@{ $orderedCNV{$patient}{$Chrom} } ) {
				my$ok=1; my$i=0; my$cnvOK=0;
				my@nextReg=($orderedCNV{$patient}{$Chrom}[$r]);
				my@tab = split(/\s/,$results{$patient}{$orderedCNV{$patient}{$Chrom}[$r]});
				while ($ok) {
					$i++;
					if (exists $regionOrder{$Chrom}[($regionIndice{$orderedCNV{$patient}{$Chrom}[$r]} + $i)]) {
						push(@nextReg, $regionOrder{$Chrom}[($regionIndice{$orderedCNV{$patient}{$Chrom}[$r]} + $i)]);
						if (exists $results{$patient}{$nextReg[$i]}) {
							my@tab2 = split(/\s/, $results{$patient}{$nextReg[$i]});
							if($tab2[0] eq $tab[0]) { $cnvOK = $i; }	##same CNV type
							else { $ok=0; }
							}
						else {
							unless (exists $Regions[$nextReg[$i]]{"Appel"}) { $ok=0; }
							}
						}
					else { $ok=0; }
					}
				if ($cnvOK >= ($minCNV-1)) {
					if ($cnvOK == 0) {
						$Result2{$patient}{$nextReg[0]} = $results{$patient}{$nextReg[0]};
						if (exists $regionIndice{$nextReg[0]}) {
							print CNV1 $Chrom.":".$Regions[$nextReg[0]]{"Borne_5P"}."-".$Regions[$nextReg[0]]{"Borne_3P"}."\t".$tab[0]."\t1\n";
							print CNV2 $Chrom."\t".$Regions[$nextReg[0]]{"Borne_5P"}."\t".$Regions[$nextReg[0]]{"Borne_3P"}."\t".$Regions[$nextReg[0]]{"label"}."\t".($regionIndice{$nextReg[0]}+1)."\t".$results{$patient}{$nextReg[0]}."\t".($Regions[$nextReg[0]]{"Borne_3P"}-$Regions[$nextReg[0]]{"Borne_5P"})." bp\n\n";
							}
						}
					else {
						print CNV1 $Chrom.":".$Regions[$nextReg[0]]{"Borne_5P"}."-".$Regions[$nextReg[$cnvOK]]{"Borne_3P"}."\t".$tab[0]."\t".($cnvOK+1)."\n";
						for (my$j=0;$j<=$cnvOK;$j++) {
							if (exists $results{$patient}{$nextReg[$j]}) {
								$Result2{$patient}{$nextReg[$j]} = $results{$patient}{$nextReg[$j]};
								if (exists $regionIndice{$nextReg[$j]})
									{ print CNV2 $Chrom."\t".$Regions[$nextReg[$j]]{"Borne_5P"}."\t".$Regions[$nextReg[$j]]{"Borne_3P"}."\t".$Regions[$nextReg[$j]]{"label"}."\t".($regionIndice{$nextReg[$j]}+1)."\t".$results{$patient}{$nextReg[$j]}."\t".($Regions[$nextReg[$j]]{"Borne_3P"}-$Regions[$nextReg[$j]]{"Borne_5P"})." bp\n"; }
								}
							else	{ 
								$Result2{$patient}{$nextReg[$j]} = "NA";
								if (exists $regionIndice{$nextReg[$j]}) {
									print CNV2 $Chrom."\t".$Regions[$nextReg[$j]]{"Borne_5P"}."\t".$Regions[$nextReg[$j]]{"Borne_3P"}."\t".$Regions[$nextReg[$j]]{"label"}."\t".($regionIndice{$nextReg[$j]}+1)."\tNA\t";
									if (exists $Regions[$nextReg[$j]]{$patient}{"Moyenne_Inter"})
										{ print CNV2 sprintf("%.3f",$Regions[$nextReg[$j]]{$patient}{"Moyenne_Inter"}); }
									else { print CNV2 "na"; }
									print CNV2 "\t".($Regions[$nextReg[$j]]{"Borne_3P"}-$Regions[$nextReg[$j]]{"Borne_5P"})." bp"."\n";
									}
								}
							}
						print CNV2 "\n";
						}
					}
				if ($i > 1) { $r += ($i-1); }
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

}

####################

sub norm {
my($norm,$h1)=@_;
my@allDepth=@$h1;
my$normDepth;
if ($norm eq "med") {
	@allDepth = sort{$a<=>$b}@allDepth;
	#odd?
	if(scalar@allDepth%2) 
		{ $normDepth = $allDepth[int(scalar@allDepth/2)]; }
	#even
	else { $normDepth = ( $allDepth[int(scalar@allDepth/2)-1] + $allDepth[int(scalar@allDepth/2)] )/2; }
	}
else {
	foreach (@allDepth)
		{ $normDepth += $_; } 
	$normDepth /= scalar@allDepth;
	}
return($normDepth);
}


####################

sub graphChr1 {

my($nGraf,$outdir,$norm,$seuil_deletion,$seuil_duplication,$file,$maxDepthGraph,$h1,$h2,$h3,$h4,$h5,$h6,$h7)= @_;
my@Files = @$h1;
my%sampleName = %$h2;
my@Regions = @$h3;
my@ChromOrder = @$h4;
my%RegionOrder = %$h5;
my%Patients = %$h6;
my%Results = %$h7;

print "cmdR, for sample $sampleName{$file}\n";

my$Nbr_Chr= scalar(keys%RegionOrder);
if ($nGraf eq "max") { $nGraf = $Nbr_Chr; }

my$maxX=0;
for my$Chrom (keys%RegionOrder) {	
	if (scalar@{ $RegionOrder{$Chrom} } > $maxX) { $maxX = scalar@{ $RegionOrder{$Chrom} }; }
	}

my$cmdR = "";
my$c=1; #chr iteration
my$n=1; #chr iteration, stepped back to 0 each time a graph is done
my$N=1; #graph iteration
my$maxGeneLab=200;
my$maxGeneSep=500;

foreach my$Chrom (@ChromOrder) { 
	
	$Chrom =~ s/^chr//;
	
	if (exists $RegionOrder{$Chrom}) {

		my$maxYsup=$seuil_duplication; my$maxYinf=$seuil_deletion;	
		foreach my$region (@{ $RegionOrder{$Chrom} }) {
			foreach my$f (@Files) {
				if (exists$Regions[$region]{$f}{"Moyenne_Inter"}) {
					if ($Regions[$region]{$f}{"Moyenne_Inter"} > $maxYsup)
						{ $maxYsup = $Regions[$region]{$f}{"Moyenne_Inter"}; }
					if ($norm eq "std") {
						if ($Regions[$region]{$f}{"Moyenne_Inter"} < $maxYinf)
							{ $maxYinf = $Regions[$region]{$f}{"Moyenne_Inter"}; }
						}
					}
				}
			}
		if ($maxDepthGraph && $maxYsup > $maxDepthGraph) { $maxYsup = $maxDepthGraph; }
		if ($norm eq "std") {
			if ($maxDepthGraph && $maxYinf < (-$maxDepthGraph)) { $maxYinf = (-$maxDepthGraph); }
			}

		if ($norm eq "std") { $cmdR .= "par(fig=c(0,1,".(1-(($n-0.05)/$nGraf)).",".(1-(($n-0.95)/$nGraf))."), new=TRUE)
	plot (c(0,0), xlim=c(0,$maxX), ylim=c($maxYinf,$maxYsup), type =\"n\", main=\"chrom: $Chrom\", xlab=\"\", ylab=\"depth_ratio_to_$norm\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n"; }
		else { $cmdR .= "par(fig=c(0,1,".(1-(($n-0.05)/$nGraf)).",".(1-(($n-0.95)/$nGraf))."), new=TRUE)
	plot (c(0,0), xlim=c(0,$maxX), ylim=c(0,$maxYsup), type =\"n\", main=\"chrom: $Chrom\", xlab=\"\", ylab=\"depth_ratio_to_$norm\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n"; }

		my$Nbr_Reg = scalar@{ $RegionOrder{$Chrom} };

		##gene separations
		my$currentGene=""; my$tmpTxt=""; my$Nbr_gene=0;
		for (my$r=0;$r<$Nbr_Reg;$r++) {
			if ($Regions[$RegionOrder{$Chrom}[$r]]{"Gene"}) {
				if ($Regions[$RegionOrder{$Chrom}[$r]]{"Gene"} ne "NA" && $Regions[$RegionOrder{$Chrom}[$r]]{"Gene"} ne $currentGene)  {
					$tmpTxt .= "abline(v=".($r+0.5).", col=\"blue\", lty = \"dotted\", lwd=2)\n";
					$currentGene = $Regions[$RegionOrder{$Chrom}[$r]]{"Gene"};
					}
				$Nbr_gene++;
				}
			}
		if ($Nbr_gene < $maxGeneSep) { $cmdR .= $tmpTxt; }

		my$Nbr_CNV=0;
		for (my$r=0;$r<$Nbr_Reg;$r++) {
			if (exists$Results{$file}{$RegionOrder{$Chrom}[$r]})
				{ $Nbr_CNV++; }
			}

		##region labels
		my@printReg=();
		if ($Nbr_Reg<$maxGeneLab) {
			##in red if invalid
			for (my$r=0;$r<$Nbr_Reg;$r++) { 
				if (defined $Regions[$RegionOrder{$Chrom}[$r]]{"Appel"}) { push(@printReg,$r) ; }
				}
			if (@printReg) {
				$cmdR .= "axis(1, at=c(";
				foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
				chop $cmdR;
				$cmdR .="), labels=c(";
				foreach my$r (@printReg) { $cmdR .= "\"".$Regions[$RegionOrder{$Chrom}[$r]]{"label"}."\","; }
				chop $cmdR;
				$cmdR .= "), col.axis=\"red\", las=2)\n";
				}
			@printReg=();
			for (my$r=0;$r<$Nbr_Reg;$r++) { 
				if(!defined $Regions[$RegionOrder{$Chrom}[$r]]{"Appel"}) { push(@printReg,$r) ; }
				}
			if (@printReg) {
				$cmdR .= "axis(1, at=c(";
				foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
				chop $cmdR;
				$cmdR .="), labels=c(";
				foreach my$r (@printReg) { $cmdR .= "\"".$Regions[$RegionOrder{$Chrom}[$r]]{"label"}."\","; }
				chop $cmdR;
				$cmdR .= "), col.axis=\"black\", las=2)\n";
				}
			}
		else {
			if ($Nbr_gene < $maxGeneLab) { #&& ($Nbr_gene+$Nbr_CNV)>=$maxGeneLab) {
				$currentGene="";
				for (my$r=0;$r<$Nbr_Reg;$r++) {
					if ($Regions[$RegionOrder{$Chrom}[$r]]{"Gene"}) {
						if ($Regions[$RegionOrder{$Chrom}[$r]]{"Gene"} ne "NA" && $Regions[$RegionOrder{$Chrom}[$r]]{"Gene"} ne $currentGene)  { 
							push(@printReg,$r); 
							$currentGene = $Regions[$RegionOrder{$Chrom}[$r]]{"Gene"};
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
				foreach my$r (@printReg) { $cmdR .= "\"".$Regions[$RegionOrder{$Chrom}[$r]]{"label"}."\","; }
				chop $cmdR;
				$cmdR .= "), col.axis=\"black\", las=2)\n";
				}
			}

		#all not target sample lines (black):
		for (my$f=0;$f<scalar@Files;$f++) {
			unless ($Files[$f] eq $file) {
				my$r1=0;
				while ($r1<$Nbr_Reg) {
					my$r2=$r1;
					while ($r2<$Nbr_Reg) {
						if (exists $Regions[$RegionOrder{$Chrom}[$r2]]{$Files[$f]}{"Moyenne_Inter"}) { $r2++;}
						else { last; }
						}
					if (($r2-1) > $r1) {
						$cmdR .= "lines( c(";
						for (my$r=$r1;$r<$r2;$r++)
							{ $cmdR .= ($r+1).","; }
						chop $cmdR;
						$cmdR .= "), c(";
						for (my$r=$r1;$r<$r2;$r++)
							{ $cmdR .= $Regions[$RegionOrder{$Chrom}[$r]]{$Files[$f]}{"Moyenne_Inter"}.","; }
						chop $cmdR;
						$cmdR .= "), type =\"l\", lwd=2, col=\"black\")\n";
						}
					elsif (($r2-1) == $r1) {
						$cmdR .= "lines( c(".($r1+1)."), c(".$Regions[$RegionOrder{$Chrom}[$r1]]{$Files[$f]}{"Moyenne_Inter"}."), type =\"p\", lwd=2, col=\"black\")\n";
						}
					$r1 = ($r2+1);
					}
				}
			}
		#target sample line (green)
		my$r1=0;
		while ($r1<$Nbr_Reg) {
			my$r2=$r1;
			while ($r2<$Nbr_Reg) {
				if (exists$Regions[$RegionOrder{$Chrom}[$r2]]{$file}{"Moyenne_Inter"}) { $r2++;}
				else { last; }
				}
			if (($r2-1) > $r1) {
				$cmdR .= "lines( c(";
				for (my$r=$r1;$r<$r2;$r++)
					{ $cmdR .= ($r+1).","; }
				chop $cmdR;
				$cmdR .= "), c(";
				for (my$r=$r1;$r<$r2;$r++)
					{ $cmdR .= $Regions[$RegionOrder{$Chrom}[$r]]{$file}{"Moyenne_Inter"}.","; }
				chop $cmdR;
				$cmdR .= "), type =\"l\", lwd=2, col=\"green\")\n";
				}
			elsif (($r2-1) == $r1) {
				$cmdR .= "lines( c(".($r1+1)."), c(".$Regions[$RegionOrder{$Chrom}[$r1]]{$file}{"Moyenne_Inter"}."), type =\"p\", lwd=2, col=\"green\")\n";
				}
			$r1 = ($r2+1);
			}
		#points for CNVs
		my$points .= "points( c(";
		for (my$r=0;$r<$Nbr_Reg;$r++) {
			if (exists$Results{$file}{$RegionOrder{$Chrom}[$r]} && exists$Regions[$RegionOrder{$Chrom}[$r]]{$file}{"Moyenne_Inter"})
				{ $points .= ($r+1).","; }
			}
		chop $points;
		$points .= "), c(";
		foreach my$region (@{ $RegionOrder{$Chrom} }) {
			if (exists$Results{$file}{$region} && exists$Regions[$region]{$file}{"Moyenne_Inter"})
				{ $points .= $Regions[$region]{$file}{"Moyenne_Inter"}.","; }
			}
		chop $points;
		$points .= "), type =\"p\", pch = 16, lwd=3, col=\"red\")\n";
		if ($Nbr_CNV) { $cmdR .= $points; }

		#thereshold lines
		$cmdR .= "abline(h=$seuil_deletion, col=\"darkgrey\", lty = \"dotted\", lwd=2)\n";
		$cmdR .= "abline(h=$seuil_duplication, col=\"darkgrey\", lty = \"dotted\", lwd=2)\n";
		if ($norm eq "std") { $cmdR .= "abline(h=0, col=\"darkgrey\", lty = \"dotted\", lwd=2)\n"; }
		else { $cmdR .= "abline(h=1, col=\"darkgrey\", lty = \"dotted\", lwd=2)\n"; }

		if ($c==$Nbr_Chr || $n==$nGraf) {
			open (CMDR, ">$outdir/$sampleName{$file}\_temp.R") || die;
			print CMDR "#!/usr/bin/env Rscript\n\n" ;
			if ($nGraf==$Nbr_Chr) { print CMDR "png(\"".$outdir."/CNV_$sampleName{$file}.png\", 1500, ".($nGraf*400).")\n
	par(mfrow=c($nGraf,1))\n"; }
			else {
				if ($N>1) { print CMDR "png(\"".$outdir."/CNV_$sampleName{$file}\_$N.png\", 1500, ".($nGraf*400).")\n
	par(mfrow=c($nGraf,1))\n"; }
				else { print CMDR "png(\"".$outdir."/CNV_$sampleName{$file}\_$N.png\", 1500, ".($n*400).")\n
	par(mfrow=c($n,1))\n"; }
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
my%results = %$h5;	#$results{$patient}{$region}


my$maxX=0;
for my$Chrom (keys%regionOrder) {	
	if (scalar@{ $regionOrder{$Chrom} } > $maxX) { $maxX = scalar@{ $regionOrder{$Chrom} }; }
	}

my$Nbr_Chr= scalar(keys%regionOrder);
print "cmdR,  for $Patients{$patient}{ID}\n";
open (CMDR, ">$outdir/$Patients{$patient}{ID}\_temp.R") || die;
print CMDR "#!/usr/bin/env Rscript\n\n" ;
print CMDR "png(\"".$outdir."/CNV_$Patients{$patient}{ID}.png\", 1500, ".($Nbr_Chr*400).")\n
par(mfrow=c($Nbr_Chr,1))\n";

my$cmdR = "";
my$c=0; #chr iteration
my$maxGeneLab=200;
my$maxGeneSep=500;

for my$Chrom (@ChrOrder) {

	my$maxYsup=$seuil_duplication; my$maxYinf=$seuil_deletion;	
	foreach my$region (@{ $regionOrder{$Chrom} }) {
		for (my$p=0;$p<scalar(keys%Patients);$p++) {
			if (exists $Regions[$region]{$p}{"Moyenne_Inter"}) {
				if ($Regions[$region]{$p}{"Moyenne_Inter"} > $maxYsup)
					{ $maxYsup = $Regions[$region]{$p}{"Moyenne_Inter"}; }
				if ($norm eq "std") {
					if ($Regions[$region]{$p}{"Moyenne_Inter"} < $maxYinf)
						{ $maxYinf = $Regions[$region]{$p}{"Moyenne_Inter"}; }
					}
				}
			}
		}
	if ($maxDepthGraph && $maxYsup > $maxDepthGraph) { $maxYsup = $maxDepthGraph; }
	if ($norm eq "std") {
		if ($maxDepthGraph && $maxYinf < (-$maxDepthGraph)) { $maxYinf = (-$maxDepthGraph); }
		}

	my$ChrName = $Chrom; $ChrName =~ s/^chr//;

	if ($norm eq "std") { $cmdR .= "par(fig=c(0,1,".(1-(($c+0.95)/$Nbr_Chr)).",".(1-(($c+0.05)/$Nbr_Chr))."), new=TRUE)
plot (c(0,0), xlim=c(0,$maxX), ylim=c($maxYinf,$maxYsup), type =\"n\", main=\"chrom: $ChrName\", xlab=\"\", ylab=\"depth_ratio_to_$norm\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n"; }
	else { $cmdR .= "par(fig=c(0,1,".(1-(($c+0.95)/$Nbr_Chr)).",".(1-(($c+0.05)/$Nbr_Chr))."), new=TRUE)
plot (c(0,0), xlim=c(0,$maxX), ylim=c(0,$maxYsup), type =\"n\", main=\"chrom: $ChrName\", xlab=\"\", ylab=\"depth_ratio_to_$norm\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n"; }

	my$Nbr_Reg = scalar@{ $regionOrder{$Chrom} };

	#gene vertical separations
	my$currentGene=""; my$tmpTxt=""; my$Nbr_gene=0;
	for (my$r=0;$r<$Nbr_Reg;$r++) {
		if ($Regions[$regionOrder{$Chrom}[$r]]{"Gene"} ne "NA" && $Regions[$regionOrder{$Chrom}[$r]]{"geneID"} ne $currentGene) {
			$tmpTxt .= "abline(v=".($r+0.5).", col=\"blue\", lty = \"dotted\", lwd=2)\n";
			$currentGene = $Regions[$regionOrder{$Chrom}[$r]]{"geneID"};
			$Nbr_gene++;
			}
		}
	if ($Nbr_gene < $maxGeneSep) { $cmdR .= $tmpTxt; }

	my$Nbr_CNV=0;
	for (my$r=0;$r<$Nbr_Reg;$r++) {
		if (exists$results{$patient}{$regionOrder{$Chrom}[$r]})
			{ $Nbr_CNV++; }
		}

	#x labels
	my@printReg=();
	if ($Nbr_Reg<$maxGeneLab) {
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
			$cmdR .= "), col.axis=\"red\", las=2)\n";
			}
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
			$cmdR .= "), col.axis=\"black\", las=2)\n";
			}
		}
	else {
		if ($Nbr_gene<$maxGeneLab) {	# && ($Nbr_gene+$Nbr_CNV)>=$maxGeneLab) {
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
			$cmdR .= "), col.axis=\"black\", las=2)\n";
			}
		}

	#all not target sample lines (black):
	for (my$p=0;$p<scalar(keys%Patients);$p++) {
		unless ($p == $patient) {

			my$r1=0;
			while ($r1<$Nbr_Reg) {
				my$r2=$r1;
				while ($r2<$Nbr_Reg) {
					if (exists$Regions[$regionOrder{$Chrom}[$r2]]{$p}{"Moyenne_Inter"}) { $r2++;}
					else { last; }
					}
				if (($r2-1) > $r1) {
					$cmdR .= "lines( c(";
					for (my$r=$r1;$r<$r2;$r++)
						{ $cmdR .= ($r+1).","; }
					chop $cmdR;
					$cmdR .= "), c(";
					for (my$r=$r1;$r<$r2;$r++)
						{ $cmdR .= $Regions[$regionOrder{$Chrom}[$r]]{$p}{"Moyenne_Inter"}.","; }
					chop $cmdR;
					$cmdR .= "), type =\"l\", lwd=2, col=\"black\")\n";
					}
				elsif (($r2-1) == $r1) {
					$cmdR .= "lines( c(".($r1+1)."), c(".$Regions[$regionOrder{$Chrom}[$r1]]{$p}{"Moyenne_Inter"}."), type =\"p\", lwd=2, col=\"black\")\n";
					}
				$r1 = ($r2+1);
				}
			}
		}
	#target sample line (green):
	my$r1=0;
	while ($r1<$Nbr_Reg) {
		my$r2=$r1;
		while ($r2<$Nbr_Reg) {
			if (exists$Regions[$regionOrder{$Chrom}[$r2]]{$patient}{"Moyenne_Inter"}) { $r2++;}
			else { last; }
			}
		if (($r2-1) > $r1) {
			$cmdR .= "lines( c(";
			for (my$r=$r1;$r<$r2;$r++)
				{ $cmdR .= ($r+1).","; }
			chop $cmdR;
			$cmdR .= "), c(";
			for (my$r=$r1;$r<$r2;$r++)
				{ $cmdR .= $Regions[$regionOrder{$Chrom}[$r]]{$patient}{"Moyenne_Inter"}.","; }
			chop $cmdR;
			$cmdR .= "), type =\"l\", lwd=2, col=\"green\")\n";
			}
		elsif (($r2-1) == $r1) {
			$cmdR .= "lines( c(".($r1+1)."), c(".$Regions[$regionOrder{$Chrom}[$r1]]{$patient}{"Moyenne_Inter"}."), type =\"p\", lwd=2, col=\"green\")\n";
			}
		$r1 = ($r2+1);
		}

	#points for CNVs
	my$points .= "points( c(";
	for (my$r=0;$r<$Nbr_Reg;$r++) {
		if (exists$results{$patient}{$regionOrder{$Chrom}[$r]} && exists$Regions[$regionOrder{$Chrom}[$r]]{$patient}{"Moyenne_Inter"})
			{ $points .= ($r+1).","; }
		}
	chop $points;
	$points .= "), c(";
	foreach my$region (@{ $regionOrder{$Chrom} }) {
		if (exists$results{$patient}{$region} && exists$Regions[$region]{$patient}{"Moyenne_Inter"})
			{ $points .= $Regions[$region]{$patient}{"Moyenne_Inter"}.","; }
		}
	chop $points;
	$points .= "), type =\"p\", pch = 16, lwd=3, col=\"red\")\n";
	if ($Nbr_CNV) { $cmdR .= $points; }

	#horizontal threshold lines:
	$cmdR .= "abline(h=$seuil_deletion, col=\"darkgrey\", lty = \"dotted\", lwd=2)\n";
	$cmdR .= "abline(h=$seuil_duplication, col=\"darkgrey\", lty = \"dotted\", lwd=2)\n";
	if ($norm eq "std") { $cmdR .= "abline(h=0, col=\"darkgrey\", lty = \"dotted\", lwd=2)\n"; }
	else { $cmdR .= "abline(h=1, col=\"darkgrey\", lty = \"dotted\", lwd=2)\n"; }

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
