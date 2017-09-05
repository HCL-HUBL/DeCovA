package Bio::NGS::HCL::DeCovA::CNV_tool;


use strict;
use warnings;

#############
#@hashSub = Bio::NGS::HCL::DeCovA::CNV_tool::CNV_detect(\@Files,\%sName2,"$outdir/CNV_analysis",$nGraf,$fichier_sexe,\%CNV_opt,$refBedLines,\@ChromOrder);
sub CNV_detect {

my($Files_r,$sampleName_r,$outdir,$nGraf,$fichier_sexe,$CNV_opt_r,$Regions_r,$ChromOrder_r) = @_;

my$center = ${$CNV_opt_r}{"center"};
my$spread = ${$CNV_opt_r}{"spread"};
my$center_test = ${$CNV_opt_r}{"center_test"};
my$spread_test = ${$CNV_opt_r}{"spread_test"};
my$seuil_del = ${$CNV_opt_r}{"seuil_del"};
my$seuil_dup = ${$CNV_opt_r}{"seuil_dup"};
my$spread_del = ${$CNV_opt_r}{"spread_del"};
my$spread_dup = ${$CNV_opt_r}{"spread_dup"};
my$range = ${$CNV_opt_r}{"range"};		##sample within med-$range*quartile1 and med+$range*quartile3 will be included for avg+std calculation
my$seuil_region = ${$CNV_opt_r}{"seuil_region"};
if ($seuil_region < (1/scalar@{$Files_r})) { $seuil_region = (1/scalar@{$Files_r}); }
my$seuil_patient = ${$CNV_opt_r}{"seuil_patient"};
if ($seuil_patient < (1/scalar@{$Regions_r})) { $seuil_patient = (1/scalar@{$Regions_r}); }
my$minCov = ${$CNV_opt_r}{"seuil_cov"};
my$minDP = ${$CNV_opt_r}{"min_DP"} ;
my$minCNV = ${$CNV_opt_r}{"min_following_CNV"};
my$maxNonCNV = ${$CNV_opt_r}{"max_Non_CNV"};
my$maxNonCNVrate = ${$CNV_opt_r}{"max_Non_CNV_rate"};
my$ratioByGender = ${$CNV_opt_r}{"ratioByGender"};
my$RefDepth = ${$CNV_opt_r}{"RefDepth"};
my$RefNoChrY = ${$CNV_opt_r}{"RefNoChrY"};
my$graphByChr = ${$CNV_opt_r}{"chromGraph"};

my@cnvFields = @{ ${$CNV_opt_r}{"fields"} };	#min, max, avg , std , med
my%cnvStats;
foreach (@cnvFields) { $cnvStats{$_} = 1; }
$cnvStats{$center} = 1;
if (defined $spread) { $cnvStats{$spread} = 1; }
if (exists $cnvStats{"std"}) { $cnvStats{"avg"} = 1;  }
if (defined $range || (defined $spread && $spread eq "Qtile") || exists $cnvStats{"Q1"} || exists $cnvStats{"Q3"}) {
	$cnvStats{"med"} = 1;
	$cnvStats{"Qtile"} = 1;
	}
my$sortDepths = 0;
if (exists $cnvStats{"med"} || exists $cnvStats{"min"} || exists $cnvStats{"max"}) { $sortDepths = 1; }

my$nCol = scalar(split(/\t/,${$Regions_r}[0]{"allLine"}));

my%Patients;

#read gender file
if ($fichier_sexe) {
	open(PATIENTS, "$fichier_sexe") or die "Fichier $fichier_sexe impossible a lire\n";
	foreach my $ligne (<PATIENTS>) {
		$ligne =~ m/(\S+)\s+(\S+)/;
		my$ok="";
		foreach my$file (@{$Files_r}) {
			if ($1 eq ${$sampleName_r}{$file}) {
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
foreach my$file (@{$Files_r}) {
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
#$bedLines[$i]{$file}{"cov"}
#$bedLines[$i]{$file}{"mean"}
#$bedLines[$i]{$file}{"tot"}
for (my $r = 0 ; $r < scalar@{$Regions_r} ; $r++) {

	my $region_a_conserver = 1;
	if ($minCov) {
		foreach my$file (@{$Files_r}) {
			if ( ${$Regions_r}[$r]{$file}{"cov"} <= $minCov) { $region_a_conserver = 0; }
			else { $region_a_conserver = 1; last; }
			}
		}
	if (defined $minDP) {	#def = 0
		my(@allDepth,$centerDepth);
		foreach (@{$Files_r}) { push(@allDepth, ${$Regions_r}[$r]{$_}{"mean"}); }
		#$centerDepth = mediane(\@allDepth);
		$centerDepth = average(\@allDepth);
		if ( $centerDepth <= $minDP) { $region_a_conserver = 0; }
		#foreach my$file (@{$Files_r}) {
		#	if ( ${$Regions_r}[$r]{$file}{"mean"} < $minDP) { $region_a_conserver = 0; }
		#	else { $region_a_conserver = 1; last; }
		#	}
		}
	# si region à conserver, calcul Profondeur_Autosomes/gonosomes
	if ($region_a_conserver) {
		foreach my$file (@{$Files_r}) {
			if (${$Regions_r}[$r]{"Chrom"} !~ /^(chr[XY]|[XY])$/) {
				$autosom_Regions++;
				$Patients{$file}{"tot_Autosomes_depth"} += ${$Regions_r}[$r]{$file}{$RefDepth};
				}
			else {
				if (${$Regions_r}[$r]{"Chrom"} =~ /^(chrX|X)$/) {
					$gonosom_Regions++;
					$Patients{$file}{"tot_chrX_depth"} += ${$Regions_r}[$r]{$file}{$RefDepth};
					}
				else {
					$Patients{$file}{"tot_chrY_depth"} += ${$Regions_r}[$r]{$file}{$RefDepth};
					}
				}
			}
		}
	# Si region mal couverte, on l'imprime dans un fichier POUBELLE
	else
		{ ${$Regions_r}[$r]{"Appel"} = "Couverture_Faible"; }

	}

#my$depthTxt="";
my$sexTxt="";
foreach my$file (@{$Files_r}) {
	##if tot bases
#	if ($RefDepth eq "tot") {		#deprecated: total sequenced bases from bam
#		#autosomes
#		my$cmd = "samtools view -F 0x4";
#		if ($mmq) { $cmd .= " -q $mmq"; }
#		$cmd .= " $file";
#		foreach my$chr (@{ ${$fai_r}{$file} }) {
#			unless ($chr =~ /^(chr[XY]|[XY])$/)
#				{ $cmd .= " $chr"; }
#			}
#		$cmd .= " | perl -ne 'END { print \"\$cnt\" } \@fields = split (/\\t/,\$_); \@bases = split (//,\$fields[9]);";
#		if ($mbq) { $cmd .= " \@quals = map { unpack(\"C*\", \$_ )} split(//,\$fields[10]); for \$i (0..\$#bases) { if ((\$bases[\$i] =~ /[ACGT]/i) && (\$quals[\$i] >= $mbq)) { \$cnt++; } }'"; }
#		else { $cmd .= " for \$i (0..\$#bases) { if (\$bases[\$i] =~ /[ACGT]/i) { \$cnt++; } }'"; }
#		print "$cmd\n";
#		$Patients{$file}{"tot_Autosomes_depth"} = `$cmd`;
#		chomp $Patients{$file}{"tot_Autosomes_depth"};
#		#chrX
#		$cmd = "samtools view -F 0x4";
#		if ($mmq) { $cmd .= " -q $mmq"; }
#		$cmd .= " $file";
#		foreach my$chr (@{ ${$fai_r}{$file} }) {
#			if ($chr =~ /^(chrX|X)$/)
#				{ $cmd .= " $chr"; }
#			}
#		$cmd .= " | perl -ne 'END { print \"\$cnt\" } \@fields = split (/\\t/,\$_); \@bases = split (//,\$fields[9]);";
#		if ($mbq) { $cmd .= " \@quals = map { unpack(\"C*\", \$_ )} split(//,\$fields[10]); for \$i (0..\$#bases) { if ((\$bases[\$i] =~ /[ACGT]/i) && (\$quals[\$i] >= $mbq)) { \$cnt++; } }'"; }
#		else { $cmd .= " for \$i (0..\$#bases) { if (\$bases[\$i] =~ /[ACGT]/i) { \$cnt++; } }'"; }
#		print "$cmd\n";
#		$Patients{$file}{"tot_chrX_depth"} = `$cmd`;
#		chomp $Patients{$file}{"tot_chrX_depth"};
#		#chrY
#		$cmd = "samtools view -F 0x4";
#		if ($mmq) { $cmd .= " -q $mmq"; }
#		$cmd .= " $file";
#		foreach my$chr (@{ ${$fai_r}{$file} }) {
#			if ($chr =~ /^(chrY|Y)$/)
#				{ $cmd .= " $chr"; }
#			}
#		$cmd .= " | perl -ne 'END { print \"\$cnt\" } \@fields = split (/\\t/,\$_); \@bases = split (//,\$fields[9]);";
#		if ($mbq) { $cmd .= " \@quals = map { unpack(\"C*\", \$_ )} split(//,\$fields[10]); for \$i (0..\$#bases) { if ((\$bases[\$i] =~ /[ACGT]/i) && (\$quals[\$i] >= $mbq)) { \$cnt++; } }'"; }
#		else { $cmd .= " for \$i (0..\$#bases) { if (\$bases[\$i] =~ /[ACGT]/i) { \$cnt++; } }'"; }
#		print "$cmd\n";
#		$Patients{$file}{"tot_chrY_depth"} = `$cmd`;
#		chomp $Patients{$file}{"tot_chrY_depth"};
#		
#		$depthTxt .= ${$sampleName_r}{$file}.":\n\ttotal sequenced bases for ${$sampleName_r}{$file} : \n\tautoZ:".$Patients{$file}{"tot_Autosomes_depth"}."\n\tchrX:".$Patients{$file}{"tot_chrX_depth"}."\n\tchrY:".$Patients{$file}{"tot_chrY_depth"}."\n";
#		}

	##gender if no sex file
	$sexTxt .= ${$sampleName_r}{$file}.":\n";
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
	if ($RefNoChrY) {
		if ($Patients{$file}{"Sexe"} eq "H") {
			$Patients{$file}{"Ref_depth"} = ($Patients{$file}{"tot_chrX_depth"} * 2) + $Patients{$file}{"tot_Autosomes_depth"}; } 
		else {
			$Patients{$file}{"Ref_depth"} = $Patients{$file}{"tot_chrX_depth"} + $Patients{$file}{"tot_Autosomes_depth"}; }
		}
	else {
		$Patients{$file}{"Ref_depth"} = $Patients{$file}{"tot_chrX_depth"} + $Patients{$file}{"tot_chrY_depth"} + $Patients{$file}{"tot_Autosomes_depth"};
		}

	}
#if ($depthTxt) { print "\nbases counts:\n$depthTxt\n"; }
print "\ngender:\n$sexTxt\n";

my$meanRef;
foreach my$file (@{$Files_r}) {
	if ($Patients{$file}{"Ref_depth"})
		{ $meanRef += $Patients{$file}{"Ref_depth"}; }
	else { $Patients{$file}{"ecarte"} = 1; }
	}
$meanRef /= scalar@{$Files_r};
foreach my$file (@{$Files_r}) {
	$Patients{$file}{"Ref_depth"} /= $meanRef; 
	print "normalization factor for ${$sampleName_r}{$file} : ".$Patients{$file}{"Ref_depth"}."\n";
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
	foreach my$file (@{$Files_r})
		{ print SORTIE ${$sampleName_r}{$file}."\t"; }
	print SORTIE "Statut_Region\n";


	# PREMIER PARCOURS BIS #
	# RECALCUL DE LA PROFONDEUR TOTALE POUR LA PONDERATION INTRA
	# CE (RE)CALCUL A LIEU SI DES REGIONS ETAIENT "MOCHES" APRES L'APPEL DE CNV
	if ($nb_parcours > 0) {

		# REINITIALISATION DE LA PROFONDEUR TOTALE PAR PATIENT
		foreach my$file (@{$Files_r}) {
			$Patients{$file}{"tot_Autosomes_depth"} = 0;
			$Patients{$file}{"Ref_Gonosomes_Patient"} = 0;
			}
		for (my $r = 0 ; $r < scalar@{$Regions_r} ; $r++) {
			if(!defined ${$Regions_r}[$r]{"Appel"}) {
				foreach my$file (@{$Files_r}) {
					if (${$Regions_r}[$r]{"Chrom"} !~ /^(chr[XY]|[XY])$/) {
						$Patients{$file}{"tot_Autosomes_depth"} += ${$Regions_r}[$r]{$file}{$RefDepth};
						}
					elsif (${$Regions_r}[$r]{"Chrom"} =~ /^(chrX|X)$/) {
						$Patients{$file}{"Ref_Gonosomes_Patient"} += ${$Regions_r}[$r]{$file}{$RefDepth};
						}
					else {
						if ((!$RefNoChrY) && ($Patients{$file}{"Sexe"} eq "H")) {
							$Patients{$file}{"Ref_Gonosomes_Patient"} += ${$Regions_r}[$r]{$file}{$RefDepth};
							}
						}
					}
				}
			}
		
		foreach my$file (@{$Files_r}) {
			if ($RefNoChrY && ($Patients{$file}{"Sexe"} eq "H")) {
				$Patients{$file}{"Ref_depth"} = ($Patients{$file}{"Ref_Gonosomes_Patient"} * 2) + $Patients{$file}{"tot_Autosomes_depth"}; } 
			else {
				$Patients{$file}{"Ref_depth"} = $Patients{$file}{"Ref_Gonosomes_Patient"} + $Patients{$file}{"tot_Autosomes_depth"}; }
			}
		foreach my$file (@{$Files_r})
			{ $meanRef += $Patients{$file}{"Ref_depth"}; }
		$meanRef /= scalar@{$Files_r};
		foreach my$file (@{$Files_r}) {
			$Patients{$file}{"Ref_depth"} /= $meanRef;
			}
		}
	print LOG "sample normalization factors :\n";
	foreach my$file (@{$Files_r}) {
		print LOG "\t${$sampleName_r}{$file} : $Patients{$file}{Ref_depth}\n";
	}
	print LOG "\n";


	# SECOND PARCOURS DES REGIONS #
	# PERMET DE PONDERER LA PROFONDEUR PAR LES AUTRES REGIONS (INTRA) ET ENTRE LES PATIENTS (INTER)
	for (my $r = 0 ; $r < scalar@{$Regions_r} ; $r++) {

		print SORTIE ${$Regions_r}[$r]{"allLine"}."\t";
		print SORTIE $r."\t";

		# SI LA REGION N'EST PAS "MOCHE"
		if (!defined ${$Regions_r}[$r]{"Appel"}) {

			$nb_regions_conservees++;
			${$Regions_r}[$r]{"nb_CNV"}{"DEL"} = 0;
			${$Regions_r}[$r]{"nb_CNV"}{"DUP"} = 0;
			my $nb_evts = 0;
			my$ratio2center = "";

			if ($ratioByGender) {

				# Pour chaque patient, un premier parcours permet la pondération intra-patient
				# (en divisant la profondeur de la region par la profondeur totale du patient)
				@{ ${$Regions_r}[$r]{"all_normByS_depths_fem"} } = ();
				@{ ${$Regions_r}[$r]{"all_normByS_depths_males"} } = ();
				foreach my$file (@{$Files_r}) {
					unless ($Patients{$file}{"ecarte"}) {
						if(${$Regions_r}[$r]{$file}{"mean"} > 0) {
							# seules les prof > 0 prises en compte pour la norm
							${$Regions_r}[$r]{$file}{"normByS_depth"} = ${$Regions_r}[$r]{$file}{"mean"} / $Patients{$file}{"Ref_depth"};
							if ($Patients{$file}{"Sexe"} eq "F") {
								push( @{ ${$Regions_r}[$r]{"all_normByS_depths_fem"} }, ${$Regions_r}[$r]{$file}{"normByS_depth"} );
							} else {
								push( @{ ${$Regions_r}[$r]{"all_normByS_depths_males"} }, ${$Regions_r}[$r]{$file}{"normByS_depth"} );
							}
						} else { ${$Regions_r}[$r]{$file}{"normByS_depth"} = 0; }
					}
				}

				# Nous calculons pour chaque region et chaque sexe la moyenne de profondeur (ponderee pour chaque patient par rapport aux autres regions)

				my@Autosomes=();
				if ($ratioByGender eq "all") {
					if(@{ ${$Regions_r}[$r]{"all_normByS_depths_fem"} }) {
						${$Regions_r}[$r]{"normByR_depth_fem"} = getRegionStats($sortDepths,$CNV_opt_r,\%cnvStats,\@{ ${$Regions_r}[$r]{"all_normByS_depths_fem"} });
					}
					if(@{ ${$Regions_r}[$r]{"all_normByS_depths_males"} }) {
						${$Regions_r}[$r]{"normByR_depth_males"} = getRegionStats($sortDepths,$CNV_opt_r,\%cnvStats,\@{ ${$Regions_r}[$r]{"all_normByS_depths_males"} });
					}
				} else {
					if (${$Regions_r}[$r]{"Chrom"} !~ /^(chr[XY]|[XY])$/) {
						if(@{ ${$Regions_r}[$r]{"all_normByS_depths_fem"} }) {
							push(@Autosomes, @{ ${$Regions_r}[$r]{"all_normByS_depths_fem"} });
						}
						if(@{ ${$Regions_r}[$r]{"all_normByS_depths_males"} }) {
							push(@Autosomes, @{ ${$Regions_r}[$r]{"all_normByS_depths_males"} });
						}
						if (@Autosomes) {
							${$Regions_r}[$r]{"normByR_depth"} = getRegionStats($sortDepths,$CNV_opt_r,\%cnvStats,\@Autosomes);
						}
					} else {
						if(@{ ${$Regions_r}[$r]{"all_normByS_depths_fem"} }) {
							${$Regions_r}[$r]{"normByR_depth_fem"} = getRegionStats($sortDepths,$CNV_opt_r,\%cnvStats,\@{ ${$Regions_r}[$r]{"all_normByS_depths_fem"} });
						}
						if(@{ ${$Regions_r}[$r]{"all_normByS_depths_males"} }) {
							${$Regions_r}[$r]{"normByR_depth_males"} = getRegionStats($sortDepths,$CNV_opt_r,\%cnvStats,\@{ ${$Regions_r}[$r]{"all_normByS_depths_males"} });
						}
					}
				}

			
				# Un second parcours permet l'appel de CNV (en comparant la profondeur pour un patient par la profondeur moyenne de la region)
				foreach my$file (@{$Files_r}) {

					if ($Patients{$file}{"ecarte"}) {
						print SORTIE "NA\t";

					} else {
						# calcul ratio sur valeur avg/med
						my($ratio2center,$ratio2spread);
						if ($ratioByGender eq "all") {
							if($Patients{$file}{"Sexe"} eq "F") {
								if (${$Regions_r}[$r]{"Chrom"} !~ m/^(chrY|Y)$/) {
									if (@{ ${$Regions_r}[$r]{"all_normByS_depths_fem"} }) {	
										($ratio2center,$ratio2spread) = depth_ratio($center,$spread,${$Regions_r}[$r]{$file}{"normByS_depth"},${$Regions_r}[$r]{"normByR_depth_fem"});
									} else {
										${$Regions_r}[$r]{"Appel"} = "No_Data";
										print SORTIE "NA\t";
									}
								} else {
									# les femmes ne sont pas considérées pour l'appel de CNV des régions du chrY
									print SORTIE "NA\t";
								}
							} else {
								if (@{ ${$Regions_r}[$r]{"all_normByS_depths_males"} }) {
									($ratio2center,$ratio2spread) = depth_ratio($center,$spread,${$Regions_r}[$r]{$file}{"normByS_depth"},${$Regions_r}[$r]{"normByR_depth_males"});
								} else {
									${$Regions_r}[$r]{"Appel"} = "No_Data";
									print SORTIE "NA\t";
								}
							}

						} else {
							# Pour les autosomes
							if (${$Regions_r}[$r]{"Chrom"} !~ m/^chr[XY]$|^[XY]$/) {
								if (@Autosomes) {
									($ratio2center,$ratio2spread) = depth_ratio($center,$spread,${$Regions_r}[$r]{$file}{"normByS_depth"},${$Regions_r}[$r]{"normByR_depth"});
								} else {
									${$Regions_r}[$r]{"Appel"} = "No_Data";
									print SORTIE "NA\t";
								}
							} else {	# Pour les gonosomes
								if($Patients{$file}{"Sexe"} eq "F") {
									if (${$Regions_r}[$r]{"Chrom"} =~ m/^(chrX|X)$/) {
										if (@{ ${$Regions_r}[$r]{"all_normByS_depths_fem"} }) {	
											($ratio2center,$ratio2spread) = depth_ratio($center,$spread,${$Regions_r}[$r]{$file}{"normByS_depth"},${$Regions_r}[$r]{"normByR_depth_fem"});
										} else {
											${$Regions_r}[$r]{"Appel"} = "No_Data";
											print SORTIE "NA\t";
										}
									} else {
										# les femmes ne sont pas considérées pour l'appel de CNV des régions du chrY
										print SORTIE "NA\t";
									}
								} else {
									if (@{ ${$Regions_r}[$r]{"all_normByS_depths_males"} }) {
										($ratio2center,$ratio2spread) = depth_ratio($center,$spread,${$Regions_r}[$r]{$file}{"normByS_depth"},${$Regions_r}[$r]{"normByR_depth_males"});
									} else {
										${$Regions_r}[$r]{"Appel"} = "No_Data";
										print SORTIE "NA\t";
									}
								}

							}
						}
						if (defined $ratio2center) { 
							${$Regions_r}[$r]{$file}{"ratio2center"} = $ratio2center;
							if (defined	$ratio2spread) {
								${$Regions_r}[$r]{$file}{"ratio2spread"} = $ratio2spread;
								print SORTIE sprintf("%.3f",$ratio2center)."(".sprintf("%.3f",$ratio2spread).")\t";
							} else {
								print SORTIE sprintf("%.3f",$ratio2center)."\t";
							}
						}

						## test CNV
						my($del,$dup) = CNV_test($center_test,$spread_test,$ratio2center,$ratio2spread,$seuil_del,$seuil_dup,$spread_del,$spread_dup);
						if ($del) {
							${$Regions_r}[$r]{"nb_CNV"}{"DEL"}++;
							$nb_evts++;
							$CNV{$file}++;
							$Results{$file}{$r} = "DEL";
						} elsif ($dup) {
							${$Regions_r}[$r]{"nb_CNV"}{"DUP"}++;
							$nb_evts++;
							$CNV{$file}++;
							$Results{$file}{$r} = "DUP";
						}
					}
				}

				my$recurrent="";
				if (@{ ${$Regions_r}[$r]{"all_normByS_depths_fem"} } && @{ ${$Regions_r}[$r]{"all_normByS_depths_males"} }) {
					if ( ($nb_evts / (scalar@{ ${$Regions_r}[$r]{"all_normByS_depths_fem"} } + scalar@{ ${$Regions_r}[$r]{"all_normByS_depths_males"} })) > $seuil_region && $nb_parcours > 0 ) { $recurrent = 1; }
				} elsif (@{ ${$Regions_r}[$r]{"all_normByS_depths_fem"} } && !@{ ${$Regions_r}[$r]{"all_normByS_depths_males"} }) {
					if ( ($nb_evts / scalar@{ ${$Regions_r}[$r]{"all_normByS_depths_fem"} }) > $seuil_region && $nb_parcours > 0 ) { $recurrent = 1; }
				} elsif (!@{ ${$Regions_r}[$r]{"all_normByS_depths_fem"} } && @{ ${$Regions_r}[$r]{"all_normByS_depths_males"} }) {
					if ( ($nb_evts / scalar@{ ${$Regions_r}[$r]{"all_normByS_depths_males"} }) > $seuil_region && $nb_parcours > 0 ) { $recurrent = 1; }
				}
				if ($recurrent) {
					print SORTIE "CNV_Recurrent";
					print LOG "Region ecartee \: ". ${$Regions_r}[$r]{"allLine"}."\n";
					${$Regions_r}[$r]{"Appel"} = "CNV_Recurrent";
					$regions_ecartees++;
					$continuer = 1;

				} else {		
					print SORTIE "OK";
				}
				print SORTIE "\n";


			} else {		#ratio with both genders whatever the chrom

				# Pour chaque patient, un premier parcours permet la pondération intra-patient
				# (en divisant la profondeur de la region par la profondeur totale recalculée du patient)
				@{ ${$Regions_r}[$r]{"all_normByS_depths"} } = ();
				foreach my$file (@{$Files_r}) {
					unless ($Patients{$file}{"ecarte"}) {
						if(${$Regions_r}[$r]{$file}{"mean"} > 0) {
							# seules les prof > 0 prises en compte pour la norm
							${$Regions_r}[$r]{$file}{"normByS_depth"} = ${$Regions_r}[$r]{$file}{"mean"} / $Patients{$file}{"Ref_depth"};
							if ( ($Patients{$file}{"Sexe"} eq "H") && (${$Regions_r}[$r]{"Chrom"} =~ /^(chrX|X)$/) )
								{ push(@{ ${$Regions_r}[$r]{"all_normByS_depths"} }, (${$Regions_r}[$r]{$file}{"normByS_depth"}*2)); }
							elsif ( ($Patients{$file}{"Sexe"} eq "F") && (${$Regions_r}[$r]{"Chrom"} =~ /^(chrY|Y)$/) )
								{ next ; }
							else { push(@{ ${$Regions_r}[$r]{"all_normByS_depths"} }, ${$Regions_r}[$r]{$file}{"normByS_depth"}); }
						} else { ${$Regions_r}[$r]{$file}{"normByS_depth"} = 0; }
					}
				}

				# Nous calculons pour chaque region la mediane/moyenne de profondeur (ponderee pour chaque patient par rapport aux autres regions)
				if(@{ ${$Regions_r}[$r]{"all_normByS_depths"} }) {
					${$Regions_r}[$r]{"normByR_depth"} = getRegionStats($sortDepths,$CNV_opt_r,\%cnvStats,\@{ ${$Regions_r}[$r]{"all_normByS_depths"} });
				}

				# Un second parcours permet l'appel de CNV (en comparant la profondeur pour un patient par la profondeur moyenne de la region)
				foreach my$file (@{$Files_r}) {

					if ($Patients{$file}{"ecarte"}) {
						print SORTIE "NA\t";

					} else {

						if (@{ ${$Regions_r}[$r]{"all_normByS_depths"} }) {

							# calcul ratio sur valeur moy/med
							my($ratio2center,$ratio2spread);
							if($Patients{$file}{"Sexe"} eq "F") {
								if (${$Regions_r}[$r]{"Chrom"} !~ m/^(chrY|Y)$/) {
									($ratio2center,$ratio2spread) = depth_ratio($center,$spread,${$Regions_r}[$r]{$file}{"normByS_depth"},${$Regions_r}[$r]{"normByR_depth"});
								} else {
									#les femmes ne sont pas considérées pour l'appel de CNV des régions du chrY
									print SORTIE "NA\t";
								}
							} else {
								#*2 for chrX
								if (${$Regions_r}[$r]{"Chrom"} =~ m/^(chrX|X)$/) {
									($ratio2center,$ratio2spread) = depth_ratio($center,$spread,(2*${$Regions_r}[$r]{$file}{"normByS_depth"}),${$Regions_r}[$r]{"normByR_depth"});
								} else {
									($ratio2center,$ratio2spread) = depth_ratio($center,$spread,${$Regions_r}[$r]{$file}{"normByS_depth"},${$Regions_r}[$r]{"normByR_depth"});
								}
							}
							if (defined $ratio2center) { 
								${$Regions_r}[$r]{$file}{"ratio2center"} = $ratio2center;
								if (defined	$ratio2spread) {
									${$Regions_r}[$r]{$file}{"ratio2spread"} = $ratio2spread;
									print SORTIE sprintf("%.3f",$ratio2center)."(".sprintf("%.3f",$ratio2spread).")\t";
								} else {
									print SORTIE sprintf("%.3f",$ratio2center)."\t";
								}
							}

							## test CNV
							my($del,$dup) = CNV_test($center_test,$spread_test,$ratio2center,$ratio2spread,$seuil_del,$seuil_dup,$spread_del,$spread_dup);
							if ($del) {
								${$Regions_r}[$r]{"nb_CNV"}{"DEL"}++;
								$nb_evts++;
								$CNV{$file}++;
								$Results{$file}{$r} = "DEL";
							} elsif ($dup) {
								${$Regions_r}[$r]{"nb_CNV"}{"DUP"}++;
								$nb_evts++;
								$CNV{$file}++;
								$Results{$file}{$r} = "DUP";
							}

						} else {
							${$Regions_r}[$r]{"Appel"} = "No_Data";
							print SORTIE "NA\t";
						}

					}

				}

				if ( @{ ${$Regions_r}[$r]{"all_normByS_depths"} } && ($nb_evts/scalar@{ ${$Regions_r}[$r]{"all_normByS_depths"} }) > $seuil_region && $nb_parcours > 0 ) {
					print SORTIE "CNV_Recurrent";
					print LOG "Region ecartee \: ". ${$Regions_r}[$r]{"allLine"}."\n";
					${$Regions_r}[$r]{"Appel"} = "CNV_Recurrent";
					$regions_ecartees++;
					$continuer = 1;

				} else {		
					print SORTIE "OK";
				}
				print SORTIE "\n";
			}

		# SI LA REGION EST "MOCHE"
		} else {
			foreach my$file (@{$Files_r}) {
				print SORTIE "NA\t";
			}
			print SORTIE ${$Regions_r}[$r]{"Appel"}."\n";
		}

	}

	print LOG "Nombre de regions ecartees lors de cette iteration \: ".$regions_ecartees."\n\n\n";	

	my $patients_ecartes = 0;

	foreach my$file (@{$Files_r}) {

		if (! $Patients{$file}{"ecarte"}) {
			my $prcent_evt;
			if (defined($CNV{$file})) {

				$prcent_evt = $CNV{$file}/$nb_regions_conservees;

				if ($prcent_evt >= $seuil_patient) {
					print LOG "Patient ecarte \: ".${$sampleName_r}{$file}."\n";
					$Patients{$file}{"ecarte"} = 1;
					$continuer = 1;
					$patients_ecartes++;			
				} else {
					$Patients{$file}{"ecarte"} = 0;
				}


			} else {
				print ${$sampleName_r}{$file}." \: aucun CNV identifie\n";
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


for my$file (@{$Files_r}) {
	mkdir "$outdir/${$sampleName_r}{$file}";
	}


##print CNV foreach sample

## $Results{$file}{$r} = "DEL";
## @{ $RegionOrder{$Chrom} } = [ regions sorted by pos ]
my(%uniqReg,%RegInArray,%regionOrder,%regionIndice);
for (my $r = 0 ; $r < scalar@{$Regions_r} ; $r++) {
	unless (exists $uniqReg{ ${$Regions_r}[$r]{"Chrom"}."-".${$Regions_r}[$r]{"Start"}."-".${$Regions_r}[$r]{"End"} }) {
		push (@{ $RegInArray{ ${$Regions_r}[$r]{"Chrom"} } }, $r);
		$uniqReg{ ${$Regions_r}[$r]{"Chrom"}."-".${$Regions_r}[$r]{"Start"}."-".${$Regions_r}[$r]{"End"} } = 1;
		}
	my@tab = split(/\t/, ${$Regions_r}[$r]{"allLine"});
	if ($tab[3]) {
		my@tab2 = split(/:|,/, $tab[3]);
		${$Regions_r}[$r]{"label"} = substr($tab2[0], 0, 25);
		${$Regions_r}[$r]{"Gene"} = $tab2[0];
		}
	else { ${$Regions_r}[$r]{"label"} = "$tab[1]-$tab[2]"; }
	}
foreach my$Chrom (keys%RegInArray) {
	##sort by region ends (in case several regions with same start) then by starts
	@{ $regionOrder{$Chrom} } = sort{${$Regions_r}[$a]{"End"}<=>${$Regions_r}[$b]{"End"}}@{ $RegInArray{$Chrom} }; 
	@{ $regionOrder{$Chrom} } = sort{${$Regions_r}[$a]{"Start"}<=>${$Regions_r}[$b]{"Start"}}@{ $RegInArray{$Chrom} }; 
	for (my$i=0;$i<scalar@{ $regionOrder{$Chrom} };$i++) { $regionIndice{ $regionOrder{$Chrom}[$i] } = $i; }
	}

## @{ $orderedCNV{$patient}{$Chrom} } = [r1, r2,...]
my%orderedCNV;
foreach my$file (keys%Results) {
	my(%uniqCNV,%CNVinArray);
	foreach my$r (sort{$a<=>$b}keys%{ $Results{$file} }) {
		if (!exists $uniqCNV{${$Regions_r}[$r]{"Chrom"}."-".${$Regions_r}[$r]{"Start"}."-".${$Regions_r}[$r]{"End"}}) { 
			push (@{ $CNVinArray{ ${$Regions_r}[$r]{"Chrom"} } }, $r);
			$uniqCNV{ ${$Regions_r}[$r]{"Chrom"}."-".${$Regions_r}[$r]{"Start"}."-".${$Regions_r}[$r]{"End"} } = 1;
			}
		}
	foreach my$Chrom (keys%CNVinArray) {
		##sort by region ends (in case several regions with same start) then by starts
		@{ $orderedCNV{$file}{$Chrom} } = sort{${$Regions_r}[$a]{"End"}<=>${$Regions_r}[$b]{"End"}}@{ $CNVinArray{$Chrom} }; 
		@{ $orderedCNV{$file}{$Chrom} } = sort{${$Regions_r}[$a]{"Start"}<=>${$Regions_r}[$b]{"Start"}}@{ $CNVinArray{$Chrom} }; 
		}
	}

unless ($minCNV) { $minCNV = 1; }
my(%Result2,%Result3,%Result4);
foreach my$file (keys%Results) {
	open (CNV1,">$outdir/${$sampleName_r}{$file}/CNV_${$sampleName_r}{$file}.summary.txt") or die "Pb lors de l'ecriture du fichier sortie $!\n";
	print CNV1 "#Region\tCNV\tclean_intervals_Nbr\tclean_ratio_to_$center\tdirty_intervals_Nbr\tdirty_ratio_to_$center\toverlapping_Samples\n";
	open (CNV2,">$outdir/${$sampleName_r}{$file}/CNV_${$sampleName_r}{$file}.allIntervals.txt") or die "Pb lors de l'ecriture du fichier sortie $!\n";
	print CNV2 "#Chrom\tStart\tEnd\tLength\tInfo\tinterval_Order\tCNV_type\tratio_to_$center\toccurences";
	if (@cnvFields) {
		foreach (@cnvFields) { 
			#if ($_ eq "norm") {
			#	if ($norm eq "std") { print CNV2 "\tavg\tstd"; }
			#	else  { print CNV2 "\t$norm"; }
			#	}
			#else { print CNV2 "\t$_"; }
			print CNV2 "\t$_";
			}
		}
	print CNV2 "\n";
	foreach my$Chromosome (@{$ChromOrder_r}) {
		my$Chrom = $Chromosome ; $Chrom =~ s/chr//i;
		if (exists $orderedCNV{$file}{$Chrom}) {
			my$r=0;		##index in @{ $orderedCNV{$file}{$Chrom} }
			while ($r < scalar@{ $orderedCNV{$file}{$Chrom} } ) {
				##merge consecutive CNVs
				my($i,$cnvOK,$nonCNVtot,$nextReg_r) = mergeConsecutiveCNV("",$maxNonCNV,$orderedCNV{$file}{$Chrom}[$r],$regionIndice{$orderedCNV{$file}{$Chrom}[$r]},\@{ $regionOrder{$Chrom} },\%{ $Results{$file} },$Regions_r);
				if ($maxNonCNVrate) {
					while (($nonCNVtot/($cnvOK+1)) > $maxNonCNVrate) {
						($i,$cnvOK,$nonCNVtot,$nextReg_r) = mergeConsecutiveCNV(($cnvOK-1),$maxNonCNV,$orderedCNV{$file}{$Chrom}[$r],$regionIndice{$orderedCNV{$file}{$Chrom}[$r]},\@{ $regionOrder{$Chrom} },\%{ $Results{$file} },$Regions_r);
						}
					}
				##overlapping samples
				my$overlapCNT=0; my%overlapSMPL=();
				for (my$j=0;$j<=$cnvOK;$j++) {
					foreach my$other (keys%Patients) {
						if ( (exists $Results{$other}{${$nextReg_r}[$j]}) && ($Results{$other}{${$nextReg_r}[$j]} eq $Results{$file}{${$nextReg_r}[0]}) && (!exists $overlapSMPL{$other}) ) { 
							$overlapCNT++;
							$overlapSMPL{$other}=1;
							}
						}
					}
				##print patient's summary and allIntervals
				if ($cnvOK >= ($minCNV-1)) {
					if ($cnvOK == 0) {
						$Result2{$file}{${$nextReg_r}[0]} = $Results{$file}{${$nextReg_r}[0]};
						$Result3{$file}{${$nextReg_r}[0]}{"type"} = $Results{$file}{${$nextReg_r}[0]};
						$Result3{$file}{${$nextReg_r}[0]}{"end"} = ${$nextReg_r}[0];
						$Result4{$file}{$Chrom}{${$Regions_r}[${$nextReg_r}[0]]{"Start"}}{${$Regions_r}[${$nextReg_r}[0]]{"End"}} = $Results{$file}{${$nextReg_r}[0]};
						if (exists $regionIndice{${$nextReg_r}[0]}) {
							print CNV1 $Chromosome.":".${$Regions_r}[${$nextReg_r}[0]]{"Start"}."-".${$Regions_r}[${$nextReg_r}[0]]{"End"}."\t".$Results{$file}{${$nextReg_r}[0]}."\t1\t".sprintf("%.3f",${$Regions_r}[${$nextReg_r}[0]]{$file}{"ratio2center"})."\t.\t.\t$overlapCNT\n";
							print CNV2 $Chromosome."\t".${$Regions_r}[${$nextReg_r}[0]]{"Start"}."\t".${$Regions_r}[${$nextReg_r}[0]]{"End"}."\t".(${$Regions_r}[${$nextReg_r}[0]]{"End"}-${$Regions_r}[${$nextReg_r}[0]]{"Start"}+1)." bp\t".${$Regions_r}[${$nextReg_r}[0]]{"label"}."\t".($regionIndice{${$nextReg_r}[0]}+1)."\t".$Results{$file}{${$nextReg_r}[0]}."\t".sprintf("%.3f",${$Regions_r}[${$nextReg_r}[0]]{$file}{"ratio2center"})."\t".${$Regions_r}[${$nextReg_r}[0]]{"nb_CNV"}{$Results{$file}{${$nextReg_r}[0]}};
							if (@cnvFields) {
								my$txt = printCNVfields($ratioByGender,\@cnvFields,\%{ ${$Regions_r}[${$nextReg_r}[0]] },\%{ $Patients{$file} });
								print CNV2 "$txt";
								}
							print CNV2 "\n\n";
							}
						}
					else {
						my@cleanCNV = (); my@dirtyCNV = ();
						for (my$j=0;$j<=$cnvOK;$j++) {
							if (exists $Results{$file}{${$nextReg_r}[$j]}) {
								$Result2{$file}{${$nextReg_r}[$j]} = $Results{$file}{${$nextReg_r}[$j]};
								$Result4{$file}{$Chrom}{${$Regions_r}[${$nextReg_r}[$j]]{"Start"}}{${$Regions_r}[${$nextReg_r}[$j]]{"End"}} = $Results{$file}{${$nextReg_r}[$j]};
								push(@cleanCNV, ${$Regions_r}[${$nextReg_r}[$j]]{$file}{"ratio2center"});
								push(@dirtyCNV, ${$Regions_r}[${$nextReg_r}[$j]]{$file}{"ratio2center"});
								if (exists $regionIndice{${$nextReg_r}[$j]}) {
									print CNV2 $Chromosome."\t".${$Regions_r}[${$nextReg_r}[$j]]{"Start"}."\t".${$Regions_r}[${$nextReg_r}[$j]]{"End"}."\t".(${$Regions_r}[${$nextReg_r}[$j]]{"End"}-${$Regions_r}[${$nextReg_r}[$j]]{"Start"}+1)." bp\t".${$Regions_r}[${$nextReg_r}[$j]]{"label"}."\t".($regionIndice{${$nextReg_r}[$j]}+1)."\t".$Results{$file}{${$nextReg_r}[$j]}."\t".sprintf("%.3f",${$Regions_r}[${$nextReg_r}[$j]]{$file}{"ratio2center"})."\t".${$Regions_r}[${$nextReg_r}[$j]]{"nb_CNV"}{$Results{$file}{${$nextReg_r}[$j]}};
									if (@cnvFields) {
										my$txt = printCNVfields($ratioByGender,\@cnvFields,\%{ ${$Regions_r}[${$nextReg_r}[$j]] },\%{ $Patients{$file} });
										print CNV2 "$txt";
										}
									print CNV2 "\n";
									}
								}
							else {
								$Result2{$file}{${$nextReg_r}[$j]} = "NA";
								$Result4{$file}{$Chrom}{${$Regions_r}[${$nextReg_r}[$j]]{"Start"}}{${$Regions_r}[${$nextReg_r}[$j]]{"End"}} = "NA";
								if (exists $regionIndice{${$nextReg_r}[$j]}) {
									print CNV2 $Chromosome."\t".${$Regions_r}[${$nextReg_r}[$j]]{"Start"}."\t".${$Regions_r}[${$nextReg_r}[$j]]{"End"}."\t".(${$Regions_r}[${$nextReg_r}[$j]]{"End"}-${$Regions_r}[${$nextReg_r}[$j]]{"Start"}+1)." bp\t".${$Regions_r}[${$nextReg_r}[$j]]{"label"}."\t".($regionIndice{${$nextReg_r}[$j]}+1)."\t";
									if (exists ${$Regions_r}[${$nextReg_r}[$j]]{"Appel"}) { print CNV2 "NA\t"; }
									else { print CNV2 "no\t"; }
									if (exists ${$Regions_r}[${$nextReg_r}[$j]]{$file}{"ratio2center"}) {
										print CNV2 sprintf("%.3f",${$Regions_r}[${$nextReg_r}[$j]]{$file}{"ratio2center"})."\t".${$Regions_r}[${$nextReg_r}[$j]]{"nb_CNV"}{$Results{$file}{${$nextReg_r}[0]}};
										push(@dirtyCNV, ${$Regions_r}[${$nextReg_r}[$j]]{$file}{"ratio2center"});
										}
									else { print CNV2 "na\tna"; }
									if (@cnvFields) {
										my$txt = printCNVfields($ratioByGender,\@cnvFields,\%{ ${$Regions_r}[${$nextReg_r}[$j]] },\%{ $Patients{$file} });
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
						print CNV1 $Chromosome.":".${$Regions_r}[${$nextReg_r}[0]]{"Start"}."-".${$Regions_r}[${$nextReg_r}[$cnvOK]]{"End"}."\t".$Results{$file}{${$nextReg_r}[0]}."\t".scalar@cleanCNV."\t".sprintf("%.3f",$cleanAverage)."\t".($cnvOK+1)."\t".sprintf("%.3f",$dirtyAverage)."\t$overlapCNT\n";
						$Result3{$file}{${$nextReg_r}[0]}{"type"} = $Results{$file}{${$nextReg_r}[0]};
						$Result3{$file}{${$nextReg_r}[0]}{"end"} = ${$nextReg_r}[$cnvOK];
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


##summary
open (OUT,">$outdir/CNV.summary.txt") or die "Pb lors de l'ecriture du fichier sortie $!\n";
print OUT "CNV analysis\n";
print OUT "\nparameters:\n";
if ($center_test) {
	print OUT "\tlevel = ratio to $center of all normalized sample depths
	intervals detected as deletion if normalized depth below $seuil_del of $center
	intervals detected as duplication if normalized depth above $seuil_dup of $center\n";
	}
if ($spread_test) {
	print OUT "\tdispersion = $spread of normalized sample depths
	intervals detected as deletion if normalized depth below $center + ($spread_del * $spread)
	intervals detected as duplication if normalized depth above $center + ($spread_dup * $spread)\n";
	}
if ($seuil_patient <1) { print OUT "\tsamples kept if NO more than ".(100*$seuil_patient)."% of intervals are CNVs\n"; }
if ($seuil_region <1) { print OUT "\tintervals kept if NO more than ".(100*$seuil_region)."% of samples are CNVs\n"; }
if ($minCov) { print OUT "\tintervals discarded if at least one sample covered less or equal to $minCov\n"; }
if (defined $minDP) { print OUT "\tntervals discarded if avg depth of all samples are less or equal to $minDP\n"; }
if ($minCNV) { print OUT "\tCNV calling requires at least $minCNV consecutive called intervals\n"; }
if ($maxNonCNV) { print OUT "\tCNV calling tolerates NO more than $maxNonCNV non-CNV inside\n"; }
if ($ratioByGender eq "all") { print OUT "\tinterval depth normalization by gender for all chromosomes\n"; }
elsif ($ratioByGender eq "gono") { print OUT "\tinterval depth normalization by gender for sex chromosomes only\n"; }
if ($RefDepth eq "mean") { print OUT "\teach sample normalized by the sum of all region mean depths\n"; }
else { print OUT "\teach sample normalized by the sum of all region total depths\n"; }
if ($RefNoChrY) { print OUT "\t(during sample normalization, depths within chrX are doubled for males, and those within chrY are skipped)\n"; }
else { print OUT "\t(during sample normalization, depths from all chr are taken, whatever the sex)\n"; }

#if ($depthTxt) { print OUT "bases counts:\n$depthTxt\n"; }
if ($sexTxt) { print OUT "\ngender determination:\n$sexTxt\n"; }

print OUT "\nintervals nber: ".scalar@{$Regions_r}."\n";
print OUT "\nintervals discarded:\n";
my$N_noData=0; my$N_CNV_Recurrent=0; my$N_lowCov=0;
for (my $r = 0 ; $r < scalar@{$Regions_r} ; $r++) {
	if(defined ${$Regions_r}[$r]{"Appel"}) {
		if (${$Regions_r}[$r]{"Appel"} eq "No_Data") { $N_noData++; }
		elsif (${$Regions_r}[$r]{"Appel"} eq "CNV_Recurrent") { $N_CNV_Recurrent++; }
		elsif (${$Regions_r}[$r]{"Appel"} eq "Couverture_Faible")  { $N_lowCov++; }
		}
	}
print OUT "\tnot enough data: $N_noData\n";
print OUT "\tCNV_Recurrent: $N_CNV_Recurrent\n";
if (defined $minDP || $minCov) { print OUT "\tlow coverage :  $N_lowCov\n"; }

print OUT "\nsamples discarded:\n";
my$N_discard=0;
foreach my$file (@{$Files_r}) {
		if ($Patients{$file}{"ecarte"}) { print OUT "\t${$sampleName_r}{$file}\n"; $N_discard++; }
		}
unless ($N_discard) { print OUT "\tnone\n"; }

print OUT "\nResults :\n
patient\tCNV\tnber\n";
foreach my$file (@{$Files_r}) {
	unless($Patients{$file}{"ecarte"}) {
		my$N_dup=0; my$N_del=0;
		foreach (keys%{ $Result3{$file} }) {
			if ($Result3{$file}{$_}{"type"} eq "DUP") { $N_dup++; }
			elsif ($Result3{$file}{$_}{"type"} eq "DEL") { $N_del++; }
			}
		print OUT "\n${$sampleName_r}{$file}: 
		DUP: $N_dup
		DEL: $N_del\n";
		}
	}
close OUT;


##print graph foreach Chrom/CNV
if ($graphByChr) {
	print "drawing graphs:\n";
	##graph by Chrom/CNV ?
	my$Nbr_Reg_max = 0;
	foreach my$Chrom (keys%regionOrder) {
		if (scalar@{ $regionOrder{$Chrom} } > $Nbr_Reg_max) {
			$Nbr_Reg_max = scalar@{ $regionOrder{$Chrom} };
			}
		}
	my$graphByCNV = 0;
	if ($Nbr_Reg_max > ${$CNV_opt_r}{"switch2graphByCNV"}) { $graphByCNV = 1; print "\tgraphs by CNV:\n"; }
	else { print "\tgraphs by chrom:\n"; }
	foreach my$file (keys%Result2)  {
		if ($graphByCNV) {
			graphByCNV1("$outdir/".${$sampleName_r}{$file},$CNV_opt_r,$file,$Files_r,$sampleName_r,$Regions_r,$ChromOrder_r,\%regionOrder,\%Patients,\%Result2,\%Result3);
			}
		else {
			graphByChr1($nGraf,"$outdir/".${$sampleName_r}{$file},$CNV_opt_r,$file,$Files_r,$sampleName_r,$Regions_r,$ChromOrder_r,\%regionOrder,\%Patients,\%Result2);
			}
		}
	}

return(\%Patients,\%Result4);

}



####################
####################
####################


#CNV_tool::CNV_reAnalyse($fichier_cov,$outdir,$nGraf,$fichier_sexe,\%CNV_opt);
sub CNV_reAnalyse
{
my($fichier_cov,$outdir,$nGraf,$fichier_sexe,$CNV_opt_r)=@_;

my$center = ${$CNV_opt_r}{"center"};
my$spread = ${$CNV_opt_r}{"spread"};
my$center_test = ${$CNV_opt_r}{"center_test"};
my$spread_test = ${$CNV_opt_r}{"spread_test"};
my$seuil_del = ${$CNV_opt_r}{"seuil_del"};
my$seuil_dup = ${$CNV_opt_r}{"seuil_dup"};
my$spread_del = ${$CNV_opt_r}{"spread_del"};
my$spread_dup = ${$CNV_opt_r}{"spread_dup"};
my$range = ${$CNV_opt_r}{"range"};		##sample within med-$range*quartile1 and med+$range*quartile3 will be included for avg+std calculation
my$seuil_region = ${$CNV_opt_r}{"seuil_region"};
my$seuil_patient = ${$CNV_opt_r}{"seuil_patient"};
my$seuil_cov = ${$CNV_opt_r}{"seuil_cov"};
my$minDP = ${$CNV_opt_r}{"min_DP"} ;
my$minCNV = ${$CNV_opt_r}{"min_following_CNV"};
my$maxNonCNV = ${$CNV_opt_r}{"max_Non_CNV"};
my$maxNonCNVrate = ${$CNV_opt_r}{"max_Non_CNV_rate"};
my$ratioByGender = ${$CNV_opt_r}{"ratioByGender"};
my$RefDepth = ${$CNV_opt_r}{"RefDepth"};
my$RefNoChrY = ${$CNV_opt_r}{"RefNoChrY"};
my$graphByChr = ${$CNV_opt_r}{"chromGraph"};

my@cnvFields = @{ ${$CNV_opt_r}{"fields"} };	#min, max, avg , std , med
my%cnvStats;
foreach (@cnvFields) { $cnvStats{$_} = 1; }
$cnvStats{$center} = 1;
if (defined $spread) { $cnvStats{$spread} = 1; }
if (exists $cnvStats{"std"}) { $cnvStats{"avg"} = 1;  }
if (defined $range || (defined $spread && $spread eq "Qtile") || exists $cnvStats{"Q1"} || exists $cnvStats{"Q3"}) {
	$cnvStats{"med"} = 1;
	$cnvStats{"Qtile"} = 1;
	}
my$sortDepths = 0;
if (exists $cnvStats{"med"} || exists $cnvStats{"min"} || exists $cnvStats{"max"}) { $sortDepths = 1; }

my%Sexe;
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
			if (($p > 0) && (!exists$idxV{($p-1)}{"tot"}) && ($RefDepth eq "tot")) { die "tot column not found in DeCovA output\n"; }
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
				if ($cov <= $seuil_cov) { $region_a_conserver = 0; }
				else { $region_a_conserver = 1; last; }
			}
		}
		if (defined $minDP) {	#def = 0
			my(@allDepth,$centerDepth);
			for (my$p=0;$p<scalar@idxP;$p++) { push(@allDepth, $INFOS[$idxV{$p}{"mean"}]); }
			#$centerDepth = mediane(\@allDepth);
			$centerDepth = average(\@allDepth);
			if ( $centerDepth <= $minDP) { $region_a_conserver = 0; }
			#for (my$p=0;$p<scalar@idxP;$p++) {
			#	if ($INFOS[$idxV{$p}{"mean"}] < $minDP) { $region_a_conserver = 0; }
			#	else { $region_a_conserver = 1; last; }
		}

		# Si la region est mal couverte, on l'imprime dans un fichier POUBELLE
		if ($region_a_conserver != 1) {
			$Regions[$region_nbr]{"Appel"} = "Couverture_Faible";

		} else {
		# Sinon, on parcourt la region une seconde fois pour recuperer les informations désirées pour chacun des patients (et donc la profondeur moyenne pour la region)
			for (my $patient = 0 ; $patient < $patient_nbr ; $patient++) {
				$Regions[$region_nbr]{$patient}{"raw_depth"} = $INFOS[$idxV{$patient}{"mean"}];
				$Regions[$region_nbr]{$patient}{"for_Tot_depth"} = $INFOS[$idxV{$patient}{$RefDepth}];	#"mean" or "tot"
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
	##gender if no sex file
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

	##Ref_Profondeur_Patient
	if ($RefNoChrY) {
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
						if ((!$RefNoChrY) && ($Patients{$patient}{"Sexe"} eq "H")) {
							$Patients{$patient}{"tot_sexChr_depth"} += $Regions[$region_nbr]{$patient}{"for_Tot_depth"};
						}
					}
				}
			}
		}

		for (my $p = 0 ; $p < $patient_nbr ; $p++) {
			if ($RefNoChrY && ($Patients{$p}{"Sexe"} eq "H")) {
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

			if ($ratioByGender) {

				# Pour chaque patient, un premier parcours permet la pondération intra-patient
				# (en divisant la profondeur de la region par la profondeur totale du patient)
				@{ $Regions[$r]{"all_normByS_depths_fem"} } = ();
				@{ $Regions[$r]{"all_normByS_depths_males"} } = ();
				for (my $patient = 0 ; $patient < $patient_nbr ; $patient++) {
					unless ($Patients{$patient}{"ecarte"}) {
						if($Regions[$r]{$patient}{"raw_depth"} > 0) {
							# seules les prof > 0 prises en compte pour la norm
							$Regions[$r]{$patient}{"normByS_depth"} = $Regions[$r]{$patient}{"raw_depth"}/$Patients{$patient}{"Ref_depth"};
							if ($Patients{$patient}{"Sexe"} eq "F") {
								push( @{ $Regions[$r]{"all_normByS_depths_fem"} }, $Regions[$r]{$patient}{"normByS_depth"} );
							} else {
								push( @{ $Regions[$r]{"all_normByS_depths_males"} }, $Regions[$r]{$patient}{"normByS_depth"} );
							}
						} else { $Regions[$r]{$patient}{"normByS_depth"} = 0; }
					}
				}

				# Nous calculons pour chaque region et chaque sexe la moy/mediane de profondeur (ponderee pour chaque patient par rapport aux autres regions)
				my@Autosomes=();
				if ($ratioByGender eq "all") {
					if(@{ $Regions[$r]{"all_normByS_depths_fem"} }) {
						$Regions[$r]{"normByR_depth_fem"} = getRegionStats($sortDepths,$CNV_opt_r,\%cnvStats,\@{ $Regions[$r]{"all_normByS_depths_fem"} });
					}
					if(@{ $Regions[$r]{"all_normByS_depths_males"} }) {
						$Regions[$r]{"normByR_depth_males"} = getRegionStats($sortDepths,$CNV_opt_r,\%cnvStats,\@{ $Regions[$r]{"all_normByS_depths_males"} });
					}
				} else {
					if ($Regions[$r]{"Chrom"} !~ m/^chr[XY]$|^[XY]$/) {
						if(@{ $Regions[$r]{"all_normByS_depths_fem"} }) {
							push(@Autosomes, @{ $Regions[$r]{"all_normByS_depths_fem"} });
						}
						if(@{ $Regions[$r]{"all_normByS_depths_males"} }) {
							push(@Autosomes, @{ $Regions[$r]{"all_normByS_depths_males"} });
						}
						if (@Autosomes) {
							$Regions[$r]{"normByR_depth"} = getRegionStats($sortDepths,$CNV_opt_r,\%cnvStats,\@Autosomes);
						}
					} else {
						if(@{ $Regions[$r]{"all_normByS_depths_fem"} }) {
							$Regions[$r]{"normByR_depth_fem"} = getRegionStats($sortDepths,$CNV_opt_r,\%cnvStats,\@{ $Regions[$r]{"all_normByS_depths_fem"} });
						}
						if(@{ $Regions[$r]{"all_normByS_depths_males"} }) {
							$Regions[$r]{"normByR_depth_males"} = getRegionStats($sortDepths,$CNV_opt_r,\%cnvStats,\@{ $Regions[$r]{"all_normByS_depths_males"} });
						}
					}
				}
				
				# Un second parcours permet l'appel de CNV (en comparant la profondeur pour un patient par la profondeur moyenne de la region)
				for (my $patient = 0 ; $patient < $patient_nbr ; $patient++) {

					if ($Patients{$patient}{"ecarte"}) {
						print SORTIE "NA\t";

					} else {

						# calcul ratio sur valeur avg/med
						my($ratio2center,$ratio2spread);
						if ($ratioByGender eq "all") {
							if($Patients{$patient}{"Sexe"} eq "F") {
								if ($Regions[$r]{"Chrom"} !~ m/^(chrY|Y)$/) {
									if (@{ $Regions[$r]{"all_normByS_depths_fem"} }) {	
										($ratio2center,$ratio2spread) = depth_ratio($center,$spread,$Regions[$r]{$patient}{"normByS_depth"},$Regions[$r]{"normByR_depth_fem"});
									} else {
										$Regions[$r]{"Appel"} = "No_Data";
										print SORTIE "NA\t";
									}
								} else {
									# les femmes ne sont pas considérées pour l'appel de CNV des régions du chrY
									print SORTIE "NA\t";
								}
							} else {
								if (@{ $Regions[$r]{"all_normByS_depths_males"} }) {
									($ratio2center,$ratio2spread) = depth_ratio($center,$spread,$Regions[$r]{$patient}{"normByS_depth"},$Regions[$r]{"normByR_depth_males"});
								} else {
									$Regions[$r]{"Appel"} = "No_Data";
									print SORTIE "NA\t";
								}
							}

						} else {
							# Pour les autosomes
							if ($Regions[$r]{"Chrom"} !~ m/^chr[XY]$|^[XY]$/) {
								if (@Autosomes) {
									($ratio2center,$ratio2spread) = depth_ratio($center,$spread,$Regions[$r]{$patient}{"normByS_depth"},$Regions[$r]{"normByR_depth"});
								} else {
									$Regions[$r]{"Appel"} = "No_Data";
									print SORTIE "NA\t";
								}
							} else {
								if($Patients{$patient}{"Sexe"} eq "F") {
									if ($Regions[$r]{"Chrom"} =~ m/^(chrX|X)$/) {
										if (@{ $Regions[$r]{"all_normByS_depths_fem"} }) {	
											($ratio2center,$ratio2spread) = depth_ratio($center,$spread,$Regions[$r]{$patient}{"normByS_depth"},$Regions[$r]{"normByR_depth_fem"});
										} else {
											$Regions[$r]{"Appel"} = "No_Data";
											print SORTIE "NA\t";
										}
									} else {
										# les femmes ne sont pas considérées pour l'appel de CNV des régions du chrY
										print SORTIE "NA\t";
									}
								} else {
									if (@{ $Regions[$r]{"all_normByS_depths_males"} }) {
										($ratio2center,$ratio2spread) = depth_ratio($center,$spread,$Regions[$r]{$patient}{"normByS_depth"},$Regions[$r]{"normByR_depth_males"});
									} else {
										$Regions[$r]{"Appel"} = "No_Data";
										print SORTIE "NA\t";
									}
								}

							}
						}
						if (defined $ratio2center) { 
							$Regions[$r]{$patient}{"ratio2center"} = $ratio2center;
							if (defined	$ratio2spread) {
								$Regions[$r]{$patient}{"ratio2spread"} = $ratio2spread;
								print SORTIE sprintf("%.3f",$ratio2center)."(".sprintf("%.3f",$ratio2spread).")\t";
							} else {
								print SORTIE sprintf("%.3f",$ratio2center)."\t";
							}
						}

						## test CNV
						my($del,$dup) = CNV_test($center_test,$spread_test,$ratio2center,$ratio2spread,$seuil_del,$seuil_dup,$spread_del,$spread_dup);
						if ($del) {
							$Regions[$r]{"nb_CNV"}{"DEL"}++;
							$nb_evts++;
							$CNV{$patient}++;
							$Results{$patient}{$r} = "DEL";
						} elsif ($dup) {
							$Regions[$r]{"nb_CNV"}{"DUP"}++;
							$nb_evts++;
							$CNV{$patient}++;
							$Results{$patient}{$r} = "DUP";
						}

					}

				}

				## test recurrence CNV dans region
				my$recurrent="";
				if (@{ $Regions[$r]{"all_normByS_depths_fem"} } && @{ $Regions[$r]{"all_normByS_depths_males"} }) {
					if ( ($nb_evts / (scalar@{ $Regions[$r]{"all_normByS_depths_fem"} } + scalar@{ $Regions[$r]{"all_normByS_depths_males"} })) > $seuil_region && $nb_parcours > 0 ) { $recurrent = 1; }
				} elsif (@{ $Regions[$r]{"all_normByS_depths_fem"} } && !@{ $Regions[$r]{"all_normByS_depths_males"} }) {
					if ( ($nb_evts / scalar@{ $Regions[$r]{"all_normByS_depths_fem"} }) > $seuil_region && $nb_parcours > 0 ) { $recurrent = 1; }
				} elsif (!@{ $Regions[$r]{"all_normByS_depths_fem"} } && @{ $Regions[$r]{"all_normByS_depths_males"} }) {
					if ( ($nb_evts / scalar@{ $Regions[$r]{"all_normByS_depths_males"} }) > $seuil_region && $nb_parcours > 0 ) { $recurrent = 1; }
				}
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



			} else {		# ! $ratioByGender

				# Pour chaque patient, un premier parcours permet la pondération intra-patient
				#(en divisant la profondeur de la region par la profondeur totale recalculée du patient)
				@{ $Regions[$r]{"all_normByS_depths"} }=();
				for (my $patient = 0 ; $patient < $patient_nbr ; $patient++) {
					unless ($Patients{$patient}{"ecarte"}) {
						if($Regions[$r]{$patient}{"raw_depth"} > 0) {
							# seules les prof > 0 prises en compte pour la norme
							$Regions[$r]{$patient}{"normByS_depth"} = $Regions[$r]{$patient}{"raw_depth"} / $Patients{$patient}{"Ref_depth"};
							if ( ($Patients{$patient}{"Sexe"} eq "H") && ($Regions[$r]{"Chrom"} =~ m/^(chrX|X)$/) ) {
								push(@{ $Regions[$r]{"all_normByS_depths"} }, ($Regions[$r]{$patient}{"normByS_depth"}*2));
							} elsif ( ($Patients{$patient}{"Sexe"} eq "F") && ($Regions[$r]{"Chrom"} =~ m/^(chrY|Y)$/) ) {
								 next ;
							} else { push(@{ $Regions[$r]{"all_normByS_depths"} }, $Regions[$r]{$patient}{"normByS_depth"}); }
						} else { $Regions[$r]{$patient}{"normByS_depth"} = 0; }
					}
				}

				# Nous calculons pour chaque region la mediane/moyenne de profondeur (ponderee pour chaque patient par rapport aux autres regions)
				if(@{ $Regions[$r]{"all_normByS_depths"} }) {
					$Regions[$r]{"normByR_depth"} = getRegionStats($sortDepths,$CNV_opt_r,\%cnvStats,\@{ $Regions[$r]{"all_normByS_depths"} });
				}

				# Un second parcours permet l'appel de CNV (en comparant la profondeur pour un patient par la profondeur moyenne de la region)
				for (my $patient = 0 ; $patient < $patient_nbr ; $patient++) {

					if ($Patients{$patient}{"ecarte"}) {
						print SORTIE "NA\t";

					} else {

						if (@{ $Regions[$r]{"all_normByS_depths"} }) {

							# calcul ratio sur valeur moy/med
							my($ratio2center,$ratio2spread);
							if($Patients{$patient}{"Sexe"} eq "F") {
								if ($Regions[$r]{"Chrom"} !~ m/^(chrY|Y)$/) {
									($ratio2center,$ratio2spread) = depth_ratio($center,$spread,$Regions[$r]{$patient}{"normByS_depth"},$Regions[$r]{"normByR_depth"});
								} else {
									#les femmes ne sont pas considérées pour l'appel de CNV des régions du chrY
									print SORTIE "NA\t";
								}
							} else {
								#*2 for chrX
								if ($Regions[$r]{"Chrom"} =~ m/^(chrX|X)$/) {
									($ratio2center,$ratio2spread) = depth_ratio($center,$spread,(2*$Regions[$r]{$patient}{"normByS_depth"}),$Regions[$r]{"normByR_depth"});
								} else {
									($ratio2center,$ratio2spread) = depth_ratio($center,$spread,$Regions[$r]{$patient}{"normByS_depth"},$Regions[$r]{"normByR_depth"});
								}
							}
							if (defined $ratio2center) { 
								$Regions[$r]{$patient}{"ratio2center"} = $ratio2center;
								if (defined	$ratio2spread) {
									$Regions[$r]{$patient}{"ratio2spread"} = $ratio2spread;
									print SORTIE sprintf("%.3f",$ratio2center)."(".sprintf("%.3f",$ratio2spread).")\t";
								} else {
									print SORTIE sprintf("%.3f",$ratio2center)."\t";
								}
							}

							## test CNV
							my($del,$dup) = CNV_test($center_test,$spread_test,$ratio2center,$ratio2spread,$seuil_del,$seuil_dup,$spread_del,$spread_dup);
							if ($del) {
								$Regions[$r]{"nb_CNV"}{"DEL"}++;
								$nb_evts++;
								$CNV{$patient}++;
								$Results{$patient}{$r} = "DEL";
							} elsif ($dup) {
								$Regions[$r]{"nb_CNV"}{"DUP"}++;
								$nb_evts++;
								$CNV{$patient}++;
								$Results{$patient}{$r} = "DUP";
							}

						} else {
							$Regions[$r]{"Appel"} = "No_Data";
							print SORTIE "NA\t";
						}

					}

				}

				## test recurrence CNV dans region
				if ( @{ $Regions[$r]{"all_normByS_depths"} } && ($nb_evts/scalar@{ $Regions[$r]{"all_normByS_depths"} }) > $seuil_region && $nb_parcours > 0 ) {
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
print "print results:\n";
foreach my$patient (keys%Results) {
	print "\t".$Patients{$patient}{"ID"}."\n";
	open (CNV1,">$outdir/CNV_$Patients{$patient}{ID}.summary.txt") or die "Pb lors de l'ecriture du fichier sortie $!\n";
	print CNV1 "#Region\tCNV\tN_clean_intervals\tclean_ratio_to_$center\tN_dirty_intervals\tdirty_ratio_to_$center\toverlapping_Samples\n";
	open (CNV2,">$outdir/CNV_$Patients{$patient}{ID}.allIntervals.txt") or die "Pb lors de l'ecriture du fichier sortie $!\n";
	print CNV2 "#Chrom\tStart\tEnd\tLength\tInfo\tinterval_Order\tCNV_type\tratio_to_$center\toccurences";
	if (@cnvFields) {
		foreach (@cnvFields) { 
			#if ($_ eq "norm") {
			#	if ($norm =~ /^std$/i) { print CNV2 "\tavg\tstd"; }
			#	else  { print CNV2 "\t$norm"; }
			#	}
			#else { print CNV2 "\t$_"; }
			print CNV2 "\t$_";
			}
		}
	print CNV2 "\n";
	foreach my$Chrom (@ChrOrder) {
		if (exists $orderedCNV{$patient}{$Chrom}) {
			my$r=0;		##index in @{ $orderedCNV{$patient}{$Chrom} }
			while ($r < scalar@{ $orderedCNV{$patient}{$Chrom} } ) {
				##merge consecutive CNVs
				my($i,$cnvOK,$nonCNVtot,$nextReg_r) = mergeConsecutiveCNV("",$maxNonCNV,$orderedCNV{$patient}{$Chrom}[$r],$regionIndice{$orderedCNV{$patient}{$Chrom}[$r]},\@{ $regionOrder{$Chrom} },\%{ $Results{$patient} },\@Regions);

				if ($maxNonCNVrate) {
					while (($nonCNVtot/($cnvOK+1)) > $maxNonCNVrate) {
						($i,$cnvOK,$nonCNVtot,$nextReg_r) = mergeConsecutiveCNV(($cnvOK-1),$maxNonCNV,$orderedCNV{$patient}{$Chrom}[$r],$regionIndice{$orderedCNV{$patient}{$Chrom}[$r]},\@{ $regionOrder{$Chrom} },\%{ $Results{$patient} },\@Regions);
						}
					}
				##overlapping samples
				my$overlapCNT=0; my%overlapSMPL=();
				for (my$j=0;$j<=$cnvOK;$j++) {
					foreach my$other (keys%Patients) {
						if ( (exists $Results{$other}{${$nextReg_r}[$j]}) && ($Results{$other}{${$nextReg_r}[$j]} eq $Results{$patient}{${$nextReg_r}[0]}) && (!exists $overlapSMPL{$other}) ) { 
							$overlapCNT++;
							$overlapSMPL{$other}=1;
							}
						}
					}
				##print patient's summary and allIntervals
				if ($cnvOK >= ($minCNV-1)) {
					if ($cnvOK == 0) {
						$Result2{$patient}{${$nextReg_r}[0]} = $Results{$patient}{${$nextReg_r}[0]};
						$Result3{$patient}{${$nextReg_r}[0]}{"type"} = $Results{$patient}{${$nextReg_r}[0]};
						$Result3{$patient}{${$nextReg_r}[0]}{"end"} = ${$nextReg_r}[0];
						if (exists $regionIndice{${$nextReg_r}[0]}) {
							print CNV1 $Chrom.":".$Regions[${$nextReg_r}[0]]{"start"}."-".$Regions[${$nextReg_r}[0]]{"end"}."\t".$Results{$patient}{${$nextReg_r}[0]}."\t1\t".sprintf("%.3f",$Regions[${$nextReg_r}[0]]{$patient}{"ratio2center"})."\t.\t.\t$overlapCNT\n";
							print CNV2 $Chrom."\t".$Regions[${$nextReg_r}[0]]{"start"}."\t".$Regions[${$nextReg_r}[0]]{"end"}."\t".($Regions[${$nextReg_r}[0]]{"end"}-$Regions[${$nextReg_r}[0]]{"start"})." bp\t".$Regions[${$nextReg_r}[0]]{"label"}."\t".($regionIndice{${$nextReg_r}[0]}+1)."\t".$Results{$patient}{${$nextReg_r}[0]}."\t".sprintf("%.3f",$Regions[${$nextReg_r}[0]]{$patient}{"ratio2center"})."\t".$Regions[${$nextReg_r}[0]]{"nb_CNV"}{$Results{$patient}{${$nextReg_r}[0]}};
							if (@cnvFields) {
								my$txt = printCNVfields($ratioByGender,\@cnvFields,\%{ $Regions[${$nextReg_r}[0]] },\%{ $Patients{$patient} });
								print CNV2 "$txt";
								}
							print CNV2 "\n\n";
							}
						}
					else {
						my@cleanCNV = (); my@dirtyCNV = ();
						for (my$j=0;$j<=$cnvOK;$j++) {
							if (exists $Results{$patient}{${$nextReg_r}[$j]}) {
								$Result2{$patient}{${$nextReg_r}[$j]} = $Results{$patient}{${$nextReg_r}[$j]};
								push(@cleanCNV, $Regions[${$nextReg_r}[$j]]{$patient}{"ratio2center"});
								push(@dirtyCNV, $Regions[${$nextReg_r}[$j]]{$patient}{"ratio2center"});
								if (exists $regionIndice{${$nextReg_r}[$j]}) {
									print CNV2 $Chrom."\t".$Regions[${$nextReg_r}[$j]]{"start"}."\t".$Regions[${$nextReg_r}[$j]]{"end"}."\t".($Regions[${$nextReg_r}[$j]]{"end"}-$Regions[${$nextReg_r}[$j]]{"start"})." bp\t".$Regions[${$nextReg_r}[$j]]{"label"}."\t".($regionIndice{${$nextReg_r}[$j]}+1)."\t".$Results{$patient}{${$nextReg_r}[$j]}."\t".sprintf("%.3f",$Regions[${$nextReg_r}[$j]]{$patient}{"ratio2center"})."\t".$Regions[${$nextReg_r}[$j]]{"nb_CNV"}{$Results{$patient}{${$nextReg_r}[$j]}};
									if (@cnvFields) {
										my$txt = printCNVfields($ratioByGender,\@cnvFields,\%{ $Regions[${$nextReg_r}[$j]] },\%{ $Patients{$patient} });
										print CNV2 "$txt";
										}
									print CNV2 "\n";
									}
								}
							else	{ 
								$Result2{$patient}{${$nextReg_r}[$j]} = "NA";
								if (exists $regionIndice{${$nextReg_r}[$j]}) {
									print CNV2 $Chrom."\t".$Regions[${$nextReg_r}[$j]]{"start"}."\t".$Regions[${$nextReg_r}[$j]]{"end"}."\t".($Regions[${$nextReg_r}[$j]]{"end"}-$Regions[${$nextReg_r}[$j]]{"start"})." bp\t".$Regions[${$nextReg_r}[$j]]{"label"}."\t".($regionIndice{${$nextReg_r}[$j]}+1)."\t";
									if (exists $Regions[${$nextReg_r}[$j]]{"Appel"}) { print CNV2 "NA\t"; }
									else { print CNV2 "no\t"; }
									if (exists $Regions[${$nextReg_r}[$j]]{$patient}{"ratio2center"}) {
										print CNV2 sprintf("%.3f",$Regions[${$nextReg_r}[$j]]{$patient}{"ratio2center"})."\t".$Regions[${$nextReg_r}[$j]]{"nb_CNV"}{$Results{$patient}{${$nextReg_r}[0]}};
										push(@dirtyCNV, $Regions[${$nextReg_r}[$j]]{$patient}{"ratio2center"});
										}
									else { print CNV2 "na\tna"; }
									if (@cnvFields) {
										my$txt = printCNVfields($ratioByGender,\@cnvFields,\%{ $Regions[${$nextReg_r}[$j]] },\%{ $Patients{$patient} });
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
						print CNV1 $Chrom.":".$Regions[${$nextReg_r}[0]]{"start"}."-".$Regions[${$nextReg_r}[$cnvOK]]{"end"}."\t".$Results{$patient}{${$nextReg_r}[0]}."\t".scalar@cleanCNV."\t".sprintf("%.3f",$cleanAverage)."\t".($cnvOK+1)."\t".sprintf("%.3f",$dirtyAverage)."\t$overlapCNT\n";
						$Result3{$patient}{${$nextReg_r}[0]}{"type"} = $Results{$patient}{${$nextReg_r}[0]};
						$Result3{$patient}{${$nextReg_r}[0]}{"end"} = ${$nextReg_r}[$cnvOK];
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


##summary
open (OUT,">$outdir/CNV.summary.txt") or die "Pb lors de l'ecriture du fichier sortie $!\n";
print OUT "CNV analysis\n";
print OUT "\nparameters:\n";
if ($center_test) {
	print OUT "\tlevel = ratio to $center of all normalized sample depths
	intervals detected as deletion if normalized depth below $seuil_del of $center
	intervals detected as duplication if normalized depth above $seuil_dup of $center\n";
	}
if ($spread_test) {
	print OUT "\tdispersion = $spread of normalized sample depths
	intervals detected as deletion if normalized depth below $center + ($spread_del * $spread)
	intervals detected as duplication if normalized depth above $center + ($spread_dup * $spread)\n";
	}
if ($seuil_patient <1) { print OUT "\tsamples kept if NO more than ".(100*$seuil_patient)."% of intervals are CNVs\n"; }
if ($seuil_region <1) { print OUT "\tintervals kept if NO more than ".(100*$seuil_region)."% of samples are CNVs\n"; }
if ($seuil_cov) { print OUT "\tintervals discarded if all samples are covered less or equal to $seuil_cov\n"; }
if (defined $minDP) { print OUT "\tintervals discarded if avg/med depth of all samples are less or equal to $minDP\n"; }
if ($minCNV) { print OUT "\tCNV calling requires at least $minCNV consecutive called intervals\n"; }
if ($maxNonCNV) { print OUT "\tCNV calling tolerates NO more than $maxNonCNV non-CNV inside\n"; }
if ($ratioByGender eq "all") { print OUT "\tinterval depth normalization by gender for all chromosomes\n"; }
elsif ($ratioByGender eq "gono") { print OUT "\tinterval depth normalization by gender for sex chromosomes only\n"; }
if ($RefDepth eq "mean") { print OUT "\teach sample normalized by the sum of all region mean depths\n"; }
else { print OUT "\teach sample normalized by the sum of all region total depths\n"; }
if ($RefNoChrY) { print OUT "\t(during sample normalization, depths within chrX are doubled for males, and those within chrY are skipped)\n"; }
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
if (defined $minDP || $seuil_cov) { print OUT "\tlow coverage :  $N_lowCov\n"; }

print OUT "\nsamples discarded:\n";
my$N_discard=0;
foreach my$patient (sort(keys%Patients)) {
		if ($Patients{$patient}{"ecarte"}) { print OUT "\t$Patients{$patient}{ID}\n"; $N_discard++; }
		}
unless ($N_discard) { print OUT "\tnone\n"; }

print OUT "\nResults :\n
patient\tCNV\tnber\n";
foreach my$patient (sort(keys%Patients)) {
	unless($Patients{$patient}{"ecarte"}) {
		my$N_dup=0; my$N_del=0;
		foreach (keys%{ $Result3{$patient} }) {
			if ($Result3{$patient}{$_}{"type"} eq "DUP") { $N_dup++; }
			elsif ($Result3{$patient}{$_}{"type"} eq "DEL") { $N_del++; }
			}
		print OUT "\n$Patients{$patient}{ID}: 
		DUP: $N_dup
		DEL: $N_del\n";
		}
	}
close OUT;


##print graph foreach Chrom/CNV
if ($graphByChr) {
	print "drawing graphs:\n";
	##graph by Chrom/CNV ?
	my$Nbr_Reg_max = 0;
	foreach my$Chrom (@ChrOrder) {
		if (scalar@{ $regionOrder{$Chrom} } > $Nbr_Reg_max) {
			$Nbr_Reg_max = scalar@{ $regionOrder{$Chrom} };
			}
		}
	my$graphByCNV = 0;
	if ($Nbr_Reg_max > ${$CNV_opt_r}{"switch2graphByCNV"}) { $graphByCNV = 1; print "\tgraphs by CNV:\n"; }
	else { print "\tgraphs by chrom:\n"; }
	foreach my$patient (keys%Result2)  {
		if ($graphByCNV) {
			graphByCNV2($outdir,$CNV_opt_r,$patient,\%Patients,\@Regions,\@ChrOrder,\%regionOrder,\%Result2,\%Result3);
			}
		else {
			graphByChr2($outdir,$CNV_opt_r,$patient,\%Patients,\@Regions,\@ChrOrder,\%regionOrder,\%Result2);
			}
		}
	}


}



####################
# $Regions[$r]{"normByR_depth"} = getRegionStats($sortDepths,$CNV_opt_r,\%cnvStats,\@Autosomes);
# $Regions[$r]{"normByR_depth"}->{"avg"}
# $$Regions[$r]{"normByR_depth"}{"avg"}
# ${$Regions}[$r]{"normByR_depth"}{"avg"}

sub getRegionStats {

my($sortValues,$CNV_opt_r,$cnvStats_r,$allDepths_r) = @_;

my%stats;

my@selectDepths = ();
if ($sortValues) { @selectDepths = sort{$a<=>$b}@{$allDepths_r}; }
else { @selectDepths = @{$allDepths_r}; }

if (exists ${$cnvStats_r}{"min"}) { $stats{"min"} = $selectDepths[0]; }
if (exists ${$cnvStats_r}{"max"}) {$stats{"max"} = $selectDepths[-1]; }

if (exists ${$cnvStats_r}{"med"}) { $stats{"med"} = mediane(\@selectDepths); }

if (exists ${$cnvStats_r}{"Qtile"} || defined ${$CNV_opt_r}{"range"}) {
	if (scalar@selectDepths > 2) {
		my (@q1Set,@q3Set);
		foreach (@selectDepths) {
			if ($_ < $stats{"med"}) { push(@q1Set,$_); }
			elsif ($_ > $stats{"med"}) { push(@q3Set,$_); }
			}
		$stats{"Q1"} = mediane(\@q1Set);
		$stats{"Q3"} = mediane(\@q3Set);
		#$IQ = $Q3 - $Q1;
		if (defined ${$CNV_opt_r}{"range"}) {
			my@currDepths = ();
			foreach (@selectDepths) {
				## > Q1-r*(med-Q1) && < Q3+r*(Q3-med)
				#if ($_ > ($stats{"Q1"}*(1+$range)-$stats{"med"}*$range) && $_ < ($stats{"Q3"}*(1+$range)-$stats{"med"}*$range)) { push(@selectDepths,$_); }
				## > med-r*(med-Q1) && < med+r*(Q3-med)
				if ($_ > ($stats{"med"}-${$CNV_opt_r}{"range"}*($stats{"med"}-$stats{"Q1"})) && $_ < ($stats{"med"}+${$CNV_opt_r}{"range"}*($stats{"Q3"}-$stats{"med"}))) {
					push(@currDepths, $_);
					}
				}
			@selectDepths = @currDepths;
			}
		}
	else { $stats{"Q1"} = $stats{"med"}; $stats{"Q3"} = $stats{"med"}; }
	if (exists ${$CNV_opt_r}{"spread"} && ${$CNV_opt_r}{"spread"} eq "Qtile") {
		#$stats{"spread_inf"} = $stats{"Q1"}*(1+${$CNV_opt_r}{"spread_del"})-$stats{"med"}*${$CNV_opt_r}{"spread_del"};		#Q1-r*(med-Q1)
		$stats{"spread_inf"} = $stats{"med"} + ${$CNV_opt_r}{"spread_del"}*($stats{"med"}-$stats{"Q1"});					#med-r*(med-Q1)
		#$stats{"spread_sup"} = $stats{"Q3"}*(1+${$CNV_opt_r}{"spread_dup"})-$stats{"med"}*${$CNV_opt_r}{"spread_dup"};		#Q3+r*(Q3-med)
		$stats{"spread_sup"} = $stats{"med"} + ${$CNV_opt_r}{"spread_dup"}*($stats{"Q3"}-$stats{"med"});					#med-r*(Q3-med)
		}
	}

if (exists ${$cnvStats_r}{"avg"}) {
	foreach (@selectDepths) { $stats{"avg"} += $_; } 
	$stats{"avg"} /= scalar@selectDepths;
	}
if (exists ${$cnvStats_r}{"std"}) {
	if (scalar@selectDepths > 1) {
		my$sqtotal = 0;
		foreach (@selectDepths) { $sqtotal += ($stats{"avg"}-$_)**2; }
		$stats{"std"} = ($sqtotal / (scalar@selectDepths-1))**0.5;
		}
	else { $stats{"std"} = 0; }
	if (exists ${$CNV_opt_r}{"spread"} && ${$CNV_opt_r}{"spread"} eq "std") {
		$stats{"spread_inf"} = $stats{"avg"} + ${$CNV_opt_r}{"spread_del"}*$stats{"std"};
		$stats{"spread_sup"} = $stats{"avg"} + ${$CNV_opt_r}{"spread_dup"}*$stats{"std"};
		}
	}

return(\%stats);

}

sub mediane {
my($data_r) = @_;
my$med;
#odd?
if(scalar@{$data_r}%2) 
	{ $med = ${$data_r}[int(scalar@{$data_r}/2)]; }
#even
else { $med = ( ${$data_r}[int(scalar@{$data_r}/2)-1] + ${$data_r}[int(scalar@{$data_r}/2)] )/2; }
return($med);
}

sub average {
my($data_r) = @_;
my$avrg = 0;
foreach (@{$data_r}) { $avrg += $_; } 
$avrg /= scalar@{$data_r};
return($avrg);
}

####################

sub depth_ratio {

my($center,$spread,$sample_depth,$depth_stats_r) = @_;		# ${$depth_stats_r}{"avg"} or $depth_stats_r->{"avg"}

my($ratio2center,$ratio2spread);

if ($spread) {
	if ($spread eq "std") {
		$ratio2center = $sample_depth / $depth_stats_r->{"avg"};
		if ($depth_stats_r->{"std"}) { $ratio2spread = ($sample_depth - $depth_stats_r->{"avg"}) / $depth_stats_r->{"std"}; }
		}
	elsif ($spread eq "Qtile") {
		$ratio2center = $sample_depth / $depth_stats_r->{"med"};
		if ($ratio2center > 1) {
			if ($depth_stats_r->{"Q3"} != $depth_stats_r->{"med"})
				{ $ratio2spread = ($sample_depth - $depth_stats_r->{"med"}) / ($depth_stats_r->{"Q3"} - $depth_stats_r->{"med"}); }
			}
		elsif ($ratio2center < 1) {
			if ($depth_stats_r->{"Q1"} != $depth_stats_r->{"med"})	
				{ $ratio2spread = ($sample_depth - $depth_stats_r->{"med"}) / ($depth_stats_r->{"med"} - $depth_stats_r->{"Q1"}); }
			}
		else { $ratio2spread = 0; }
		}
	}
else {
	$ratio2center = $sample_depth / $depth_stats_r->{$center};
	}

return($ratio2center,$ratio2spread);

}

####################

sub CNV_test {
my($center_test, $spread_test, $ratio2center, $ratio2spread, $seuil_del, $seuil_dup, $spread_del, $spread_dup) = @_;
my($del,$dup);
if ($center_test && $spread_test) {
	if (defined $ratio2center && defined $ratio2spread) {
		if ( $ratio2center == 0 || ($ratio2center < $seuil_del && $ratio2spread < $spread_del) ) { $del++; }
		elsif ($ratio2center > $seuil_dup && $ratio2spread > $spread_dup) { $dup++; }
		}
	}
elsif (!$spread_test) {
	if (defined $ratio2center) {
		if ($ratio2center < $seuil_del)  { $del++; }
		elsif ($ratio2center > $seuil_dup) { $dup++; }
		}
	}
elsif (!$center_test) {
	if (defined $ratio2spread) {
		if ($ratio2center == 0 || $ratio2spread < $spread_del)  { $del++; }
		elsif ($ratio2spread > $spread_dup) { $dup++; }
		}
	}
return($del,$dup);
}

####################
##my$txt = printCNVfields($ratioByGender,\@cnvFields,\%{ $Regions[$nextReg[$j]] },\%{ $Patients{$patient} });

sub printCNVfields {
my($ratioByGender,$cnvFields_r,$Regions_r,$Patients_r) = @_;
my$txt = "";
foreach my$val (@{$cnvFields_r}) {
	if ( !$ratioByGender || ($ratioByGender && ($ratioByGender eq "gono" && ${$Regions_r}{"Chrom"} !~ m/^chr[XY]$|^[XY]$/)) ) {
		if (defined ${$Regions_r}{"normByR_depth"}->{$val}) { $txt .= "\t".sprintf("%.1f",${$Regions_r}{"normByR_depth"}->{$val}); }
		else { $txt .= "\tna"; }
		}
	else {
		if (${$Patients_r}{"Sexe"} eq "F") {
			if (defined ${$Regions_r}{"normByR_depth_fem"}->{$val}) { $txt .= "\t".sprintf("%.1f",${$Regions_r}{"normByR_depth_fem"}->{$val}); }
			else { $txt .= "\tna"; }
			}
		else {
			if (defined ${$Regions_r}{"normByR_depth_males"}->{$val}) { $txt .= "\t".sprintf("%.1f",${$Regions_r}{"normByR_depth_males"}->{$val}); }
			else { $txt .= "\tna"; }
			}
		}

	}
return ($txt);
}


####################

sub mergeConsecutiveCNV {
my($maxI,$maxNonCNV,$orderedCNV,$regionIndice,$regionOrderRef,$ResultsRef,$Regions_r)=@_;
#my@regionOrder = @$h1;
#my%Results = %$h2;
#my@Regions = @$h3;
my$ok=1; my$i=0; my$cnvOK=0; my$nonCNV=0; my$nonCNVtot=0;
my@nextReg=($orderedCNV);
while ($ok) {
	$i++;
	if (exists ${$regionOrderRef}[($regionIndice + $i)]) {
		push(@nextReg, ${$regionOrderRef}[($regionIndice + $i)]);
		if (exists ${$ResultsRef}{$nextReg[$i]}) {
			## same CNV type ?
			if(${$ResultsRef}{$nextReg[$i]} eq ${$ResultsRef}{$nextReg[0]}) {
				$cnvOK = $i;
				#$nonCNVtot += $nonCNV;
				$nonCNV=0;		##to reset
				}
			else { $ok=0; }
			}
		else {
			if (!exists ${$Regions_r}[$nextReg[$i]]{"Appel"}) { 
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

sub graphByChr1 {

my($nGraf,$outdir,$CNV_opt_r,$file,$Files_r,$sampleName_r,$Regions_r,$ChromOrder_r,$regionOrderRef,$PatientsRef,$ResultsRef)= @_;

print "\tcmdR, for sample ${$sampleName_r}{$file}\n";

my$center = ${$CNV_opt_r}{"center"};
my$spread = ${$CNV_opt_r}{"spread"};
my$seuil_del = ${$CNV_opt_r}{"seuil_del"};
my$seuil_dup = ${$CNV_opt_r}{"seuil_dup"};
my$spread_del = ${$CNV_opt_r}{"spread_del"};
my$spread_dup = ${$CNV_opt_r}{"spread_dup"};
my$maxDepthGraph = ${$CNV_opt_r}{"maxDepthGraph"};

##all in 1 sheet:
#my$Nbr_Chr= scalar(keys%{$regionOrderRef});
#if ($nGraf eq "max") { $nGraf = $Nbr_Chr; }

##all in 1 sheet:
#my$maxX=0;
#foreach my$Chrom (keys%{$regionOrderRef}) {	
#	if (scalar@{ ${$regionOrderRef}{$Chrom} } > $maxX) { $maxX = scalar@{ ${$regionOrderRef}{$Chrom} }; }
#	}

##all in 1 sheet:
#my$cmdR = "";
#my$c=1; #chr iteration
#my$n=1; #chr iteration, stepped back to 0 each time a graph is done
#my$N=1; #graph iteration

##1 sheet / chr :
open (CMDR, ">$outdir/${$sampleName_r}{$file}\_temp.R") || die;
print CMDR "#!/usr/bin/env Rscript\n\n" ;
print CMDR "pdf(\"$outdir/CNV_".${$sampleName_r}{$file}.".pdf\", width=11.69, height=4,135)\n";

foreach my$Chrom (@{$ChromOrder_r}) {

	$Chrom =~ s/^chr//i;
	
	if (exists ${$regionOrderRef}{$Chrom}) {

		my$cmdR = "";

		my$Nbr_Reg = scalar@{ ${$regionOrderRef}{$Chrom} };
		##1 sheet / chr :
		my$maxX = $Nbr_Reg;

		my$maxYsup=$seuil_dup; #my$maxYinf=$seuil_del;	
		foreach my$region (@{ ${$regionOrderRef}{$Chrom} }) {
			foreach my$f (@{$Files_r}) {
				if (exists ${$Regions_r}[$region]{$f}{"ratio2center"}) {
					if (${$Regions_r}[$region]{$f}{"ratio2center"} > $maxYsup)
						{ $maxYsup = ${$Regions_r}[$region]{$f}{"ratio2center"}; }
					#if ($normGraf eq "std") {
					#	if (${$Regions_r}[$region]{$f}{"ratio2center"} < $maxYinf)
					#		{ $maxYinf = ${$Regions_r}[$region]{$f}{"ratio2center"}; }
					#	}
					}
				}
			}
		if ($maxDepthGraph && $maxYsup > $maxDepthGraph) { $maxYsup = $maxDepthGraph; }
		#if ($normGraf eq "std") {
		#	if ($maxDepthGraph && $maxYinf < (-$maxDepthGraph)) { $maxYinf = (-$maxDepthGraph); }
		#	}

		#if ($normGraf eq "std") {
			##all in 1 sheet:
			#$cmdR .= "par(fig=c(0,1,".(1-(($n-0.05)/$nGraf)).",".(1-(($n-0.95)/$nGraf))."), new=TRUE)
#plot (c(0,0), xlim=c(0,$maxX), ylim=c($maxYinf,$maxYsup), type =\"n\", main=\"chrom: $Chrom\", xlab=\"\", ylab=\"depth_ratio_to_$normGraf\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n";
			##1 sheet / chr :
		#	$cmdR .= "plot (c(0,0), xlim=c(0,$maxX), ylim=c($maxYinf,$maxYsup), type =\"n\", main=\"chrom: $Chrom\", xlab=\"\", ylab=\"depth_ratio_to_$normGraf\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n";
		#	}

		##all in 1 sheet:
		#$cmdR .= "par(fig=c(0,1,".(1-(($n-0.05)/$nGraf)).",".(1-(($n-0.95)/$nGraf))."), new=TRUE)
#plot (c(0,0), xlim=c(0,$maxX), ylim=c(0,$maxYsup), type =\"n\", main=\"chrom: $Chrom\", xlab=\"\", ylab=\"depth_ratio_to_$normGraf\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n";
		##1 sheet / chr :
		$cmdR .= "plot (c(0,0), xlim=c(0,$maxX), ylim=c(0,$maxYsup), type =\"n\", main=\"chrom: $Chrom\", xlab=\"\", ylab=\"depth_ratio_to_$center\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n";

		##gene separations
		my$currentGene=""; my$tmpTxt=""; my$Nbr_gene=0;
		for (my$r=0;$r<$Nbr_Reg;$r++) {
			if (${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"}) {
				if (${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"} ne "NA" && ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"} ne $currentGene)  {
					$tmpTxt .= "abline(v=".($r+0.5).", col=\"blue\", lty = \"dotted\", lwd=1)\n";
					$currentGene = ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"};
					$Nbr_gene++;
					}
				}
			}
		if ($Nbr_gene < ${$CNV_opt_r}{"maxGeneSep"}) { $cmdR .= $tmpTxt; }

		my$Nbr_CNV=0;
		for (my$r=0;$r<$Nbr_Reg;$r++) {
			if (exists ${$ResultsRef}{$file}{${$regionOrderRef}{$Chrom}[$r]})
				{ $Nbr_CNV++; }
			}

		##region labels
		my@printReg=();
		if ($maxX < ${$CNV_opt_r}{"maxGeneLab"}) {
			##in grey if invalid
			for (my$r=0;$r<$Nbr_Reg;$r++) { 
				if (defined ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"Appel"}) { push(@printReg,$r) ; }
				}
			if (@printReg) {
				$cmdR .= "axis(1, at=c(";
				foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
				chop $cmdR;
				$cmdR .="), labels=c(";
				foreach my$r (@printReg) { $cmdR .= "\"".${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"label"}."\","; }
				chop $cmdR;
				$cmdR .= "), col.axis=\"darkgrey\", las=2";
				if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
				else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_Reg))." )\n"; }
				}
			##in black if valid
			@printReg=();
			for (my$r=0;$r<$Nbr_Reg;$r++) { 
				if(!defined ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"Appel"}) { push(@printReg,$r) ; }
				}
			if (@printReg) {
				$cmdR .= "axis(1, at=c(";
				foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
				chop $cmdR;
				$cmdR .="), labels=c(";
				foreach my$r (@printReg) { $cmdR .= "\"".${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"label"}."\","; }
				chop $cmdR;
				$cmdR .= "), col.axis=\"black\", las=2";
				if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
				else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_Reg))." )\n"; }
				}
			}
		else {
			##in grey if invalid; only ticks
			for (my$r=0;$r<$Nbr_Reg;$r++) { 
				if (defined ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"Appel"}) { push(@printReg,$r) ; }
				}
			if (@printReg && scalar@printReg < ${$CNV_opt_r}{"maxGeneSep"}) {
				$cmdR .= "axis(1, at=c(";
				foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
				chop $cmdR;
				$cmdR .= "), labels = FALSE, col.ticks=\"darkgrey\")\n";
				}
			##in black if valid
			@printReg=();
			if ($Nbr_gene < ${$CNV_opt_r}{"maxGeneLab"}) { #&& ($Nbr_gene+$Nbr_CNV)>=${$CNV_opt_r}{"maxGeneLab"}) {
				$currentGene="";
				for (my$r=0;$r<$Nbr_Reg;$r++) {
					if (${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"}) {
						if (${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"} ne "NA" && ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"} ne $currentGene)  { 
							push(@printReg,$r); 
							$currentGene = ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"};
							}
						}
					}
				}
			#elsif ($Nbr_gene<${$CNV_opt_r}{"maxGeneLab"} && ($Nbr_gene+$Nbr_CNV)<${$CNV_opt_r}{"maxGeneLab"}) {
			#	for (my$r=0;$r<$Nbr_Reg;$r++) {
			#		if ( (${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"} && ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"} ne "NA" && ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"} ne $currentGene) || (exists ${$ResultsRef}{$file}{${$regionOrderRef}{$Chrom}[$r]}) ) {
			#			push(@printReg,$r);
			#			$currentGene = ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"};
			#			}
			#		}
			#	}
			if (@printReg) {
				$cmdR .= "axis(1, at=c(";
				foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
				chop $cmdR;
				$cmdR .="), labels=c(";
				foreach my$r (@printReg) { $cmdR .= "\"".${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"}."\","; }
				chop $cmdR;
				$cmdR .= "), col.axis=\"black\", las=2";
				if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
				else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_gene))." )\n"; }
				}
			##in red if CNV; only ticks
			@printReg=();
			for (my$r=0;$r<$Nbr_Reg;$r++) {
				if (exists ${$ResultsRef}{$file}{${$regionOrderRef}{$Chrom}[$r]}) { push(@printReg,$r) ; }
				}
			if (@printReg && scalar@printReg < ${$CNV_opt_r}{"maxGeneSep"}) {
				$cmdR .= "axis(1, at=c(";
				foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
				chop $cmdR;
				$cmdR .= "), labels = FALSE, col.ticks=\"red\")\n";
				}
			
			
			}

		##all not target sample lines (grey):
		for (my$f=0;$f<scalar@{$Files_r};$f++) {
			unless (${$Files_r}[$f] eq $file) {
				my$r1=0;
				while ($r1<$Nbr_Reg) {
					my$r2=$r1;
					while ($r2<$Nbr_Reg) {
						if (exists ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r2]]{${$Files_r}[$f]}{"ratio2center"}) { $r2++;}
						else { last; }
						}
					if (($r2-1) > $r1) {
						$cmdR .= "lines( c(";
						for (my$r=$r1;$r<$r2;$r++)
							{ $cmdR .= ($r+1).","; }
						chop $cmdR;
						$cmdR .= "), c(";
						for (my$r=$r1;$r<$r2;$r++)
							{ $cmdR .= ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{${$Files_r}[$f]}{"ratio2center"}.","; }
						chop $cmdR;
						$cmdR .= "), type =\"l\", lwd=1, col=\"darkgrey\")\n";
						}
					elsif (($r2-1) == $r1) {
						$cmdR .= "lines( c(".($r1+1)."), c(".${$Regions_r}[${$regionOrderRef}{$Chrom}[$r1]]{${$Files_r}[$f]}{"ratio2center"}."), type =\"p\", lwd=1, col=\"darkgrey\")\n";
						}
					$r1 = ($r2+1);
					}
				}
			}

		##threshold lines (black):
		if (${$CNV_opt_r}{"spread_test"}) {
			foreach my$gender ("normByR_depth","normByR_depth_fem","normByR_depth_males") {
				my$r1=0;
				while ($r1<$Nbr_Reg) {
					my$r2=$r1;
					while ($r2<$Nbr_Reg) {
						if (${$Regions_r}[${$regionOrderRef}{$Chrom}[$r2]]{$gender}->{$center}) { $r2++;}
						else { last; }
						}
					if (($r2-1) > $r1) {
						foreach my$fold ($spread_dup,$spread_del) {
							$cmdR .= "lines( c(";
							for (my$r=$r1;$r<$r2;$r++)
								{ $cmdR .= ($r+1).","; }
							chop $cmdR;
							$cmdR .= "), c(";
							for (my$r=$r1;$r<$r2;$r++) {
								if ( ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{$gender}->{$center} )
									{ $cmdR .= (1+$fold*(${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{$gender}->{$spread} / ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{$gender}->{$center})).","; }
								else 	{ $cmdR .= "1,"; }
								}
							chop $cmdR;
							$cmdR .= "), type =\"l\", lwd=1, col=\"black\")\n";
							}
						}
					elsif (($r2-1) == $r1) {
						$cmdR .= "lines( c(".($r1+1)."), c(".(1+(${$Regions_r}[${$regionOrderRef}{$Chrom}[$r1]]{$gender}->{$spread} / ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r1]]{$gender}->{$center}))."), type =\"p\", lwd=1, col=\"black\")\n";
						$cmdR .= "lines( c(".($r1+1)."), c(".(1-(${$Regions_r}[${$regionOrderRef}{$Chrom}[$r1]]{$gender}->{$spread} / ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r1]]{$gender}->{$center}))."), type =\"p\", lwd=1, col=\"black\")\n";
						}
					$r1 = ($r2+1);
					}
				}
			}

	if (${$CNV_opt_r}{"center_test"}) {
			$cmdR .= "abline(h=$seuil_del, col=\"black\", lty = \"dashed\", lwd=1)\n";
			$cmdR .= "abline(h=$seuil_dup, col=\"black\", lty = \"dashed\", lwd=1)\n";
			}

		##target sample line (green)
		my$r1=0;
		while ($r1<$Nbr_Reg) {
			my$r2=$r1;
			while ($r2<$Nbr_Reg) {
				if (exists ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r2]]{$file}{"ratio2center"}) { $r2++;}
				else { last; }
				}
			if (($r2-1) > $r1) {
				$cmdR .= "lines( c(";
				for (my$r=$r1;$r<$r2;$r++)
					{ $cmdR .= ($r+1).","; }
				chop $cmdR;
				$cmdR .= "), c(";
				for (my$r=$r1;$r<$r2;$r++)
					{ $cmdR .= ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{$file}{"ratio2center"}.","; }
				chop $cmdR;
				$cmdR .= "), type =\"l\", lwd=1, col=\"green\")\n";
				}
			elsif (($r2-1) == $r1) {
				$cmdR .= "lines( c(".($r1+1)."), c(".${$Regions_r}[${$regionOrderRef}{$Chrom}[$r1]]{$file}{"ratio2center"}."), type =\"p\", lwd=1, col=\"green\")\n";
				}
			$r1 = ($r2+1);
			}
		##points for CNVs
		my$points .= "points( c(";
		for (my$r=0;$r<$Nbr_Reg;$r++) {
			if (exists ${$ResultsRef}{$file}{${$regionOrderRef}{$Chrom}[$r]} && exists ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{$file}{"ratio2center"})
				{ $points .= ($r+1).","; }
			}
		chop $points;
		$points .= "), c(";
		foreach my$region (@{ ${$regionOrderRef}{$Chrom} }) {
			if (exists ${$ResultsRef}{$file}{$region} && exists ${$Regions_r}[$region]{$file}{"ratio2center"})
				{ $points .= ${$Regions_r}[$region]{$file}{"ratio2center"}.","; }
			}
		chop $points;
		$points .= "), type =\"p\", pch = 16, lwd=2, col=\"red\")\n";
		if ($Nbr_CNV) { $cmdR .= $points; }

		#if ($normGraf eq "std") { $cmdR .= "abline(h=0, col=\"darkgrey\", lty = \"dashed\", lwd=1)\n"; }
		#else { $cmdR .= "abline(h=1, col=\"darkgrey\", lty = \"dashed\", lwd=1)\n"; }
		$cmdR .= "abline(h=1, col=\"darkgrey\", lty = \"dashed\", lwd=1)\n";

		##all in 1 sheet:
		#if ($c==$Nbr_Chr || $n==$nGraf) {
		#	open (CMDR, ">$outdir/${$sampleName_r}{$file}\_temp.R") || die;
		#	print CMDR "#!/usr/bin/env Rscript\n\n" ;
		#	if ($nGraf==$Nbr_Chr) { print CMDR "pdf(\"".$outdir."/CNV_${$sampleName_r}{$file}.pdf\", width=11.69, height=".($nGraf*3).")\npar(mfrow=c($nGraf,1))\n"; }
		#	else {
		#		if ($N>1) { print CMDR "pdf(\"".$outdir."/CNV_${$sampleName_r}{$file}\_$N.pdf\", width=11.69, height=".($nGraf*3).")\npar(mfrow=c($nGraf,1))\n"; }
		#		else { print CMDR "pdf(\"".$outdir."/CNV_${$sampleName_r}{$file}\_$N.pdf\", width=11.69, height=".($n*3).")\npar(mfrow=c($n,1))\n"; }
		#		}
		#	print CMDR "$cmdR";
		#	print CMDR "title(main=\"sample: ${$sampleName_r}{$file}";
		#	if (${$PatientsRef}{$file}{"ecarte"}) { print CMDR " (invalid)\", col.main=\"red\""; }
		#	else { print CMDR "\""; }
		#	print CMDR ", outer=TRUE, line=-2, cex.main=2)\n";
		#	print CMDR "dev.off()\nquit(save=\"no\")\n";
		#	close CMDR;
		#	system "Rscript $outdir/${$sampleName_r}{$file}\_temp.R";
		#	unlink "$outdir/${$sampleName_r}{$file}\_temp.R";
		#	$cmdR="";
		#	$n=0;
		#	$N++;
		#	}
		#$c++;$n++;

		##1 sheet / chr :
		print CMDR "$cmdR";

		}	
	}
##1 sheet / chr :
#print CMDR "title(main=\"sample: ${$sampleName_r}{$file}";
#if (${$PatientsRef}{$file}{"ecarte"}) { print CMDR " (invalid)\", col.main=\"red\""; }
#else { print CMDR "\""; }
#print CMDR ", outer=TRUE, line=-2, cex.main=2)\n";
print CMDR "dev.off()\nquit(save=\"no\")\n";
close CMDR;
system "Rscript $outdir/${$sampleName_r}{$file}\_temp.R";
unlink "$outdir/${$sampleName_r}{$file}\_temp.R";

}





####################
sub graphByChr2 {

my($outdir,$CNV_opt_r,$patient,$PatientsRef,$Regions_r,$ChrOrderRef,$regionOrderRef,$ResultsRef)= @_;

print "\tcmdR, for ${$PatientsRef}{$patient}{ID}\n";

my$center = ${$CNV_opt_r}{"center"};
my$spread = ${$CNV_opt_r}{"spread"};
my$seuil_del = ${$CNV_opt_r}{"seuil_del"};
my$seuil_dup = ${$CNV_opt_r}{"seuil_dup"};
my$spread_del = ${$CNV_opt_r}{"spread_del"};
my$spread_dup = ${$CNV_opt_r}{"spread_dup"};
my$maxDepthGraph = ${$CNV_opt_r}{"maxDepthGraph"};

##all in 1 sheet:
#my$maxX=0;
#foreach my$Chrom (keys%{$regionOrderRef}) {	
#	if (scalar@{ $regionOrderRef->{$Chrom} } > $maxX) { $maxX = scalar@{ $regionOrderRef->{$Chrom} }; }		#or @{ ${$regionOrderRef}{$Chrom} }?
#	}

##all in 1 sheet:
#my$Nbr_Chr= scalar(keys%{$regionOrderRef});

open (CMDR, ">$outdir/${$PatientsRef}{$patient}{ID}\_temp.R") || die;
print CMDR "#!/usr/bin/env Rscript\n\n" ;
##all in 1 sheet:
#print CMDR "pdf(\"".$outdir."/CNV_".${$PatientsRef}{$patient}{"ID"}.".pdf\", width=11.69, height=".($Nbr_Chr*3).")\n
#par(mfrow=c($Nbr_Chr,1))\n";		#A4 size print measures 21.0 x 29.7cm, 8.27 x 11.69 inches

##1 sheet / chr
print CMDR "pdf(\"".$outdir."/CNV_".${$PatientsRef}{$patient}{"ID"}.".pdf\", width=11.69, height=4,135)\n";


for my$Chrom (@{$ChrOrderRef}) {	

	my$cmdR = "";
	my$Nbr_Reg = scalar@{ $regionOrderRef->{$Chrom} };
	##1 sheet / chr :
	my$maxX = $Nbr_Reg;

	my$maxYsup=$seuil_dup;
	foreach my$region (@{ $regionOrderRef->{$Chrom} }) {
		for (my$p=0;$p<scalar(keys%{$PatientsRef});$p++) {
			if (exists ${$Regions_r}[$region]{$p}{"ratio2center"}) {
				if (${$Regions_r}[$region]{$p}{"ratio2center"} > $maxYsup)
					{ $maxYsup = ${$Regions_r}[$region]{$p}{"ratio2center"}; }
				}
			}
		}
	if ($maxDepthGraph && $maxYsup > $maxDepthGraph) { $maxYsup = $maxDepthGraph; }

	my$ChrName = $Chrom; $ChrName =~ s/^chr//;

	##all in 1 sheet:
#		$cmdR .= "par(fig=c(0,1,".(1-(($c+0.95)/$Nbr_Chr)).",".(1-(($c+0.05)/$Nbr_Chr))."), new=TRUE)
#plot (c(0,0), xlim=c(0,$maxX), ylim=c(0,$maxYsup), type =\"n\", main=\"chrom: $ChrName\", xlab=\"\", ylab=\"depth_ratio_to_$normGraf\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n";
	$cmdR .= "plot (c(0,0), xlim=c(0,$maxX), ylim=c(0,$maxYsup), type =\"n\", main=\"chrom: $ChrName\", xlab=\"\", ylab=\"depth_ratio_to_$center\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n";

	#gene vertical separations
	my$currentGene=""; my$tmpTxt=""; my$Nbr_gene=0;
	for (my$r=0;$r<$Nbr_Reg;$r++) {
		if (${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"} ne "NA" && ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"geneID"} ne $currentGene) {
			$tmpTxt .= "abline(v=".($r+0.5).", col=\"blue\", lty = \"dotted\", lwd=1)\n";
			$currentGene = ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"geneID"};
			$Nbr_gene++;
			}
		}
	if ($Nbr_gene < ${$CNV_opt_r}{"maxGeneSep"}) { $cmdR .= $tmpTxt; }

	my$Nbr_CNV=0;
	for (my$r=0;$r<$Nbr_Reg;$r++) {
		if (exists ${$ResultsRef}{$patient}{${$regionOrderRef}{$Chrom}[$r]})
			{ $Nbr_CNV++; }
		}

	#x labels
	my@printReg=();
	if ($maxX < ${$CNV_opt_r}{"maxGeneLab"}) {
		##in grey if invalid
		for (my$r=0;$r<$Nbr_Reg;$r++) { 
			if (defined ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"Appel"}) { push(@printReg,$r) ; }
			}
		if (@printReg) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
			chop $cmdR;
			$cmdR .="), labels=c(";
			foreach my$r (@printReg) { $cmdR .= "\"".${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"label"}."\","; }
			chop $cmdR;
			$cmdR .= "), col.axis=\"darkgrey\", las=2";
			if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
			else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_Reg))." )\n"; }
			}
		##in black if valid
		@printReg=();
		for (my$r=0;$r<$Nbr_Reg;$r++) { 
			if(!defined ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"Appel"}) { push(@printReg,$r) ; }
			}
		if (@printReg) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
			chop $cmdR;
			$cmdR .="), labels=c(";
			foreach my$r (@printReg) { $cmdR .= "\"".${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"label"}."\","; }
			chop $cmdR;
			$cmdR .= "), col=\"black\", las=2";
			if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
			else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_Reg))." )\n"; }
			}
		}
	else {
		##in grey if invalid; only ticks
		for (my$r=0;$r<$Nbr_Reg;$r++) { 
			if (defined ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"Appel"}) { push(@printReg,$r) ; }
			}
		if (@printReg && scalar@printReg < ${$CNV_opt_r}{"maxGeneSep"}) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
			chop $cmdR;
			$cmdR .= "), labels = FALSE, col.ticks=\"darkgrey\")\n";
			}
		##in black if valid
		@printReg=();
		if ($Nbr_gene < ${$CNV_opt_r}{"maxGeneLab"}) {	# && ($Nbr_gene+$Nbr_CNV)>=${$CNV_opt_r}{"maxGeneLab"}) {
			$currentGene="";
			for (my$r=0;$r<$Nbr_Reg;$r++) {
				if (${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"} ne "NA" && ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"geneID"} ne $currentGene) {
					push(@printReg,$r);
					$currentGene = ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"geneID"};
					}
				}
			}
		#elsif ($Nbr_gene<${$CNV_opt_r}{"maxGeneLab"} && ($Nbr_gene+$Nbr_CNV)<${$CNV_opt_r}{"maxGeneLab"}) {
		#	for (my$r=0;$r<$Nbr_Reg;$r++) {
		#		if ( (${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"} ne "NA" && ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"geneID"} ne $currentGene) || (exists ${$ResultsRef}{$patient}{${$regionOrderRef}{$Chrom}[$r]}) ) {
		#			push(@printReg,$r);
		#			$currentGene = ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"geneID"};
		#			}
		#		}
		#	}
		if (@printReg) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
			chop $cmdR;
			$cmdR .="), labels=c(";
			foreach my$r (@printReg) { $cmdR .= "\"".${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{"label"}."\","; }
			chop $cmdR;
			$cmdR .= "), col.axis=\"black\", las=2";
			if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
			else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_gene))." )\n"; }
			}
		##in red if CNV; only ticks
		@printReg=();
		for (my$r=0;$r<$Nbr_Reg;$r++) {
			if (exists ${$ResultsRef}{$patient}{${$regionOrderRef}{$Chrom}[$r]}) { push(@printReg,$r) ; }
			}
		if (@printReg && scalar@printReg < ${$CNV_opt_r}{"maxGeneSep"}) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
			chop $cmdR;
			$cmdR .= "), labels = FALSE, col.ticks=\"red\")\n";
			}
		}

	#all not-target sample lines (grey):
	for (my$p=0;$p<scalar(keys%{$PatientsRef});$p++) {
		unless ($p == $patient) {

			my$r1=0;
			while ($r1<$Nbr_Reg) {
				my$r2=$r1;
				while ($r2<$Nbr_Reg) {
					if (exists ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r2]]{$p}{"ratio2center"}) { $r2++;}
					else { last; }
					}
				if (($r2-1) > $r1) {
					$cmdR .= "lines( c(";
					for (my$r=$r1;$r<$r2;$r++)
						{ $cmdR .= ($r+1).","; }
					chop $cmdR;
					$cmdR .= "), c(";
					for (my$r=$r1;$r<$r2;$r++)
						{ $cmdR .= ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{$p}{"ratio2center"}.","; }
					chop $cmdR;
					$cmdR .= "), type =\"l\", lwd=1, col=\"darkgrey\")\n";
					}
				elsif (($r2-1) == $r1) {
					$cmdR .= "lines( c(".($r1+1)."), c(".${$Regions_r}[${$regionOrderRef}{$Chrom}[$r1]]{$p}{"ratio2center"}."), type =\"p\", lwd=1, col=\"darkgrey\")\n";
					}
				$r1 = ($r2+1);
				}
			}
		}

	##threshold lines (black):
	if (${$CNV_opt_r}{"spread_test"}) {
		foreach my$gender ("normByR_depth","normByR_depth_fem","normByR_depth_males") {
			my$r1=0;
			while ($r1<$Nbr_Reg) {
				my$r2=$r1;
				while ($r2<$Nbr_Reg) {
					if (${$Regions_r}[${$regionOrderRef}{$Chrom}[$r2]]{$gender}->{$center}) { $r2++; }
					else { last; }
					}
				if (($r2-1) > $r1) {
					foreach my$fold ($spread_dup,$spread_del) {
						$cmdR .= "lines( c(";
						for (my$r=$r1;$r<$r2;$r++)
							{ $cmdR .= ($r+1).","; }
						chop $cmdR;
						$cmdR .= "), c(";
						for (my$r=$r1;$r<$r2;$r++) {
							if ( ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{$gender}->{$center} )
								{ $cmdR .= (1+$fold*(${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{$gender}->{$spread} / ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{$gender}->{$center})).","; }
							else 	{ $cmdR .= "1,"; }
							}
						chop $cmdR;
						$cmdR .= "), type =\"l\", lwd=1, col=\"black\")\n";
						}
					}
				elsif (($r2-1) == $r1) {
					$cmdR .= "lines( c(".($r1+1)."), c(".(1+(${$Regions_r}[${$regionOrderRef}{$Chrom}[$r1]]{$gender}->{$spread} / ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r1]]{$gender}->{$center}))."), type =\"p\", lwd=1, col=\"black\")\n";
					$cmdR .= "lines( c(".($r1+1)."), c(".(1-(${$Regions_r}[${$regionOrderRef}{$Chrom}[$r1]]{$gender}->{$spread} / ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r1]]{$gender}->{$center}))."), type =\"p\", lwd=1, col=\"black\")\n";
					}
				$r1 = ($r2+1);
				}
			}
		}

	if (${$CNV_opt_r}{"center_test"}) {
		$cmdR .= "abline(h=$seuil_del, col=\"black\", lty = \"dashed\", lwd=1)\n";
		$cmdR .= "abline(h=$seuil_dup, col=\"black\", lty = \"dashed\", lwd=1)\n";
		}

	#target sample line (green):
	my$r1=0;
	while ($r1<$Nbr_Reg) {
		my$r2=$r1;
		while ($r2<$Nbr_Reg) {
			if (exists ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r2]]{$patient}{"ratio2center"}) { $r2++;}
			else { last; }
			}
		if (($r2-1) > $r1) {
			$cmdR .= "lines( c(";
			for (my$r=$r1;$r<$r2;$r++)
				{ $cmdR .= ($r+1).","; }
			chop $cmdR;
			$cmdR .= "), c(";
			for (my$r=$r1;$r<$r2;$r++)
				{ $cmdR .= ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{$patient}{"ratio2center"}.","; }
			chop $cmdR;
			$cmdR .= "), type =\"l\", lwd=1, col=\"green\")\n";
			}
		elsif (($r2-1) == $r1) {
			$cmdR .= "lines( c(".($r1+1)."), c(".${$Regions_r}[${$regionOrderRef}{$Chrom}[$r1]]{$patient}{"ratio2center"}."), type =\"p\", lwd=1, col=\"green\")\n";
			}
		$r1 = ($r2+1);
		}

	#points for CNVs
	my$points .= "points( c(";
	for (my$r=0;$r<$Nbr_Reg;$r++) {
		if (exists ${$ResultsRef}{$patient}{${$regionOrderRef}{$Chrom}[$r]} && exists ${$Regions_r}[${$regionOrderRef}{$Chrom}[$r]]{$patient}{"ratio2center"})
			{ $points .= ($r+1).","; }
		}
	chop $points;
	$points .= "), c(";
	foreach my$region (@{ ${$regionOrderRef}{$Chrom} }) {
		if (exists ${$ResultsRef}{$patient}{$region} && exists ${$Regions_r}[$region]{$patient}{"ratio2center"})
			{ $points .= ${$Regions_r}[$region]{$patient}{"ratio2center"}.","; }
		}
	chop $points;
	$points .= "), type =\"p\", pch = 16, lwd=2, col=\"red\")\n";
	if ($Nbr_CNV) { $cmdR .= $points; }

	$cmdR .= "abline(h=1, col=\"darkgrey\", lty = \"dashed\", lwd=1)\n";

	print CMDR "$cmdR";

	}

#print CMDR "title(main=\"sample: ${$PatientsRef}{$patient}{ID}";
#if (${$PatientsRef}{$patient}{"ecarte"}) { print CMDR " (invalid)\", col.main=\"red\""; }
#else { print CMDR "\""; }
#print CMDR ", outer=TRUE, line=-2, cex.main=2)\n";
print CMDR "dev.off()\nquit(save=\"no\")\n";
close CMDR;
system "Rscript $outdir/${$PatientsRef}{$patient}{ID}\_temp.R";
unlink "$outdir/${$PatientsRef}{$patient}{ID}\_temp.R";

}



####################

sub graphByCNV1 {

my($outdir,$CNV_opt_r,$file,$Files_r,$sampleName_r,$Regions_r,$ChromOrder_r,$regionOrderRef,$PatientsRef,$Result2_r,$Result3_r)= @_;

print "\tcmdR, for sample ${$sampleName_r}{$file}\n";

my$center = ${$CNV_opt_r}{"center"};
my$spread = ${$CNV_opt_r}{"spread"};
my$seuil_del = ${$CNV_opt_r}{"seuil_del"};
my$seuil_dup = ${$CNV_opt_r}{"seuil_dup"};
my$spread_del = ${$CNV_opt_r}{"spread_del"};
my$spread_dup = ${$CNV_opt_r}{"spread_dup"};
my$ext = ${$CNV_opt_r}{"graphCNVpadding"};
my$maxDepthGraph = ${$CNV_opt_r}{"maxDepthGraph"};
my$ploidy = ${$CNV_opt_r}{"ploidy"};

open (CMDR, ">$outdir/${$sampleName_r}{$file}\_temp.R") || die;
print CMDR "#!/usr/bin/env Rscript\n\n" ;
print CMDR "pdf(\"$outdir/CNV_".${$sampleName_r}{$file}.".pdf\", width=11.69, height=4,135)\n";

foreach my$CNV (sort{$a<=>$b}keys%{ ${$Result3_r}{$file} }) {

	my$cmdR = "";

	my$Chrom = ${$Regions_r}[$CNV]{"Chrom"};
	$Chrom =~ s/^chr//i;

	my$firstR = $CNV;
	for (my$r=($CNV - $ext);$r<$CNV;$r++) {
		if (defined ${$Regions_r}[$r] && ${$Regions_r}[$r]{"Chrom"} eq $Chrom) {
			$firstR = $r; last;
			}
		}
	my$CNVend = ${$Result3_r}{$file}{$CNV}{"end"};
	my$lastR = $CNVend;
	for (my$r=($CNVend + $ext);$r>=$CNVend;$r--) {
		if (defined ${$Regions_r}[$r] && ${$Regions_r}[$r]{"Chrom"} eq $Chrom) {
			$lastR = $r; last;
			}
		}
	my$Nbr_Reg = $lastR - $firstR + 1;

	my$maxYsup=$seuil_dup;
	foreach my$region (@{ ${$regionOrderRef}{$Chrom} }) {
		foreach my$f (@{$Files_r}) {
			if (exists ${$Regions_r}[$region]{$f}{"ratio2center"}) {
				if (${$Regions_r}[$region]{$f}{"ratio2center"} > $maxYsup)
					{ $maxYsup = ${$Regions_r}[$region]{$f}{"ratio2center"}; }
				}
			}
		}
	$maxYsup *= $ploidy;
	$maxDepthGraph *= $ploidy;
	if (defined $maxDepthGraph && $maxYsup > $maxDepthGraph) { $maxYsup = $maxDepthGraph; }

	if ($ploidy == 1) {
		$cmdR .= "plot (c(0,0), xlim=c(0,$Nbr_Reg), ylim=c(0,$maxYsup), type =\"n\", main=\"".${$Result3_r}{$file}{$CNV}{"type"}.": $Chrom:".${$Regions_r}[$CNV]{"start"}."-".${$Regions_r}[$CNVend]{"end"}."\", xlab=\"\", ylab=\"depth_ratio_to_$center\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n";
		}
	else {
		$cmdR .= "plot (c(0,0), xlim=c(0,$Nbr_Reg), ylim=c(0,$maxYsup), type =\"n\", main=\"".${$Result3_r}{$file}{$CNV}{"type"}.": $Chrom:".${$Regions_r}[$CNV]{"start"}."-".${$Regions_r}[$CNVend]{"end"}."\", xlab=\"\", ylab=\"copy_nbr\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n";
		}

	##gene separations
	my$currentGene=""; my$tmpTxt=""; my$Nbr_gene=0;
	for my$r ($firstR..$lastR) {
		if (${$Regions_r}[$r]{"Gene"}) {
			if (${$Regions_r}[$r]{"Gene"} ne "NA" && ${$Regions_r}[$r]{"Gene"} ne $currentGene)  {
				$tmpTxt .= "abline(v=".($r-$firstR+0.5).", col=\"blue\", lty = \"dotted\", lwd=1)\n";
				$currentGene = ${$Regions_r}[$r]{"Gene"};
				$Nbr_gene++;
				}
			}
		}
	if ($Nbr_gene < ${$CNV_opt_r}{"maxGeneSep"}) { $cmdR .= $tmpTxt; }

	##region labels
	my@printReg=();
	if ($Nbr_Reg < ${$CNV_opt_r}{"maxGeneLab"}) {
		##in grey if invalid
		for my$r ($firstR..$lastR){ 
			if (defined ${$Regions_r}[$r]{"Appel"}) { push(@printReg,$r) ; }
			}
		if (@printReg) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r-$firstR+1).","; }
			chop $cmdR;
			$cmdR .="), labels=c(";
			foreach my$r (@printReg) { $cmdR .= "\"".${$Regions_r}[$r]{"label"}."\","; }
			chop $cmdR;
			$cmdR .= "), col.axis=\"darkgrey\", las=2";
			if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
			else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_Reg))." )\n"; }
			}
		##in black if valid
		@printReg=();
		for my$r ($firstR..$lastR) { 
			if(!defined ${$Regions_r}[$r]{"Appel"}) { push(@printReg,$r) ; }
			}
		if (@printReg) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r-$firstR+1).","; }
			chop $cmdR;
			$cmdR .="), labels=c(";
			foreach my$r (@printReg) { $cmdR .= "\"".${$Regions_r}[$r]{"label"}."\","; }
			chop $cmdR;
			$cmdR .= "), col.axis=\"black\", las=2";
			if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
			else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_Reg))." )\n"; }
			}
		}
	else {
		##in grey if invalid; only ticks
		for my$r ($firstR..$lastR) { 
			if (defined ${$Regions_r}[$r]{"Appel"}) { push(@printReg,$r) ; }
			}
		if (@printReg && scalar@printReg < ${$CNV_opt_r}{"maxGeneSep"}) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r-$firstR+1).","; }
			chop $cmdR;
			$cmdR .= "), labels = FALSE, col.ticks=\"darkgrey\")\n";
			}
		##in black if valid
		@printReg=();
		if ($Nbr_gene < ${$CNV_opt_r}{"maxGeneLab"}) {
			$currentGene="";
			for my$r ($firstR..$lastR) {
				if (${$Regions_r}[$r]{"Gene"}) {
					if (${$Regions_r}[$r]{"Gene"} ne "NA" && ${$Regions_r}[$r]{"Gene"} ne $currentGene)  { 
						push(@printReg,$r); 
						$currentGene = ${$Regions_r}[$r]{"Gene"};
						}
					}
				}
			}
		if (@printReg) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r-$firstR+1).","; }
			chop $cmdR;
			$cmdR .="), labels=c(";
			foreach my$r (@printReg) { $cmdR .= "\"".${$Regions_r}[$r]{"Gene"}."\","; }
			chop $cmdR;
			$cmdR .= "), col.axis=\"black\", las=2";
			if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
			else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_gene))." )\n"; }
			}
		##in red if CNV; only ticks
		@printReg=();
		for (my$r=0;$r<$Nbr_Reg;$r++) {
			if (exists ${$Result2_r}{$file}{$r}) { push(@printReg,$r) ; }
			}
		if (@printReg && scalar@printReg < ${$CNV_opt_r}{"maxGeneSep"}) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r-$firstR+1).","; }
			chop $cmdR;
			$cmdR .= "), labels = FALSE, col.ticks=\"red\")\n";
			}


		}

	##all not target sample lines (grey):
	for (my$f=0;$f<scalar@{$Files_r};$f++) {
		unless (${$Files_r}[$f] eq $file) {
			my$r1 = $firstR;
			while ($r1 <= $lastR) {
				my$r2=$r1;
				while ($r2 <= $lastR) {
					if (exists ${$Regions_r}[$r2]{${$Files_r}[$f]}{"ratio2center"}) { $r2++; }
					else { last; }
					}
				if (($r2-1) > $r1) {
					$cmdR .= "lines( c(";
					for (my$r=$r1;$r<$r2;$r++)
						{ $cmdR .= ($r-$firstR+1).","; }
					chop $cmdR;
					$cmdR .= "), c(";
					for (my$r=$r1;$r<$r2;$r++)
						{ $cmdR .= ($ploidy*${$Regions_r}[$r]{${$Files_r}[$f]}{"ratio2center"}).","; }
					chop $cmdR;
					$cmdR .= "), type =\"l\", lwd=1, col=\"darkgrey\")\n";
					}
				elsif (($r2-1) == $r1) {
					$cmdR .= "lines( c(".($r1-$firstR+1)."), c(".($ploidy*${$Regions_r}[$r1]{${$Files_r}[$f]}{"ratio2center"})."), type =\"p\", lwd=1, col=\"darkgrey\")\n";
					}
				$r1 = ($r2+1);
				}
			}
		}

	##threshold lines (black):
	if (${$CNV_opt_r}{"spread_test"}) {
		foreach my$gender ("normByR_depth","normByR_depth_fem","normByR_depth_males") {
			if ($Chrom =~ /^Y$/i && $gender eq "normByR_depth_fem") { next; }
			my$r1 = $firstR;
			while ($r1 <= $lastR) {
				my$r2=$r1;
				while ($r2 <= $lastR) {
					if (${$Regions_r}[$r2]{$gender}->{$center}) { $r2++;}
					else { last; }
					}
				if (($r2-1) > $r1) {
					foreach my$lim ("spread_sup","spread_inf") {
						$cmdR .= "lines( c(";
						for (my$r=$r1;$r<$r2;$r++)
							{ $cmdR .= ($r-$firstR+1).","; }
						chop $cmdR;
						$cmdR .= "), c(";
						for (my$r=$r1;$r<$r2;$r++) {
							if ( ${$Regions_r}[$r]{$gender}->{$center} )
								{ $cmdR .= ($ploidy*(${$Regions_r}[$r]{$gender}->{$lim} / ${$Regions_r}[$r]{$gender}->{$center})).","; }
							else 	{ $cmdR .= "$ploidy,"; }
							}
						chop $cmdR;
						$cmdR .= "), type =\"l\", lwd=1, col=\"black\")\n";
						}
					}
				elsif (($r2-1) == $r1) {
					if ( ${$Regions_r}[$r1]{$gender}->{$center} ) {
						$cmdR .= "lines( c(".($r1-$firstR+1)."), c(".($ploidy*(${$Regions_r}[$r1]{$gender}->{"spread_inf"} / ${$Regions_r}[$r1]{$gender}->{$center}))."), type =\"p\", lwd=1, col=\"black\")\n";
						$cmdR .= "lines( c(".($r1-$firstR+1)."), c(".($ploidy*(${$Regions_r}[$r1]{$gender}->{"spread_sup"} / ${$Regions_r}[$r1]{$gender}->{$center}))."), type =\"p\", lwd=1, col=\"black\")\n";
						}
					else 	{ $cmdR .= "$ploidy,"; }
					}
				$r1 = ($r2+1);
				}
			}
		}

	if (${$CNV_opt_r}{"center_test"}) {
		$cmdR .= "abline(h=".($ploidy*$seuil_del).", col=\"black\", lty = \"dashed\", lwd=1)\n";
		$cmdR .= "abline(h=".($ploidy*$seuil_dup).", col=\"black\", lty = \"dashed\", lwd=1)\n";
		}

	##target sample line (green)
	my$r1 = $firstR;
	while ($r1 <= $lastR) {
		my$r2=$r1;
		while ($r2 <= $lastR) {
			if (exists ${$Regions_r}[$r2]{$file}{"ratio2center"}) { $r2++;}
			else { last; }
			}
		if (($r2-1) > $r1) {
			$cmdR .= "lines( c(";
			for (my$r=$r1;$r<$r2;$r++)
				{ $cmdR .= ($r-$firstR+1).","; }
			chop $cmdR;
			$cmdR .= "), c(";
			for (my$r=$r1;$r<$r2;$r++)
				{ $cmdR .= ($ploidy*${$Regions_r}[$r]{$file}{"ratio2center"}).","; }
			chop $cmdR;
			$cmdR .= "), type =\"l\", lwd=1, col=\"green\")\n";
			}
		elsif (($r2-1) == $r1) {
			$cmdR .= "lines( c(".($r1-$firstR+1)."), c(".($ploidy*${$Regions_r}[$r1]{$file}{"ratio2center"})."), type =\"p\", lwd=1, col=\"green\")\n";
			}
		$r1 = ($r2+1);
		}

	##points for CNVs
	my$points .= "points( c(";
	for my$r ($firstR..$lastR) {
		if (exists ${$Result2_r}{$file}{$r} && exists ${$Regions_r}[$r]{$file}{"ratio2center"})
			{ $points .= ($r-$firstR+1).","; }
		}
	chop $points;
	$points .= "), c(";
	for my$r ($firstR..$lastR){
		if (exists ${$Result2_r}{$file}{$r} && exists ${$Regions_r}[$r]{$file}{"ratio2center"})
			{ $points .= ($ploidy*${$Regions_r}[$r]{$file}{"ratio2center"}).","; }
		}
	chop $points;
	$points .= "), type =\"p\", pch = 16, lwd=2, col=\"red\")\n";
	$cmdR .= $points;

	##center line
	$cmdR .= "abline(h=$ploidy, col=\"darkgrey\", lty = \"dashed\", lwd=1)\n";
	##center depth txt
	my$centerDepth = "";
	if (exists ${$Regions_r}[$CNV]{"normByR_depth"}) {
		$centerDepth = centerDepth($CNV,$CNVend,"normByR_depth",$center,$Regions_r);
		$cmdR .= "text(x = 0, y = $ploidy, labels = \"$center depth: $centerDepth X\", adj = c(0.2,-0.5), cex = ".(log(5)/log($Nbr_Reg)).", col = \"darkgray\")\n";
		}
	else {
		if ($Chrom !~ /^Y$/i) {
			$centerDepth = centerDepth($CNV,$CNVend,"normByR_depth_fem",$center,$Regions_r);
			$cmdR .= "text(x = 0, y = $ploidy, labels = \"$center depth (F): $centerDepth X\", adj = c(0.2,-0.5), cex = ".(log(5)/log($Nbr_Reg)).", col = \"darkgray\")\n";
			}
		$centerDepth = centerDepth($CNV,$CNVend,"normByR_depth_males",$center,$Regions_r);
		$cmdR .= "text(x = 0, y = $ploidy, labels = \"$center depth (M): $centerDepth X\", adj = c(0.2,1.5), cex = ".(log(5)/log($Nbr_Reg)).", col = \"darkgray\")\n";
		}

	print CMDR "$cmdR";
	
	}

##1 sheet / chr :
#print CMDR "title(main=\"sample: ${$sampleName_r}{$file}";
#if (${$PatientsRef}{$file}{"ecarte"}) { print CMDR " (invalid)\", col.main=\"red\""; }
#else { print CMDR "\""; }
#print CMDR ", outer=TRUE, line=-2, cex.main=2)\n";
print CMDR "dev.off()\nquit(save=\"no\")\n";
close CMDR;
system "Rscript $outdir/${$sampleName_r}{$file}\_temp.R";
unlink "$outdir/${$sampleName_r}{$file}\_temp.R";

}


####################
sub graphByCNV2 {

my($outdir,$CNV_opt_r,$patient,$PatientsRef,$Regions_r,$ChrOrderRef,$regionOrderRef,$Result2_r,$Result3_r)= @_;

my$center = ${$CNV_opt_r}{"center"};
my$spread = ${$CNV_opt_r}{"spread"};
my$seuil_del = ${$CNV_opt_r}{"seuil_del"};
my$seuil_dup = ${$CNV_opt_r}{"seuil_dup"};
my$spread_del = ${$CNV_opt_r}{"spread_del"};
my$spread_dup = ${$CNV_opt_r}{"spread_dup"};
my$ext = ${$CNV_opt_r}{"graphCNVpadding"};
my$maxDepthGraph = ${$CNV_opt_r}{"maxDepthGraph"};
my$ploidy = ${$CNV_opt_r}{"ploidy"};

print "\tcmdR, for ${$PatientsRef}{$patient}{ID}\n";
open (CMDR, ">$outdir/${$PatientsRef}{$patient}{ID}\_temp.R") || die;
print CMDR "#!/usr/bin/env Rscript\n\n" ;
print CMDR "pdf(\"".$outdir."/CNV_".${$PatientsRef}{$patient}{"ID"}.".pdf\", width=11.69, height=4.135)\n";

foreach my$CNV (sort{$a<=>$b}keys%{ ${$Result3_r}{$patient} }) {

	my$cmdR = "";

	my$Chrom = ${$Regions_r}[$CNV]{"Chrom"};
	my$ChrName = $Chrom; $ChrName =~ s/^chr//;

	my$firstR = $CNV;
	for (my$r=($CNV - $ext);$r<$CNV;$r++) {
		if (defined ${$Regions_r}[$r] && ${$Regions_r}[$r]{"Chrom"} eq $Chrom) {
			$firstR = $r; last;
			}
		}
	my$CNVend = ${$Result3_r}{$patient}{$CNV}{"end"};
	my$lastR = $CNVend;
	for (my$r=($CNVend + $ext);$r>=$CNVend;$r--) {
		if (defined ${$Regions_r}[$r] && ${$Regions_r}[$r]{"Chrom"} eq $Chrom) {
			$lastR = $r; last;
			}
		}
	my$Nbr_Reg = $lastR - $firstR + 1;

	my$maxYsup=$seuil_dup;
	for my$r ($firstR..$lastR) {
		for (my$p=0;$p<scalar(keys%{$PatientsRef});$p++) {
			if (exists ${$Regions_r}[$r]{$p}{"ratio2center"}) {
				if (${$Regions_r}[$r]{$p}{"ratio2center"} > $maxYsup)
					{ $maxYsup = ${$Regions_r}[$r]{$p}{"ratio2center"}; }
				}
			}
		}
	$maxYsup *= $ploidy;
	$maxDepthGraph *= $ploidy;
	if (defined $maxDepthGraph && $maxYsup > $maxDepthGraph) { $maxYsup = $maxDepthGraph; }

	if ($ploidy == 1) {
		$cmdR .= "plot (c(0,0), xlim=c(0,$Nbr_Reg), ylim=c(0,$maxYsup), type =\"n\", main=\"".${$Result3_r}{$patient}{$CNV}{"type"}.": $Chrom:".${$Regions_r}[$CNV]{"start"}."-".${$Regions_r}[$CNVend]{"end"}."\", xlab=\"\", ylab=\"depth_ratio_to_$center\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n";
		}
	else {
		$cmdR .= "plot (c(0,0), xlim=c(0,$Nbr_Reg), ylim=c(0,$maxYsup), type =\"n\", main=\"".${$Result3_r}{$patient}{$CNV}{"type"}.": $Chrom:".${$Regions_r}[$CNV]{"start"}."-".${$Regions_r}[$CNVend]{"end"}."\", xlab=\"\", ylab=\"copy_nbr\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n";
		}

	#gene vertical separations
	my$currentGene=""; my$tmpTxt=""; my$Nbr_gene=0;
	for my$r ($firstR..$lastR) {
		if (${$Regions_r}[$r]{"Gene"} ne "NA" && ${$Regions_r}[$r]{"geneID"} ne $currentGene) {
			$tmpTxt .= "abline(v=".($r-$firstR+0.5).", col=\"blue\", lty = \"dotted\", lwd=1)\n";
			$currentGene = ${$Regions_r}[$r]{"geneID"};
			$Nbr_gene++;
			}
		}
	if ($Nbr_gene < ${$CNV_opt_r}{"maxGeneSep"}) { $cmdR .= $tmpTxt; }

	#x labels
	my@printReg=();
	if ($Nbr_Reg < ${$CNV_opt_r}{"maxGeneLab"}) {
		##in grey if invalid
		for my$r ($firstR..$lastR) { 
			if (defined ${$Regions_r}[$r]{"Appel"}) { push(@printReg,$r) ; }
			}
		if (@printReg) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r-$firstR+1).","; }
			chop $cmdR;
			$cmdR .="), labels=c(";
			foreach my$r (@printReg) { $cmdR .= "\"".${$Regions_r}[$r]{"label"}."\","; }
			chop $cmdR;
			$cmdR .= "), col.axis=\"darkgrey\", las=2";
			if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
			else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_Reg))." )\n"; }
			}
		##in black if valid
		@printReg=();
		for my$r ($firstR..$lastR) { 
			if(!defined ${$Regions_r}[$r]{"Appel"}) { push(@printReg,$r) ; }
			}
		if (@printReg) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r-$firstR+1).","; }
			chop $cmdR;
			$cmdR .="), labels=c(";
			foreach my$r (@printReg) { $cmdR .= "\"".${$Regions_r}[$r]{"label"}."\","; }
			chop $cmdR;
			$cmdR .= "), col=\"black\", las=2";
			if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
			else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_Reg))." )\n"; }
			}
		}
	else {
		##in grey if invalid; only ticks
		for my$r ($firstR..$lastR) { 
			if (defined ${$Regions_r}[$r]{"Appel"}) { push(@printReg,$r) ; }
			}
		if (@printReg && scalar@printReg < ${$CNV_opt_r}{"maxGeneSep"}) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r-$firstR+1).","; }
			chop $cmdR;
			$cmdR .= "), labels = FALSE, col.ticks=\"darkgrey\")\n";
			}
		##in black if valid
		@printReg=();
		if ($Nbr_gene < ${$CNV_opt_r}{"maxGeneLab"}) {	# && ($Nbr_gene+$Nbr_CNV)>=${$CNV_opt_r}{"maxGeneLab"}) {
			$currentGene="";
			for my$r ($firstR..$lastR) {
				if (${$Regions_r}[$r]{"Gene"} ne "NA" && ${$Regions_r}[$r]{"geneID"} ne $currentGene) {
					push(@printReg,$r);
					$currentGene = ${$Regions_r}[$r]{"geneID"};
					}
				}
			}
		if (@printReg) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r-$firstR+1).","; }
			chop $cmdR;
			$cmdR .="), labels=c(";
			foreach my$r (@printReg) { $cmdR .= "\"".${$Regions_r}[$r]{"label"}."\","; }
			chop $cmdR;
			$cmdR .= "), col.axis=\"black\", las=2";
			if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
			else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_gene))." )\n"; }
			}
		##in red if CNV; only ticks
		@printReg=();
		for my$r ($firstR..$lastR)  {
			if (exists ${$Result2_r}{$patient}{$r}) { push(@printReg,$r) ; }
			}
		if (@printReg && scalar@printReg < ${$CNV_opt_r}{"maxGeneSep"}) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r-$firstR+1).","; }
			chop $cmdR;
			$cmdR .= "), labels = FALSE, col.ticks=\"red\")\n";
			}
		}

	#all not-target sample lines (grey):
	for (my$p=0;$p<scalar(keys%{$PatientsRef});$p++) {
		unless ($p == $patient) {
			my$r1 = $firstR;
			while ($r1 <= $lastR) {
				my$r2=$r1;
				while ($r2 <= $lastR) {
					if (exists ${$Regions_r}[$r2]{$p}{"ratio2center"}) { $r2++;}
					else { last; }
					}
				if (($r2-1) > $r1) {
					$cmdR .= "lines( c(";
					for (my$r=$r1;$r<$r2;$r++)
						{ $cmdR .= ($r-$firstR+1).","; }
					chop $cmdR;
					$cmdR .= "), c(";
					for (my$r=$r1;$r<$r2;$r++)
						{ $cmdR .= ($ploidy*${$Regions_r}[$r]{$p}{"ratio2center"}).","; }
					chop $cmdR;
					$cmdR .= "), type =\"l\", lwd=1, col=\"darkgrey\")\n";
					}
				elsif (($r2-1) == $r1) {
					$cmdR .= "lines( c(".($r1-$firstR+1)."), c(".($ploidy*${$Regions_r}[$r1]{$p}{"ratio2center"})."), type =\"p\", lwd=1, col=\"darkgrey\")\n";
					}
				$r1 = ($r2+1);
				}
			}
		}

	##threshold lines (black):
	if (${$CNV_opt_r}{"spread_test"}) {
		foreach my$gender ("normByR_depth","normByR_depth_fem","normByR_depth_males") {
			if ($ChrName =~ /^Y$/i && $gender eq "normByR_depth_fem") { next; }
			my$r1 = $firstR;
			while ($r1 <= $lastR) {
				my$r2=$r1;
				while ($r2 <= $lastR) {
					if (${$Regions_r}[$r2]{$gender}->{$center}) { $r2++; }
					else { last; }
					}
				if (($r2-1) > $r1) {
					foreach my$lim ("spread_sup","spread_inf") {
						$cmdR .= "lines( c(";
						for (my$r=$r1;$r<$r2;$r++)
							{ $cmdR .= ($r-$firstR+1).","; }
						chop $cmdR;
						$cmdR .= "), c(";
						for (my$r=$r1;$r<$r2;$r++) {
							if ( ${$Regions_r}[$r]{$gender}->{$center} )
								{ $cmdR .= ($ploidy*(${$Regions_r}[$r]{$gender}->{$lim} / ${$Regions_r}[$r]{$gender}->{$center})).","; }
							else 	{ $cmdR .= "$ploidy,"; }
							}
						chop $cmdR;
						$cmdR .= "), type =\"l\", lwd=1, col=\"black\")\n";
						}
					}
				elsif (($r2-1) == $r1) {
					if ( ${$Regions_r}[$r1]{$gender}->{$center} ) {
						$cmdR .= "lines( c(".($r1-$firstR+1)."), c(".($ploidy*(${$Regions_r}[$r1]{$gender}->{"spread_inf"} / ${$Regions_r}[$r1]{$gender}->{$center}))."), type =\"p\", lwd=1, col=\"black\")\n";
						$cmdR .= "lines( c(".($r1-$firstR+1)."), c(".($ploidy*(${$Regions_r}[$r1]{$gender}->{"spread_sup"} / ${$Regions_r}[$r1]{$gender}->{$center}))."), type =\"p\", lwd=1, col=\"black\")\n";
						}
					else 	{ $cmdR .= "$ploidy,"; }
					}
				$r1 = ($r2+1);
				}
			}
		}
	if (${$CNV_opt_r}{"center_test"}) {
		$cmdR .= "abline(h=".($ploidy*$seuil_del).", col=\"black\", lty = \"dashed\", lwd=1)\n";
		$cmdR .= "abline(h=".($ploidy*$seuil_dup).", col=\"black\", lty = \"dashed\", lwd=1)\n";
		}

	#target sample line (green):
	my$r1 = $firstR;
	while ($r1 <= $lastR) {
		my$r2=$r1;
		while ($r2 <= $lastR) {
			if (exists ${$Regions_r}[$r2]{$patient}{"ratio2center"}) { $r2++;}
			else { last; }
			}
		if (($r2-1) > $r1) {
			$cmdR .= "lines( c(";
			for (my$r=$r1;$r<$r2;$r++)
				{ $cmdR .= ($r-$firstR+1).","; }
			chop $cmdR;
			$cmdR .= "), c(";
			for (my$r=$r1;$r<$r2;$r++)
				{ $cmdR .= ($ploidy*${$Regions_r}[$r]{$patient}{"ratio2center"}).","; }
			chop $cmdR;
			$cmdR .= "), type =\"l\", lwd=1, col=\"green\")\n";
			}
		elsif (($r2-1) == $r1) {
			$cmdR .= "lines( c(".($r1-$firstR+1)."), c(".($ploidy*${$Regions_r}[$r1]{$patient}{"ratio2center"})."), type =\"p\", lwd=1, col=\"green\")\n";
			}
		$r1 = ($r2+1);
		}

	#points for CNVs
	my$points .= "points( c(";
	for my$r ($firstR..$lastR) {
		if (exists ${$Result2_r}{$patient}{$r} && exists ${$Regions_r}[$r]{$patient}{"ratio2center"})
			{ $points .= ($r-$firstR+1).","; }
		}
	chop $points;
	$points .= "), c(";
	for my$r ($firstR..$lastR) {
		if (exists ${$Result2_r}{$patient}{$r} && exists ${$Regions_r}[$r]{$patient}{"ratio2center"})
			{ $points .= ($ploidy*${$Regions_r}[$r]{$patient}{"ratio2center"}).","; }
		}
	chop $points;
	$points .= "), type =\"p\", pch = 16, lwd=2, col=\"red\")\n";
	$cmdR .= $points;

	## center line
	$cmdR .= "abline(h=$ploidy, col=\"darkgrey\", lty = \"dashed\", lwd=1)\n";

	##center depth txt
	my$centerDepth;
	if (exists ${$Regions_r}[$CNV]{"normByR_depth"}->{$center}) {
		$centerDepth = centerDepth($CNV,$CNVend,"normByR_depth",$center,$Regions_r);
		$cmdR .= "text(x = 0, y = $ploidy, labels = \"$center depth: $centerDepth X\", adj = c(0.2,-0.5), cex = ".(log(5)/log($Nbr_Reg)).", col = \"darkgray\")\n";
		}
	else {
		if ($ChrName !~ /^Y$/i) {
			$centerDepth = centerDepth($CNV,$CNVend,"normByR_depth_fem",$center,$Regions_r);
			$cmdR .= "text(x = 0, y = $ploidy, labels = \"$center depth (F): $centerDepth X\", adj = c(0.2,-0.5), cex = ".(log(5)/log($Nbr_Reg)).", col = \"darkgray\")\n";
			}
		$centerDepth = centerDepth($CNV,$CNVend,"normByR_depth_males",$center,$Regions_r);
		$cmdR .= "text(x = 0, y = $ploidy, labels = \"$center depth (M): $centerDepth X\", adj = c(0.2,1.5), cex = ".(log(5)/log($Nbr_Reg)).", col = \"darkgray\")\n";
		}

	print CMDR "$cmdR";

	}

#print CMDR "title(main=\"sample: ${$PatientsRef}{$patient}{ID}";
#if (${$PatientsRef}{$patient}{"ecarte"}) { print CMDR " (invalid)\", col.main=\"red\""; }
#else { print CMDR "\""; }
#print CMDR ", outer=TRUE, line=-2, cex.main=2)\n";
print CMDR "dev.off()\nquit(save=\"no\")\n";
close CMDR;
system "Rscript $outdir/${$PatientsRef}{$patient}{ID}\_temp.R";
unlink "$outdir/${$PatientsRef}{$patient}{ID}\_temp.R";

}

####################

sub centerDepth {
my($start,$end,$gender,$center,$Regions_r) = @_;
my$centerDepth = 0;
my$n = 0;
for my$r ($start..$end)  {
	if (exists ${$Regions_r}[$r]{$gender}->{$center}) {
		$centerDepth += ${$Regions_r}[$r]{$gender}->{$center};
		$n++;
		}
	if ($n) { $centerDepth /= $n; }
	#else {print "$start\t$gender\t".${$Regions_r}[$start]{$gender}->{$center}."\n".${$Regions_r}[$start]{"Chrom"}.":".${$Regions_r}[$start]{"start"}."\n";}
	}
return(sprintf("%.1f",$centerDepth));
}

####################




1;

