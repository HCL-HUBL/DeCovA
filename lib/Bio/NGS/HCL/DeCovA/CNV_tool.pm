package Bio::NGS::HCL::DeCovA::CNV_tool;


use strict;
use warnings;


#@hashSub = Bio::NGS::HCL::DeCovA::CNV_tool::CNV_detect(\@Files,\%sName2,"$outdir/CNV_analysis",$nGraf,$fichier_sexe,\%CNV_opt,$mbq,$mmq,$refBedLines,\@ChromOrder,\%fai);
sub CNV_detect
{
my($FilesRef,$sampleNameRef,$outdir,$nGraf,$fichier_sexe,$CNV_optRef,$mbq,$mmq,$RegionsRef,$ChromOrderRef,$faiRef)=@_;
#my@Files = @$h1;
#my%sampleName = %$h2;
#my@Regions = @$h4;
#my@ChromOrder = @$h5;
#my%fai = %$h6;
#my%CNV_opt = %{$CNV_optRef};
my$norm = ${$CNV_optRef}{"norm"};
my$normByGender = ${$CNV_optRef}{"normByGender"};
my$Ref = ${$CNV_optRef}{"RefDepth"};
my$RefByGender = ${$CNV_optRef}{"RefByGender"};
my$seuil_region = ${$CNV_optRef}{"seuil_region"};
if ($seuil_region < (1/scalar@{$FilesRef})) { $seuil_region = (1/scalar@{$FilesRef}); }
my$seuil_patient = ${$CNV_optRef}{"seuil_patient"};
if ($seuil_patient < (1/scalar@{$RegionsRef})) { $seuil_patient = (1/scalar@{$RegionsRef}); }
my$minCov = ${$CNV_optRef}{"seuil_cov"};
my$seuil_deletion = ${$CNV_optRef}{"seuil_deletion"};
my$seuil_duplication = ${$CNV_optRef}{"seuil_duplication"};
my$minCNV = ${$CNV_optRef}{"min_following_CNV"};
my$minDP = ${$CNV_optRef}{"min_DP"} ;
my$maxNonCNV = ${$CNV_optRef}{"max_Non_CNV"};
my$maxNonCNVrate = ${$CNV_optRef}{"max_Non_CNV_rate"};
my$graphByChr = ${$CNV_optRef}{"chromGraph"};
my@cnvFields = @{ ${$CNV_optRef}{"fields"} };
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

my$nCol = scalar(split(/\t/,${$RegionsRef}[0]{"allLine"}));

my%Patients;

#read gender file
if ($fichier_sexe) {
	open(PATIENTS, "$fichier_sexe") or die "Fichier $fichier_sexe impossible a lire\n";
	foreach my $ligne (<PATIENTS>) {
		$ligne =~ m/(\S+)\s+(\S+)/;
		my$ok="";
		foreach my$file (@{$FilesRef}) {
			if ($1 eq ${$sampleNameRef}{$file}) {
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
foreach my$file (@{$FilesRef}) {
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
for (my $r = 0 ; $r < scalar@{$RegionsRef} ; $r++) {
	my $region_a_conserver = 1;
	if ($minCov) {
		foreach my$file (@{$FilesRef}) {
			if ( ${$RegionsRef}[$r]{$file}{"Cov"} >= $minCov)
				{ $region_a_conserver = 0; last; }
			}
		}
	if ($minDP) {
		my@allDepth = ();
		foreach my$file (@{$FilesRef}) {
			push(@allDepth, ${$RegionsRef}[$r]{$file}{"Mean"});
			}
		my$normDepth = norm($norm,\@allDepth);
		if ( $normDepth < $minDP)
			{ $region_a_conserver = 0; }
		}
	# si region à conserver, calcul Profondeur_Autosomes/gonosomes
	if ($region_a_conserver) {
		foreach my$file (@{$FilesRef}) {
			if (${$RegionsRef}[$r]{"Chrom"} !~ /^(chr[XY]|[XY])$/) {
				$autosom_Regions++;
				if ($Ref eq "mean") { $Patients{$file}{"tot_Autosomes_depth"} += ${$RegionsRef}[$r]{$file}{"Mean"}; }
				}
			else {
				if (${$RegionsRef}[$r]{"Chrom"} =~ /^(chrX|X)$/) {
					$gonosom_Regions++;
					if ($Ref eq "mean") { $Patients{$file}{"tot_chrX_depth"} += ${$RegionsRef}[$r]{$file}{"Mean"}; }
					}
				else {
					if ($Ref eq "mean") { $Patients{$file}{"tot_chrY_depth"} += ${$RegionsRef}[$r]{$file}{"Mean"}; }
					}
				}
			}
		}
	# Si region mal couverte, on l'imprime dans un fichier POUBELLE
	else
		{ ${$RegionsRef}[$r]{"Appel"} = "Couverture_Faible"; }
	}

my$depthTxt="";
my$sexTxt="";
foreach my$file (@{$FilesRef}) {

	##if tot bases
	if ($Ref eq "tot") {
		#autosomes
		my$cmd = "samtools view -F 0x4";
		if ($mmq) { $cmd .= " -q $mmq"; }
		$cmd .= " $file";
		foreach my$chr (@{ ${$faiRef}{$file} }) {
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
		foreach my$chr (@{ ${$faiRef}{$file} }) {
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
		foreach my$chr (@{ ${$faiRef}{$file} }) {
			if ($chr =~ /^(chrY|Y)$/)
				{ $cmd .= " $chr"; }
			}
		$cmd .= " | perl -ne 'END { print \"\$cnt\" } \@fields = split (/\\t/,\$_); \@bases = split (//,\$fields[9]);";
		if ($mbq) { $cmd .= " \@quals = map { unpack(\"C*\", \$_ )} split(//,\$fields[10]); for \$i (0..\$#bases) { if ((\$bases[\$i] =~ /[ACGT]/i) && (\$quals[\$i] >= $mbq)) { \$cnt++; } }'"; }
		else { $cmd .= " for \$i (0..\$#bases) { if (\$bases[\$i] =~ /[ACGT]/i) { \$cnt++; } }'"; }
		print "$cmd\n";
		$Patients{$file}{"tot_chrY_depth"} = `$cmd`;
		chomp $Patients{$file}{"tot_chrY_depth"};
		
		$depthTxt .= ${$sampleNameRef}{$file}.":\n\ttotal sequenced bases for ${$sampleNameRef}{$file} : \n\tautoZ:".$Patients{$file}{"tot_Autosomes_depth"}."\n\tchrX:".$Patients{$file}{"tot_chrX_depth"}."\n\tchrY:".$Patients{$file}{"tot_chrY_depth"}."\n";
		}

	##gender if no sex file
	$sexTxt .= ${$sampleNameRef}{$file}.":\n";
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
foreach my$file (@{$FilesRef}) {
	if ($Patients{$file}{"Ref_depth"})
		{ $meanRef += $Patients{$file}{"Ref_depth"}; }
	else { $Patients{$file}{"ecarte"} = 1; }
	}
$meanRef /= scalar@{$FilesRef};
foreach my$file (@{$FilesRef}) {
	$Patients{$file}{"Ref_depth"} /= $meanRef; 
	print "normalization factor for ${$sampleNameRef}{$file} : ".$Patients{$file}{"Ref_depth"}."\n";
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
	foreach my$file (@{$FilesRef})
		{ print SORTIE ${$sampleNameRef}{$file}."\t"; }
	print SORTIE "Statut_Region\n";


	# PREMIER PARCOURS BIS #
	# RECALCUL DE LA PROFONDEUR TOTALE POUR LA PONDERATION INTRA
	# CE (RE)CALCUL A LIEU SI DES REGIONS ETAIENT "MOCHES" APRES L'APPEL DE CNV
	if ($nb_parcours > 0) {

		if ($Ref eq "mean") {
		# REINITIALISATION DE LA PROFONDEUR TOTALE PAR PATIENT
			foreach my$file (@{$FilesRef}) {
				$Patients{$file}{"tot_Autosomes_depth"} = 0;
				$Patients{$file}{"Ref_Gonosomes_Patient"} = 0;
				}
			for (my $r = 0 ; $r < scalar@{$RegionsRef} ; $r++) {
				if(!defined ${$RegionsRef}[$r]{"Appel"}) {
					foreach my$file (@{$FilesRef}) {
						if (${$RegionsRef}[$r]{"Chrom"} !~ /^(chr[XY]|[XY])$/) {
							$Patients{$file}{"tot_Autosomes_depth"} += ${$RegionsRef}[$r]{$file}{"Mean"};
							}
						elsif (${$RegionsRef}[$r]{"Chrom"} =~ /^(chrX|X)$/) {
							$Patients{$file}{"Ref_Gonosomes_Patient"} += ${$RegionsRef}[$r]{$file}{"Mean"};
							}
						else {
							if ((!$RefByGender) && ($Patients{$file}{"Sexe"} eq "H")) {
								$Patients{$file}{"Ref_Gonosomes_Patient"} += ${$RegionsRef}[$r]{$file}{"Mean"};
								}
							}
						}
					}
				}
			
			foreach my$file (@{$FilesRef}) {
				if ($RefByGender && ($Patients{$file}{"Sexe"} eq "H")) {
					$Patients{$file}{"Ref_depth"} = ($Patients{$file}{"Ref_Gonosomes_Patient"} * 2) + $Patients{$file}{"tot_Autosomes_depth"}; } 
				else {
					$Patients{$file}{"Ref_depth"} = $Patients{$file}{"Ref_Gonosomes_Patient"} + $Patients{$file}{"tot_Autosomes_depth"}; }
				}
			foreach my$file (@{$FilesRef})
				{ $meanRef += $Patients{$file}{"Ref_depth"}; }
			$meanRef /= scalar@{$FilesRef};
			foreach my$file (@{$FilesRef}) {
				$Patients{$file}{"Ref_depth"} /= $meanRef;
				}
			}
		}
	print LOG "sample normalization factors :\n";
	foreach my$file (@{$FilesRef}) {
		print LOG "\t${$sampleNameRef}{$file} : $Patients{$file}{Ref_depth}\n";
	}
	print LOG "\n";


	# SECOND PARCOURS DES REGIONS #
	# PERMET DE PONDERER LA PROFONDEUR PAR LES AUTRES REGIONS (INTRA) ET ENTRE LES PATIENTS (INTER)
	for (my $r = 0 ; $r < scalar@{$RegionsRef} ; $r++) {

		print SORTIE ${$RegionsRef}[$r]{"allLine"}."\t";
		print SORTIE $r."\t";

		# SI LA REGION N'EST PAS "MOCHE"
		if (!defined ${$RegionsRef}[$r]{"Appel"}) {

			$nb_regions_conservees++;
			${$RegionsRef}[$r]{"nb_CNV"}{"DEL"} = 0;
			${$RegionsRef}[$r]{"nb_CNV"}{"DUP"} = 0;
			my $nb_evts = 0;
			my$prof_Moyenne_Inter = "";

			if ($normByGender) {

				@{ ${$RegionsRef}[$r]{"all_normByS_depths_fem"} } = ();
				@{ ${$RegionsRef}[$r]{"all_normByS_depths_males"} } = ();
				# Pour chaque patient, un premier parcours permet la pondération intra-patient (en divisant la profondeur de la region par la profondeur totale du patient)
				foreach my$file (@{$FilesRef}) {

					if ($Patients{$file}{"Sexe"} eq "F" && !($Patients{$file}{"ecarte"})) {

						# Le controle permet de diviser uniquement par le nombre de patients restant apres l filtre, et non par l'ensemble des patients
						# Si la region concerne un gonosome et le patient est une femme, la référence est celle qui correspond aux autosomes + au chromosome X
						${$RegionsRef}[$r]{$file}{"normByS_depth"} = ${$RegionsRef}[$r]{$file}{"Mean"} / $Patients{$file}{"Ref_depth"};
						# On sauvegarde la profondeur ponderee intra patient, pour determiner ensuite la moyenne de ces profondeurs ponderees
						push( @{ ${$RegionsRef}[$r]{"all_normByS_depths_fem"} }, ${$RegionsRef}[$r]{$file}{"normByS_depth"} );


					} elsif ($Patients{$file}{"Sexe"} eq "H" && !($Patients{$file}{"ecarte"})) {

						# Le controle permet de divisier uniquement par le nombre de patients restant apres l filtre, et non par l'ensemble des patients
						# Si la region concerne un gonosome et le patient est un homme, la référence est celle qui correspond uniquement au chromosome X
						${$RegionsRef}[$r]{$file}{"normByS_depth"} = ${$RegionsRef}[$r]{$file}{"Mean"} / $Patients{$file}{"Ref_depth"};
						# On sauvegarde la profondeur ponderee intra patient, pour determiner ensuite la moyenne de ces profondeurs ponderees
						push( @{ ${$RegionsRef}[$r]{"all_normByS_depths_males"} }, ${$RegionsRef}[$r]{$file}{"normByS_depth"} );

					}

				}

				# Nous calculons pour chaque region et chaque sexe la moyenne de profondeur (ponderee pour chaque patient par rapport aux autres regions)
				if ($normByGender eq "all") {
					if(@{ ${$RegionsRef}[$r]{"all_normByS_depths_fem"} }) {
						foreach (@cnvVal) 
							{ ${$RegionsRef}[$r]{"normByR_depth_fem"}{$_} = norm($_,\@{ ${$RegionsRef}[$r]{"all_normByS_depths_fem"} }); }
						if ($norm eq "std") {
							my$sqtotal = 0;
							foreach (@{ ${$RegionsRef}[$r]{"all_normByS_depths_fem"} })
								{ $sqtotal += (${$RegionsRef}[$r]{"normByR_depth_fem"}{"moy"}-$_)**2; }
							${$RegionsRef}[$r]{"normByR_depth_fem"}{"std"} = ($sqtotal / (scalar@{ ${$RegionsRef}[$r]{"all_normByS_depths_fem"} }-1))**0.5;
						}
					}
					if(@{ ${$RegionsRef}[$r]{"all_normByS_depths_males"} }) {
						foreach (@cnvVal)
							{ ${$RegionsRef}[$r]{"normByR_depth_males"}{$_} = norm({$_},\@{ ${$RegionsRef}[$r]{"all_normByS_depths_males"} }); }
						if ($norm eq "std") {
							my$sqtotal = 0;
							foreach (@{ ${$RegionsRef}[$r]{"all_normByS_depths_males"} })
								{ $sqtotal += (${$RegionsRef}[$r]{"normByR_depth_males"}{"moy"}-$_)**2; }
							${$RegionsRef}[$r]{"normByR_depth_males"}{"std"} = ($sqtotal / (scalar@{ ${$RegionsRef}[$r]{"all_normByS_depths_males"} }-1))**0.5;
						}
					}


				} else {
					if (${$RegionsRef}[$r]{"Chrom"} !~ /^(chr[XY]|[XY])$/) {
						my@Autosomes=();
						if(@{ ${$RegionsRef}[$r]{"all_normByS_depths_fem"} }) {
							push(@Autosomes, @{ ${$RegionsRef}[$r]{"all_normByS_depths_fem"} });
						}
						if(@{ ${$RegionsRef}[$r]{"all_normByS_depths_males"} }) {
							push(@Autosomes, @{ ${$RegionsRef}[$r]{"all_normByS_depths_males"} });
						}
						if (@Autosomes) {
							foreach (@cnvVal) {
								${$RegionsRef}[$r]{"normByR_depth"}{$_} = norm($_,\@Autosomes);
								}
							if ($norm eq "std") {
								my$sqtotal = 0;
								foreach (@Autosomes)
									{ $sqtotal += (${$RegionsRef}[$r]{"normByR_depth"}{"moy"}-$_)**2; }
								${$RegionsRef}[$r]{"normByR_depth"}{"std"} = ($sqtotal / (scalar@Autosomes-1))**0.5;
							}
						}
					} else {
						if(@{ ${$RegionsRef}[$r]{"all_normByS_depths_fem"} }) {
							foreach (@cnvVal) {
								${$RegionsRef}[$r]{"normByR_depth_fem"}{$_} = norm($_,\@{ ${$RegionsRef}[$r]{"all_normByS_depths_fem"} });
							}
							if ($norm eq "std") {
								my$sqtotal = 0;
								foreach (@{ ${$RegionsRef}[$r]{"all_normByS_depths_fem"} })
									{ $sqtotal += (${$RegionsRef}[$r]{"normByR_depth_fem"}{"moy"}-$_)**2; }
								${$RegionsRef}[$r]{"normByR_depth_fem"}{"std"} = ($sqtotal / (scalar@{ ${$RegionsRef}[$r]{"all_normByS_depths_fem"} }-1))**0.5;
							}
						}
						if(@{ ${$RegionsRef}[$r]{"all_normByS_depths_males"} }) {
							foreach (@cnvVal) {
								${$RegionsRef}[$r]{"normByR_depth_males"}{$_} = norm($_,\@{ ${$RegionsRef}[$r]{"all_normByS_depths_males"} });
							}
							if ($norm eq "std") {
								my$sqtotal = 0;
								foreach (@{ ${$RegionsRef}[$r]{"all_normByS_depths_males"} })
									{ $sqtotal += (${$RegionsRef}[$r]{"normByR_depth_males"}{"moy"}-$_)**2; }
								${$RegionsRef}[$r]{"normByR_depth_males"}{"std"} = ($sqtotal / (scalar@{ ${$RegionsRef}[$r]{"all_normByS_depths_males"} }-1))**0.5;
							}
						}
					}
				}

			
				# Un second parcours permet l'appel de CNV (en comparant la profondeur pour un patient par la profondeur moyenne de la region)
				foreach my$file (@{$FilesRef}) {

					if ($norm eq "std") { $prof_Moyenne_Inter = 0; }
					else  { $prof_Moyenne_Inter = 1; }

					if(! $Patients{$file}{"ecarte"}) {

						if ($normByGender eq "all") {

							if($Patients{$file}{"Sexe"} eq "F") {
								if (${$RegionsRef}[$r]{"normByR_depth_fem"}{$norm}) {
									if (${$RegionsRef}[$r]{"Chrom"} !~ /^(chrY|Y)$/) {
										if ($norm eq "std") {
											${$RegionsRef}[$r]{$file}{"depth_ratio"}{"moy"} = ${$RegionsRef}[$r]{$file}{"normByS_depth"}/${$RegionsRef}[$r]{"normByR_depth_fem"}{"moy"};
											${$RegionsRef}[$r]{$file}{"depth_ratio"}{"std"} = (${$RegionsRef}[$r]{$file}{"normByS_depth"}-${$RegionsRef}[$r]{"normByR_depth_fem"}{"moy"})/${$RegionsRef}[$r]{"normByR_depth_fem"}{"std"};
										} else {
											${$RegionsRef}[$r]{$file}{"depth_ratio"}{$norm} = ${$RegionsRef}[$r]{$file}{"normByS_depth"}/${$RegionsRef}[$r]{"normByR_depth_fem"}{$norm};
										}
										$prof_Moyenne_Inter = ${$RegionsRef}[$r]{$file}{"depth_ratio"}{$norm};
										print SORTIE "$prof_Moyenne_Inter\t";
									} else {
										# En faisant la distinction entre chromosomes X et Y, les femmes ne sont plus considérées pour l'appel de CNV des régions du dernier
										print SORTIE "NA\t";
									}
								} else {
									${$RegionsRef}[$r]{"Appel"} = "No_Data";
									print SORTIE "NA\t";
								}

							} elsif ($Patients{$file}{"Sexe"} eq "H") {
								if (${$RegionsRef}[$r]{"normByR_depth_males"}{$norm}) {
									if ($norm eq "std") {
										${$RegionsRef}[$r]{$file}{"depth_ratio"}{"moy"} = ${$RegionsRef}[$r]{$file}{"normByS_depth"}/${$RegionsRef}[$r]{"normByR_depth_males"}{"moy"};
										${$RegionsRef}[$r]{$file}{"depth_ratio"}{"std"} = (${$RegionsRef}[$r]{$file}{"normByS_depth"}-${$RegionsRef}[$r]{"normByR_depth_males"}{"moy"})/${$RegionsRef}[$r]{"normByR_depth_males"}{"std"};
									} else {
										${$RegionsRef}[$r]{$file}{"depth_ratio"}{$norm} = ${$RegionsRef}[$r]{$file}{"normByS_depth"}/${$RegionsRef}[$r]{"normByR_depth_males"}{$norm};
									}
									$prof_Moyenne_Inter = ${$RegionsRef}[$r]{$file}{"depth_ratio"}{$norm};
									print SORTIE "$prof_Moyenne_Inter\t";
								} else {
									${$RegionsRef}[$r]{"Appel"} = "No_Data";
									print SORTIE "NA\t";
								}
							}

						} else {
							# Pour les autosomes
							if (${$RegionsRef}[$r]{"Chrom"} !~ /^(chr[XY]|[XY])$/) {
								if (${$RegionsRef}[$r]{"normByR_depth"}{$norm}) {
									if ($norm eq "std") {
										${$RegionsRef}[$r]{$file}{"depth_ratio"}{"moy"} = ${$RegionsRef}[$r]{$file}{"normByS_depth"}/${$RegionsRef}[$r]{"normByR_depth"}{"moy"};
										${$RegionsRef}[$r]{$file}{"depth_ratio"}{"std"} = (${$RegionsRef}[$r]{$file}{"normByS_depth"}-${$RegionsRef}[$r]{"normByR_depth"}{"moy"})/${$RegionsRef}[$r]{"normByR_depth"}{"std"};
									} else {
										${$RegionsRef}[$r]{$file}{"depth_ratio"}{$norm} = ${$RegionsRef}[$r]{$file}{"normByS_depth"}/${$RegionsRef}[$r]{"normByR_depth"}{$norm};
									}
									$prof_Moyenne_Inter = ${$RegionsRef}[$r]{$file}{"depth_ratio"}{$norm};
									print SORTIE "$prof_Moyenne_Inter\t";
								} else {
									${$RegionsRef}[$r]{"Appel"} = "No_Data";
									print SORTIE "NA\t";
								}
							} else {
							# Pour les gonosomes
								if($Patients{$file}{"Sexe"} eq "F") {
									if (${$RegionsRef}[$r]{"normByR_depth_fem"}{$norm}) {
										if (${$RegionsRef}[$r]{"Chrom"} =~ /^(chrX|X)$/) {
											if ($norm eq "std") {
												${$RegionsRef}[$r]{$file}{"depth_ratio"}{"moy"} = ${$RegionsRef}[$r]{$file}{"normByS_depth"}/${$RegionsRef}[$r]{"normByR_depth_fem"}{"moy"};
												${$RegionsRef}[$r]{$file}{"depth_ratio"}{"std"} = (${$RegionsRef}[$r]{$file}{"normByS_depth"}-${$RegionsRef}[$r]{"normByR_depth_fem"}{"moy"})/${$RegionsRef}[$r]{"normByR_depth_fem"}{"std"};
											} else {
												${$RegionsRef}[$r]{$file}{"depth_ratio"}{$norm} = ${$RegionsRef}[$r]{$file}{"normByS_depth"}/${$RegionsRef}[$r]{"normByR_depth_fem"}{$norm};
											}
											$prof_Moyenne_Inter = ${$RegionsRef}[$r]{$file}{"depth_ratio"}{$norm};
											print SORTIE "$prof_Moyenne_Inter\t";
										} else {
											# En faisant la distinction entre chromosomes X et Y, les femmes ne sont plus considérées pour l'appel de CNV des régions du dernier
											print SORTIE "NA\t";
										}
									} else {
										${$RegionsRef}[$r]{"Appel"} = "No_Data";
										print SORTIE "NA\t";
									}
								} elsif ($Patients{$file}{"Sexe"} eq "H") {
									if (${$RegionsRef}[$r]{"normByR_depth_males"}{$norm}) {
										if ($norm eq "std") {
											${$RegionsRef}[$r]{$file}{"depth_ratio"}{"moy"} = ${$RegionsRef}[$r]{$file}{"normByS_depth"}/${$RegionsRef}[$r]{"normByR_depth_males"}{"moy"};
											${$RegionsRef}[$r]{$file}{"depth_ratio"}{"std"} = (${$RegionsRef}[$r]{$file}{"normByS_depth"}-${$RegionsRef}[$r]{"normByR_depth_males"}{"moy"})/${$RegionsRef}[$r]{"normByR_depth_males"}{"std"};
										} else {
											${$RegionsRef}[$r]{$file}{"depth_ratio"}{$norm} = ${$RegionsRef}[$r]{$file}{"normByS_depth"}/${$RegionsRef}[$r]{"normByR_depth_males"}{$norm};
										}
										$prof_Moyenne_Inter = ${$RegionsRef}[$r]{$file}{"depth_ratio"}{$norm};
										print SORTIE "$prof_Moyenne_Inter\t";
									} else {
										${$RegionsRef}[$r]{"Appel"} = "No_Data";
										print SORTIE "NA\t";
									}
								}
							}
						}
			
						#print $Regions[$r]{$patient}{"depth_ratio"}."\t";
						if ($prof_Moyenne_Inter < $seuil_deletion) {
							${$RegionsRef}[$r]{"nb_CNV"}{"DEL"}++;
							$nb_evts++;
							$CNV{$file}++;
							$Results{$file}{$r} = "DEL";
						} elsif ($prof_Moyenne_Inter > $seuil_duplication) {
							${$RegionsRef}[$r]{"nb_CNV"}{"DUP"}++;
							$nb_evts++;
							$CNV{$file}++;
							$Results{$file}{$r} = "DUP";
						}

					} else {
						print SORTIE "NA\t";
					}
				}

				my$recurrent="";
				if(@{ ${$RegionsRef}[$r]{"all_normByS_depths_fem"} } && @{ ${$RegionsRef}[$r]{"all_normByS_depths_males"} } && (($nb_evts/(scalar@{ ${$RegionsRef}[$r]{"all_normByS_depths_fem"} } + scalar@{ ${$RegionsRef}[$r]{"all_normByS_depths_males"} })) > $seuil_region) && ($nb_parcours > 0)) { $recurrent = 1; }
				elsif(@{ ${$RegionsRef}[$r]{"all_normByS_depths_fem"} } && !@{ ${$RegionsRef}[$r]{"all_normByS_depths_males"} } && (($nb_evts/scalar@{ ${$RegionsRef}[$r]{"all_normByS_depths_fem"} }) > $seuil_region) && ($nb_parcours > 0)) { $recurrent = 1; }
				elsif(!@{ ${$RegionsRef}[$r]{"all_normByS_depths_fem"} } && @{ ${$RegionsRef}[$r]{"all_normByS_depths_males"} } && (($nb_evts/scalar@{ ${$RegionsRef}[$r]{"all_normByS_depths_males"} }) > $seuil_region) && ($nb_parcours > 0)) { $recurrent = 1; }
				if ($recurrent) {
					print SORTIE "CNV_Recurrent";
					print LOG "Region ecartee \: ". ${$RegionsRef}[$r]{"allLine"}."\n";
					${$RegionsRef}[$r]{"Appel"} = "CNV_Recurrent";
					$regions_ecartees++;
					$continuer = 1;

				} else {		
					print SORTIE "OK";
				}
				print SORTIE "\n";


			} else {

				@{ ${$RegionsRef}[$r]{"all_normByS_depths"} }=();
				# Pour chaque patient, un premier parcours permet la pondération intra-patient (en divisant la profondeur de la region par la profondeur totale recalculée du patient)
				foreach my$file (@{$FilesRef}) {
					if(! $Patients{$file}{"ecarte"}) {
						${$RegionsRef}[$r]{$file}{"normByS_depth"} = ${$RegionsRef}[$r]{$file}{"Mean"} / $Patients{$file}{"Ref_depth"};
						if ( ($Patients{$file}{"Sexe"} eq "H") && (${$RegionsRef}[$r]{"Chrom"} =~ /^(chrX|X)$/) )
							{ push(@{ ${$RegionsRef}[$r]{"all_normByS_depths"} }, (${$RegionsRef}[$r]{$file}{"normByS_depth"}*2)); }
						elsif ( ($Patients{$file}{"Sexe"} eq "F") && (${$RegionsRef}[$r]{"Chrom"} =~ /^(chrY|Y)$/) )
							{ next ; }
						else { push(@{ ${$RegionsRef}[$r]{"all_normByS_depths"} }, ${$RegionsRef}[$r]{$file}{"normByS_depth"}); }
					}
				}
				# Nous calculons pour chaque region la mediane/moyenne de profondeur (ponderee pour chaque patient par rapport aux autres regions)
				if(@{ ${$RegionsRef}[$r]{"all_normByS_depths"} }) {
					foreach (@cnvVal) { 
						${$RegionsRef}[$r]{"normByR_depth"}{$_} = norm($_,\@{ ${$RegionsRef}[$r]{"all_normByS_depths"} });
						}
					if ($norm eq "std") {
						my$sqtotal = 0;
						foreach (@{ ${$RegionsRef}[$r]{"all_normByS_depths"} })
							{ $sqtotal += (${$RegionsRef}[$r]{"normByR_depth"}{"moy"}-$_)**2; }
						${$RegionsRef}[$r]{"normByR_depth"}{"std"} = ($sqtotal / (scalar@{ ${$RegionsRef}[$r]{"all_normByS_depths"} }-1))**0.5;
					}
				}


				# Un second parcours permet l'appel de CNV (en comparant la profondeur pour un patient par la profondeur moyenne de la region)
				foreach my$file (@{$FilesRef}) {

					if ($norm eq "std") { $prof_Moyenne_Inter = 0; }
					else  { $prof_Moyenne_Inter = 1; }

					if(! $Patients{$file}{"ecarte"}) {

						if (${$RegionsRef}[$r]{"normByR_depth"}{$norm}) {
							if($Patients{$file}{"Sexe"} eq "F") {
								if (${$RegionsRef}[$r]{"Chrom"} !~ /^(chrY|Y)$/) {
									if ($norm eq "std") {
										${$RegionsRef}[$r]{$file}{"depth_ratio"}{"moy"} = ${$RegionsRef}[$r]{$file}{"normByS_depth"}/${$RegionsRef}[$r]{"normByR_depth"}{"moy"};
										${$RegionsRef}[$r]{$file}{"depth_ratio"}{"std"} = (${$RegionsRef}[$r]{$file}{"normByS_depth"}-${$RegionsRef}[$r]{"normByR_depth"}{"moy"})/${$RegionsRef}[$r]{"normByR_depth"}{"std"};
									} else {
										${$RegionsRef}[$r]{$file}{"depth_ratio"}{$norm} = ${$RegionsRef}[$r]{$file}{"normByS_depth"}/${$RegionsRef}[$r]{"normByR_depth"}{$norm};
									}
									$prof_Moyenne_Inter = ${$RegionsRef}[$r]{$file}{"depth_ratio"}{$norm};
									print SORTIE "$prof_Moyenne_Inter\t";
								} else {
									#les femmes ne sont pas considérées pour l'appel de CNV des régions du Y
									print SORTIE "NA\t";
								}
							} else {
								#*2 for chrX
								if (${$RegionsRef}[$r]{"Chrom"} =~ /^(chrX|X)$/) {
									if ($norm eq "std") {
										${$RegionsRef}[$r]{$file}{"depth_ratio"}{"moy"} = ${$RegionsRef}[$r]{$file}{"normByS_depth"}*2/${$RegionsRef}[$r]{"normByR_depth"}{"moy"};
										${$RegionsRef}[$r]{$file}{"depth_ratio"}{"std"} = (${$RegionsRef}[$r]{$file}{"normByS_depth"}*2-${$RegionsRef}[$r]{"normByR_depth"}{"moy"})/${$RegionsRef}[$r]{"normByR_depth"}{"std"};
									} else {
										${$RegionsRef}[$r]{$file}{"depth_ratio"}{$norm} = ${$RegionsRef}[$r]{$file}{"normByS_depth"}*2/${$RegionsRef}[$r]{"normByR_depth"}{$norm};
									}
								} else {
									if ($norm eq "std") {
										${$RegionsRef}[$r]{$file}{"depth_ratio"}{"moy"} = ${$RegionsRef}[$r]{$file}{"normByS_depth"}/${$RegionsRef}[$r]{"normByR_depth"}{"moy"};
										${$RegionsRef}[$r]{$file}{"depth_ratio"}{"std"} = (${$RegionsRef}[$r]{$file}{"normByS_depth"}-${$RegionsRef}[$r]{"normByR_depth"}{"moy"})/${$RegionsRef}[$r]{"normByR_depth"}{"std"};
									} else {
										${$RegionsRef}[$r]{$file}{"depth_ratio"}{$norm} = ${$RegionsRef}[$r]{$file}{"normByS_depth"}/${$RegionsRef}[$r]{"normByR_depth"}{$norm};
									}
								}
								$prof_Moyenne_Inter = ${$RegionsRef}[$r]{$file}{"depth_ratio"}{$norm};
								print SORTIE "$prof_Moyenne_Inter\t";
							}
							if ($prof_Moyenne_Inter < $seuil_deletion) {
								${$RegionsRef}[$r]{"nb_CNV"}{"DEL"}++;
								$nb_evts++;
								$CNV{$file}++;
								$Results{$file}{$r} = "DEL";
							} elsif ($prof_Moyenne_Inter > $seuil_duplication) {
								${$RegionsRef}[$r]{"nb_CNV"}{"DUP"}++;
								$nb_evts++;
								$CNV{$file}++;
								$Results{$file}{$r} = "DUP";
							}
						} else {
							${$RegionsRef}[$r]{"Appel"} = "No_Data";
							#$prof_Moyenne_Inter = 1;
							print SORTIE "NA\t";
						}
					} else {
						print SORTIE "NA\t";
					}

				}

				if(@{ ${$RegionsRef}[$r]{"all_normByS_depths"} } && (($nb_evts/scalar@{ ${$RegionsRef}[$r]{"all_normByS_depths"} }) > $seuil_region) && ($nb_parcours > 0)) {
					print SORTIE "CNV_Recurrent";
					print LOG "Region ecartee \: ". ${$RegionsRef}[$r]{"allLine"}."\n";
					${$RegionsRef}[$r]{"Appel"} = "CNV_Recurrent";
					$regions_ecartees++;
					$continuer = 1;

				} else {		
					print SORTIE "OK";
				}
				print SORTIE "\n";
			}

		# SI LA REGION EST "MOCHE"
		} else {
			foreach my$file (@{$FilesRef}) {
				print SORTIE "NA\t";
			}
			print SORTIE ${$RegionsRef}[$r]{"Appel"}."\n";
		}

	}

	print LOG "Nombre de regions ecartees lors de cette iteration \: ".$regions_ecartees."\n\n\n";	

	my $patients_ecartes = 0;

	foreach my$file (@{$FilesRef}) {

		if (! $Patients{$file}{"ecarte"}) {
			my $prcent_evt;
			if (defined($CNV{$file})) {

				$prcent_evt = $CNV{$file}/$nb_regions_conservees;

				if ($prcent_evt >= $seuil_patient) {
					print LOG "Patient ecarte \: ".${$sampleNameRef}{$file}."\n";
					$Patients{$file}{"ecarte"} = 1;
					$continuer = 1;
					$patients_ecartes++;			
				} else {
					$Patients{$file}{"ecarte"} = 0;
				}


			} else {
				print ${$sampleNameRef}{$file}." \: aucun CNV identifie\n";
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


for my$file (@{$FilesRef}) {
	mkdir "$outdir/${$sampleNameRef}{$file}";
	}


##print CNV foreach sample

## $Results{$file}{$r} = "DEL";
## @{ $RegionOrder{$Chrom} } = [ regions sorted by pos ]
my(%uniqReg,%RegInArray,%regionOrder,%regionIndice);
for (my $r = 0 ; $r < scalar@{$RegionsRef} ; $r++) {
	unless (exists $uniqReg{ ${$RegionsRef}[$r]{"Chrom"}."-".${$RegionsRef}[$r]{"Start"}."-".${$RegionsRef}[$r]{"End"} }) {
		push (@{ $RegInArray{ ${$RegionsRef}[$r]{"Chrom"} } }, $r);
		$uniqReg{ ${$RegionsRef}[$r]{"Chrom"}."-".${$RegionsRef}[$r]{"Start"}."-".${$RegionsRef}[$r]{"End"} } = 1;
		}
	my@tab = split(/\t/, ${$RegionsRef}[$r]{"allLine"});
	if ($tab[3]) {
		my@tab2 = split(/:|,/, $tab[3]);
		${$RegionsRef}[$r]{"label"} = substr($tab2[0], 0, 25);
		${$RegionsRef}[$r]{"Gene"} = $tab2[0];
		}
	else { ${$RegionsRef}[$r]{"label"} = "$tab[1]-$tab[2]"; }
	}
foreach my$Chrom (keys%RegInArray) {
	##sort by region ends (in case several regions with same start) then by starts
	@{ $regionOrder{$Chrom} } = sort{${$RegionsRef}[$a]{"End"}<=>${$RegionsRef}[$b]{"End"}}@{ $RegInArray{$Chrom} }; 
	@{ $regionOrder{$Chrom} } = sort{${$RegionsRef}[$a]{"Start"}<=>${$RegionsRef}[$b]{"Start"}}@{ $RegInArray{$Chrom} }; 
	for (my$i=0;$i<scalar@{ $regionOrder{$Chrom} };$i++) { $regionIndice{ $regionOrder{$Chrom}[$i] } = $i; }
	}

## @{ $orderedCNV{$patient}{$Chrom} } = [r1, r2,...]
my%orderedCNV;
foreach my$file (keys%Results) {
	my(%uniqCNV,%CNVinArray);
	foreach my$r (sort{$a<=>$b}keys%{ $Results{$file} }) {
		if (!exists $uniqCNV{${$RegionsRef}[$r]{"Chrom"}."-".${$RegionsRef}[$r]{"Start"}."-".${$RegionsRef}[$r]{"End"}}) { 
			push (@{ $CNVinArray{ ${$RegionsRef}[$r]{"Chrom"} } }, $r);
			$uniqCNV{ ${$RegionsRef}[$r]{"Chrom"}."-".${$RegionsRef}[$r]{"Start"}."-".${$RegionsRef}[$r]{"End"} } = 1;
			}
		}
	foreach my$Chrom (keys%CNVinArray) {
		##sort by region ends (in case several regions with same start) then by starts
		@{ $orderedCNV{$file}{$Chrom} } = sort{${$RegionsRef}[$a]{"End"}<=>${$RegionsRef}[$b]{"End"}}@{ $CNVinArray{$Chrom} }; 
		@{ $orderedCNV{$file}{$Chrom} } = sort{${$RegionsRef}[$a]{"Start"}<=>${$RegionsRef}[$b]{"Start"}}@{ $CNVinArray{$Chrom} }; 
		}
	}

unless ($minCNV) { $minCNV = 1; }
my(%Result2,%Result3,%Result4);
foreach my$file (keys%Results) {
	open (CNV1,">$outdir/${$sampleNameRef}{$file}/CNV_${$sampleNameRef}{$file}.summary.txt") or die "Pb lors de l'ecriture du fichier sortie $!\n";
	print CNV1 "#Region\tCNV\tclean_intervals_Nbr\tclean_ratio_to_$norm\tdirty_intervals_Nbr\tdirty_ratio_to_$norm\toverlapping_Samples\n";
	open (CNV2,">$outdir/${$sampleNameRef}{$file}/CNV_${$sampleNameRef}{$file}.allIntervals.txt") or die "Pb lors de l'ecriture du fichier sortie $!\n";
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
	foreach my$Chromosome (@{$ChromOrderRef}) {
		my$Chrom = $Chromosome ; $Chrom =~ s/chr//i;
		if (exists $orderedCNV{$file}{$Chrom}) {
			my$r=0;		##index in @{ $orderedCNV{$file}{$Chrom} }
			while ($r < scalar@{ $orderedCNV{$file}{$Chrom} } ) {
				##merge consecutive CNVs
				my@refSub = mergeConsecutiveCNV("",$maxNonCNV,$orderedCNV{$file}{$Chrom}[$r],$regionIndice{$orderedCNV{$file}{$Chrom}[$r]},\@{ $regionOrder{$Chrom} },\%{ $Results{$file} },$RegionsRef);
				my$i=$refSub[0]; my$cnvOK=$refSub[1]; my$nonCNVtot=$refSub[2];
				my@nextReg = @{ $refSub[3] };
				if ($maxNonCNVrate) {
					while (($nonCNVtot/($cnvOK+1)) > $maxNonCNVrate) {
						@refSub = mergeConsecutiveCNV(($cnvOK-1),$maxNonCNV,$orderedCNV{$file}{$Chrom}[$r],$regionIndice{$orderedCNV{$file}{$Chrom}[$r]},\@{ $regionOrder{$Chrom} },\%{ $Results{$file} },$RegionsRef);
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
						$Result3{$file}{$nextReg[0]}{"type"} = $Results{$file}{$nextReg[0]};
						$Result3{$file}{$nextReg[0]}{"end"} = $nextReg[0];
						$Result4{$file}{$Chrom}{${$RegionsRef}[$nextReg[0]]{"Start"}}{${$RegionsRef}[$nextReg[0]]{"End"}} = $Results{$file}{$nextReg[0]};
						if (exists $regionIndice{$nextReg[0]}) {
							print CNV1 $Chromosome.":".${$RegionsRef}[$nextReg[0]]{"Start"}."-".${$RegionsRef}[$nextReg[0]]{"End"}."\t".$Results{$file}{$nextReg[0]}."\t1\t".sprintf("%.3f",${$RegionsRef}[$nextReg[0]]{$file}{"depth_ratio"}{$norm})."\t.\t.\t$overlapCNT\n";
							print CNV2 $Chromosome."\t".${$RegionsRef}[$nextReg[0]]{"Start"}."\t".${$RegionsRef}[$nextReg[0]]{"End"}."\t".(${$RegionsRef}[$nextReg[0]]{"End"}-${$RegionsRef}[$nextReg[0]]{"Start"}+1)." bp\t".${$RegionsRef}[$nextReg[0]]{"label"}."\t".($regionIndice{$nextReg[0]}+1)."\t".$Results{$file}{$nextReg[0]}."\t".sprintf("%.3f",${$RegionsRef}[$nextReg[0]]{$file}{"depth_ratio"}{$norm})."\t".${$RegionsRef}[$nextReg[0]]{"nb_CNV"}{$Results{$file}{$nextReg[0]}};
							if (@cnvFields) {
								my$txt = printCNVfields($normByGender,\@cnvFields,\%{ ${$RegionsRef}[$nextReg[0]] },\%{ $Patients{$file} });
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
								$Result4{$file}{$Chrom}{${$RegionsRef}[$nextReg[$j]]{"Start"}}{${$RegionsRef}[$nextReg[$j]]{"End"}} = $Results{$file}{$nextReg[$j]};
								push(@cleanCNV, ${$RegionsRef}[$nextReg[$j]]{$file}{"depth_ratio"}{$norm});
								push(@dirtyCNV, ${$RegionsRef}[$nextReg[$j]]{$file}{"depth_ratio"}{$norm});
								if (exists $regionIndice{$nextReg[$j]}) {
									print CNV2 $Chromosome."\t".${$RegionsRef}[$nextReg[$j]]{"Start"}."\t".${$RegionsRef}[$nextReg[$j]]{"End"}."\t".(${$RegionsRef}[$nextReg[$j]]{"End"}-${$RegionsRef}[$nextReg[$j]]{"Start"}+1)." bp\t".${$RegionsRef}[$nextReg[$j]]{"label"}."\t".($regionIndice{$nextReg[$j]}+1)."\t".$Results{$file}{$nextReg[$j]}."\t".sprintf("%.3f",${$RegionsRef}[$nextReg[$j]]{$file}{"depth_ratio"}{$norm})."\t".${$RegionsRef}[$nextReg[$j]]{"nb_CNV"}{$Results{$file}{$nextReg[$j]}};
									if (@cnvFields) {
										my$txt = printCNVfields($normByGender,\@cnvFields,\%{ ${$RegionsRef}[$nextReg[$j]] },\%{ $Patients{$file} });
										print CNV2 "$txt";
										}
									print CNV2 "\n";
									}
								}
							else {
								$Result2{$file}{$nextReg[$j]} = "NA";
								$Result4{$file}{$Chrom}{${$RegionsRef}[$nextReg[$j]]{"Start"}}{${$RegionsRef}[$nextReg[$j]]{"End"}} = "NA";
								if (exists $regionIndice{$nextReg[$j]}) {
									print CNV2 $Chromosome."\t".${$RegionsRef}[$nextReg[$j]]{"Start"}."\t".${$RegionsRef}[$nextReg[$j]]{"End"}."\t".(${$RegionsRef}[$nextReg[$j]]{"End"}-${$RegionsRef}[$nextReg[$j]]{"Start"}+1)." bp\t".${$RegionsRef}[$nextReg[$j]]{"label"}."\t".($regionIndice{$nextReg[$j]}+1)."\t";
									if (exists ${$RegionsRef}[$nextReg[$j]]{"Appel"}) { print CNV2 "NA\t"; }
									else { print CNV2 "no\t"; }
									if (exists ${$RegionsRef}[$nextReg[$j]]{$file}{"depth_ratio"}{$norm}) {
										print CNV2 sprintf("%.3f",${$RegionsRef}[$nextReg[$j]]{$file}{"depth_ratio"}{$norm})."\t".${$RegionsRef}[$nextReg[$j]]{"nb_CNV"}{$Results{$file}{$nextReg[0]}};
										push(@dirtyCNV, ${$RegionsRef}[$nextReg[$j]]{$file}{"depth_ratio"}{$norm});
										}
									else { print CNV2 "na\tna"; }
									if (@cnvFields) {
										my$txt = printCNVfields($normByGender,\@cnvFields,\%{ ${$RegionsRef}[$nextReg[$j]] },\%{ $Patients{$file} });
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
						print CNV1 $Chromosome.":".${$RegionsRef}[$nextReg[0]]{"Start"}."-".${$RegionsRef}[$nextReg[$cnvOK]]{"End"}."\t".$Results{$file}{$nextReg[0]}."\t".scalar@cleanCNV."\t".sprintf("%.3f",$cleanAverage)."\t".($cnvOK+1)."\t".sprintf("%.3f",$dirtyAverage)."\t$overlapCNT\n";
						$Result3{$file}{$nextReg[0]}{"type"} = $Results{$file}{$nextReg[0]};
						$Result3{$file}{$nextReg[0]}{"end"} = $nextReg[$cnvOK];
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


##print graph foreach Chrom/CNV
if ($graphByChr) {
	##graph by Chrom/CNV ?
	my$Nbr_Reg_max = 0;
	foreach my$Chrom (keys%regionOrder) {
		if (scalar@{ $regionOrder{$Chrom} } > $Nbr_Reg_max) {
			$Nbr_Reg_max = scalar@{ $regionOrder{$Chrom} };
			}
		}

	if ($Nbr_Reg_max > ${$CNV_optRef}{"switch2graphByCNV"}) {
		for my$file (@{$FilesRef}) {
			graphByCNV1("$outdir/".${$sampleNameRef}{$file},$CNV_optRef,$file,$FilesRef,$sampleNameRef,$RegionsRef,$ChromOrderRef,\%regionOrder,\%Patients,\%Result2,\%Result3);
			}
		}
	else {
		for my$file (@{$FilesRef}) {
			graphByChr1($nGraf,"$outdir/".${$sampleNameRef}{$file},$CNV_optRef,$file,$FilesRef,$sampleNameRef,$RegionsRef,$ChromOrderRef,\%regionOrder,\%Patients,\%Result2);
			}
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

print OUT "\nintervals nber: ".scalar@{$RegionsRef}."\n";
print OUT "\nintervals discarded:\n";
my$N_noData=0; my$N_CNV_Recurrent=0; my$N_lowCov=0;
for (my $r = 0 ; $r < scalar@{$RegionsRef} ; $r++) {
	if(defined ${$RegionsRef}[$r]{"Appel"}) {
		if (${$RegionsRef}[$r]{"Appel"} eq "No_Data") { $N_noData++; }
		elsif (${$RegionsRef}[$r]{"Appel"} eq "CNV_Recurrent") { $N_CNV_Recurrent++; }
		elsif (${$RegionsRef}[$r]{"Appel"} eq "Couverture_Faible")  { $N_lowCov++; }
		}
	}
print OUT "\tnot enough data: $N_noData\n";
print OUT "\tCNV_Recurrent: $N_CNV_Recurrent\n";
if ($minDP || $minCov) { print OUT "\tlow coverage :  $N_lowCov\n"; }

print OUT "\nsamples discarded:\n";
my$N_discard=0;
foreach my$file (@{$FilesRef}) {
		if ($Patients{$file}{"ecarte"}) { print OUT "\t${$sampleNameRef}{$file}\n"; $N_discard++; }
		}
unless ($N_discard) { print OUT "\tnone\n"; }

print OUT "\nResults :\n
patient\tCNV\tnber\n";
foreach my$file (@{$FilesRef}) {
	unless($Patients{$file}{"ecarte"}) {
		my$N_dup=0; my$N_del=0;
		foreach (keys%{ $Result3{$file} }) {
			if ($Result3{$file}{$_}{"type"} eq "DUP") { $N_dup++; }
			elsif ($Result3{$file}{$_}{"type"} eq "DEL") { $N_del++; }
			}
		print OUT "\n${$sampleNameRef}{$file}: 
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
my($fichier_cov,$outdir,$nGraf,$fichier_sexe,$CNV_optRef)=@_;
#my%CNV_opt = %$h1;
my$norm = ${$CNV_optRef}{"norm"};
my$normByGender = ${$CNV_optRef}{"normByGender"};
my$Ref = ${$CNV_optRef}{"RefDepth"};
my$RefByGender = ${$CNV_optRef}{"RefByGender"};
my$seuil_region = ${$CNV_optRef}{"seuil_region"};
my$seuil_patient = ${$CNV_optRef}{"seuil_patient"};
my$seuil_cov = ${$CNV_optRef}{"seuil_cov"};
my$seuil_deletion = ${$CNV_optRef}{"seuil_deletion"};
my$seuil_duplication = ${$CNV_optRef}{"seuil_duplication"};
my$minCNV = ${$CNV_optRef}{"min_following_CNV"};
my$minDP = ${$CNV_optRef}{"min_DP"} ;
my$maxNonCNV = ${$CNV_optRef}{"max_Non_CNV"};
my$maxNonCNVrate = ${$CNV_optRef}{"max_Non_CNV_rate"};
my$graphByChr = ${$CNV_optRef}{"chromGraph"};
my@cnvFields = @{ ${$CNV_optRef}{"fields"} };
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
						$Result3{$patient}{$nextReg[0]}{"type"} = $Results{$patient}{$nextReg[0]};
						$Result3{$patient}{$nextReg[0]}{"end"} = $nextReg[0];
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
						$Result3{$patient}{$nextReg[0]}{"type"} = $Results{$patient}{$nextReg[0]};
						$Result3{$patient}{$nextReg[0]}{"end"} = $nextReg[$cnvOK];
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


##print graph foreach Chrom/CNV
if ($graphByChr) {
	##graph by Chrom/CNV ?
	my$Nbr_Reg_max = 0;
	foreach my$Chrom (@ChrOrder) {
		if (scalar@{ $regionOrder{$Chrom} } > $Nbr_Reg_max) {
			$Nbr_Reg_max = scalar@{ $regionOrder{$Chrom} };
			}
		}

	if ($Nbr_Reg_max > ${$CNV_optRef}{"switch2graphByCNV"}) {
		foreach my$patient (keys%Patients) {
			graphByCNV2($outdir,$CNV_optRef,$patient,\%Patients,\@Regions,\@ChrOrder,\%regionOrder,\%Result2,\%Result3);
			}
		}
	else {
		foreach my$patient (keys%Patients) {
			graphByChr2($outdir,$CNV_optRef,$patient,\%Patients,\@Regions,\@ChrOrder,\%regionOrder,\%Result2);
			}
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


}


####################

sub norm {
my($norm,$allDepthRef)=@_;
#my@allDepth=@$h1;
my(@sortDepth,$normDepth);
if ($norm eq "med" || $norm eq "min" || $norm eq "max") {
	@sortDepth = sort{$a<=>$b}@{$allDepthRef};
	if ($norm eq "min") { $normDepth = $sortDepth[0]; }
	elsif ($norm eq "max") { $normDepth = $sortDepth[-1]; }
	else {
		#odd?
		if(scalar@sortDepth%2) 
			{ $normDepth = $sortDepth[int(scalar@sortDepth/2)]; }
		#even
		else { $normDepth = ( $sortDepth[int(scalar@sortDepth/2)-1] + $sortDepth[int(scalar@sortDepth/2)] )/2; }
		}
	}

else {
	foreach (@{$allDepthRef})
		{ $normDepth += $_; } 
	$normDepth /= scalar@{$allDepthRef};
	}

return($normDepth);
}


####################

sub printCNVfields {
my($normByGender,$cnvFieldsRef,$RegionsRef,$PatientsRef) = @_;
#my@cnvFields = @$h1;
#my%Regions = %$h2;
#my%Patient = %$h3;
my$txt = "";
foreach my$val (@{$cnvFieldsRef}) {
	if (${$RegionsRef}{"normByR_depth"}{$val}) {
		if (!$normByGender || ($normByGender && ($normByGender eq "gono" && ${$RegionsRef}{"Chrom"} !~ m/^chr[XY]$|^[XY]$/))) {
			$txt .= "\t".sprintf("%.1f",${$RegionsRef}{"normByR_depth"}{$val});
			}
		else {
			if (${$PatientsRef}{"Sexe"} eq "F") {
				$txt .= "\t".sprintf("%.1f",${$RegionsRef}{"normByR_depth_fem"}{$val});
				}
			else {
				$txt .= "\t".sprintf("%.1f",${$RegionsRef}{"normByR_depth_males"}{$val});
				}
			}
		}
	else { $txt .= "\tna"; }
	}
return ($txt);
}


####################

sub mergeConsecutiveCNV {
my($maxI,$maxNonCNV,$orderedCNV,$regionIndice,$regionOrderRef,$ResultsRef,$RegionsRef)=@_;
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
			if (!exists ${$RegionsRef}[$nextReg[$i]]{"Appel"}) { 
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

my($nGraf,$outdir,$CNV_optRef,$file,$FilesRef,$sampleNameRef,$RegionsRef,$ChromOrderRef,$regionOrderRef,$PatientsRef,$ResultsRef)= @_;

print "cmdR, for sample ${$sampleNameRef}{$file}\n";

my$norm = ${$CNV_optRef}{"norm"};
my$normGraf = $norm;
if ($norm eq "std") { $normGraf = "moy"; }

my$maxDepthGraph = ${$CNV_optRef}{"maxDepthGraph"};
my$seuil_deletion = ${$CNV_optRef}{"seuil_deletion"};
my$seuil_duplication = ${$CNV_optRef}{"seuil_duplication"};

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
open (CMDR, ">$outdir/${$sampleNameRef}{$file}\_temp.R") || die;
print CMDR "#!/usr/bin/env Rscript\n\n" ;
print CMDR "pdf(\"$outdir/CNV_".${$sampleNameRef}{$file}.".pdf\", width=11.69, height=4,135)\n";

foreach my$Chrom (@{$ChromOrderRef}) {

	$Chrom =~ s/^chr//i;
	
	if (exists ${$regionOrderRef}{$Chrom}) {

		my$cmdR = "";

		my$Nbr_Reg = scalar@{ ${$regionOrderRef}{$Chrom} };
		##1 sheet / chr :
		my$maxX = $Nbr_Reg;

		my$maxYsup=$seuil_duplication; my$maxYinf=$seuil_deletion;	
		foreach my$region (@{ ${$regionOrderRef}{$Chrom} }) {
			foreach my$f (@{$FilesRef}) {
				if (exists ${$RegionsRef}[$region]{$f}{"depth_ratio"}{$normGraf}) {
					if (${$RegionsRef}[$region]{$f}{"depth_ratio"}{$normGraf} > $maxYsup)
						{ $maxYsup = ${$RegionsRef}[$region]{$f}{"depth_ratio"}{$normGraf}; }
					if ($normGraf eq "std") {
						if (${$RegionsRef}[$region]{$f}{"depth_ratio"}{$normGraf} < $maxYinf)
							{ $maxYinf = ${$RegionsRef}[$region]{$f}{"depth_ratio"}{$normGraf}; }
						}
					}
				}
			}
		if ($maxDepthGraph && $maxYsup > $maxDepthGraph) { $maxYsup = $maxDepthGraph; }
		if ($normGraf eq "std") {
			if ($maxDepthGraph && $maxYinf < (-$maxDepthGraph)) { $maxYinf = (-$maxDepthGraph); }
			}

		if ($normGraf eq "std") {
			##all in 1 sheet:
			#$cmdR .= "par(fig=c(0,1,".(1-(($n-0.05)/$nGraf)).",".(1-(($n-0.95)/$nGraf))."), new=TRUE)
#plot (c(0,0), xlim=c(0,$maxX), ylim=c($maxYinf,$maxYsup), type =\"n\", main=\"chrom: $Chrom\", xlab=\"\", ylab=\"depth_ratio_to_$normGraf\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n";
			##1 sheet / chr :
			$cmdR .= "plot (c(0,0), xlim=c(0,$maxX), ylim=c($maxYinf,$maxYsup), type =\"n\", main=\"chrom: $Chrom\", xlab=\"\", ylab=\"depth_ratio_to_$normGraf\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n";
			}
		else {
			##all in 1 sheet:
			#$cmdR .= "par(fig=c(0,1,".(1-(($n-0.05)/$nGraf)).",".(1-(($n-0.95)/$nGraf))."), new=TRUE)
#plot (c(0,0), xlim=c(0,$maxX), ylim=c(0,$maxYsup), type =\"n\", main=\"chrom: $Chrom\", xlab=\"\", ylab=\"depth_ratio_to_$normGraf\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n";
			##1 sheet / chr :
			$cmdR .= "plot (c(0,0), xlim=c(0,$maxX), ylim=c(0,$maxYsup), type =\"n\", main=\"chrom: $Chrom\", xlab=\"\", ylab=\"depth_ratio_to_$normGraf\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n";
			}

		##gene separations
		my$currentGene=""; my$tmpTxt=""; my$Nbr_gene=0;
		for (my$r=0;$r<$Nbr_Reg;$r++) {
			if (${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"}) {
				if (${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"} ne "NA" && ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"} ne $currentGene)  {
					$tmpTxt .= "abline(v=".($r+0.5).", col=\"blue\", lty = \"dotted\", lwd=1)\n";
					$currentGene = ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"};
					$Nbr_gene++;
					}
				}
			}
		if ($Nbr_gene < ${$CNV_optRef}{"maxGeneSep"}) { $cmdR .= $tmpTxt; }

		my$Nbr_CNV=0;
		for (my$r=0;$r<$Nbr_Reg;$r++) {
			if (exists ${$ResultsRef}{$file}{${$regionOrderRef}{$Chrom}[$r]})
				{ $Nbr_CNV++; }
			}

		##region labels
		my@printReg=();
		if ($maxX < ${$CNV_optRef}{"maxGeneLab"}) {
			##in grey if invalid
			for (my$r=0;$r<$Nbr_Reg;$r++) { 
				if (defined ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"Appel"}) { push(@printReg,$r) ; }
				}
			if (@printReg) {
				$cmdR .= "axis(1, at=c(";
				foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
				chop $cmdR;
				$cmdR .="), labels=c(";
				foreach my$r (@printReg) { $cmdR .= "\"".${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"label"}."\","; }
				chop $cmdR;
				$cmdR .= "), col.axis=\"darkgrey\", las=2";
				if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
				else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_Reg))." )\n"; }
				}
			##in black if valid
			@printReg=();
			for (my$r=0;$r<$Nbr_Reg;$r++) { 
				if(!defined ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"Appel"}) { push(@printReg,$r) ; }
				}
			if (@printReg) {
				$cmdR .= "axis(1, at=c(";
				foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
				chop $cmdR;
				$cmdR .="), labels=c(";
				foreach my$r (@printReg) { $cmdR .= "\"".${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"label"}."\","; }
				chop $cmdR;
				$cmdR .= "), col.axis=\"black\", las=2";
				if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
				else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_Reg))." )\n"; }
				}
			}
		else {
			##in grey if invalid; only ticks
			for (my$r=0;$r<$Nbr_Reg;$r++) { 
				if (defined ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"Appel"}) { push(@printReg,$r) ; }
				}
			if (@printReg && scalar@printReg < ${$CNV_optRef}{"maxGeneSep"}) {
				$cmdR .= "axis(1, at=c(";
				foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
				chop $cmdR;
				$cmdR .= "), labels = FALSE, col.ticks=\"darkgrey\")\n";
				}
			##in black if valid
			@printReg=();
			if ($Nbr_gene < ${$CNV_optRef}{"maxGeneLab"}) { #&& ($Nbr_gene+$Nbr_CNV)>=${$CNV_optRef}{"maxGeneLab"}) {
				$currentGene="";
				for (my$r=0;$r<$Nbr_Reg;$r++) {
					if (${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"}) {
						if (${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"} ne "NA" && ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"} ne $currentGene)  { 
							push(@printReg,$r); 
							$currentGene = ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"};
							}
						}
					}
				}
			#elsif ($Nbr_gene<${$CNV_optRef}{"maxGeneLab"} && ($Nbr_gene+$Nbr_CNV)<${$CNV_optRef}{"maxGeneLab"}) {
			#	for (my$r=0;$r<$Nbr_Reg;$r++) {
			#		if ( (${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"} && ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"} ne "NA" && ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"} ne $currentGene) || (exists ${$ResultsRef}{$file}{${$regionOrderRef}{$Chrom}[$r]}) ) {
			#			push(@printReg,$r);
			#			$currentGene = ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"};
			#			}
			#		}
			#	}
			if (@printReg) {
				$cmdR .= "axis(1, at=c(";
				foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
				chop $cmdR;
				$cmdR .="), labels=c(";
				foreach my$r (@printReg) { $cmdR .= "\"".${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"}."\","; }
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
			if (@printReg && scalar@printReg < ${$CNV_optRef}{"maxGeneSep"}) {
				$cmdR .= "axis(1, at=c(";
				foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
				chop $cmdR;
				$cmdR .= "), labels = FALSE, col.ticks=\"red\")\n";
				}
			
			
			}

		##all not target sample lines (grey):
		for (my$f=0;$f<scalar@{$FilesRef};$f++) {
			unless (${$FilesRef}[$f] eq $file) {
				my$r1=0;
				while ($r1<$Nbr_Reg) {
					my$r2=$r1;
					while ($r2<$Nbr_Reg) {
						if (exists ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r2]]{${$FilesRef}[$f]}{"depth_ratio"}{$normGraf}) { $r2++;}
						else { last; }
						}
					if (($r2-1) > $r1) {
						$cmdR .= "lines( c(";
						for (my$r=$r1;$r<$r2;$r++)
							{ $cmdR .= ($r+1).","; }
						chop $cmdR;
						$cmdR .= "), c(";
						for (my$r=$r1;$r<$r2;$r++)
							{ $cmdR .= ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{${$FilesRef}[$f]}{"depth_ratio"}{$normGraf}.","; }
						chop $cmdR;
						$cmdR .= "), type =\"l\", lwd=1, col=\"darkgrey\")\n";
						}
					elsif (($r2-1) == $r1) {
						$cmdR .= "lines( c(".($r1+1)."), c(".${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r1]]{${$FilesRef}[$f]}{"depth_ratio"}{$normGraf}."), type =\"p\", lwd=1, col=\"darkgrey\")\n";
						}
					$r1 = ($r2+1);
					}
				}
			}

		##threshold lines (black):
		if ($norm eq "std" && $normGraf eq "moy") {
			foreach my$gender ("normByR_depth","normByR_depth_fem","normByR_depth_males") {
				my$r1=0;
				while ($r1<$Nbr_Reg) {
					my$r2=$r1;
					while ($r2<$Nbr_Reg) {
						if (${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r2]]{$gender}{"moy"}) { $r2++;}
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
								if ( ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{$gender}{"moy"} )
									{ $cmdR .= (1+$fold*(${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{$gender}{"std"} / ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{$gender}{"moy"})).","; }
								else 	{ $cmdR .= "1,"; }
								}
							chop $cmdR;
							$cmdR .= "), type =\"l\", lwd=1, col=\"black\")\n";
							}
						}
					elsif (($r2-1) == $r1) {
						$cmdR .= "lines( c(".($r1+1)."), c(".(1+(${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r1]]{$gender}{"std"} / ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r1]]{$gender}{"moy"}))."), type =\"p\", lwd=1, col=\"black\")\n";
						$cmdR .= "lines( c(".($r1+1)."), c(".(1-(${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r1]]{$gender}{"std"} / ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r1]]{$gender}{"moy"}))."), type =\"p\", lwd=1, col=\"black\")\n";
						}
					$r1 = ($r2+1);
					}
				}
			}
		else {
			$cmdR .= "abline(h=$seuil_deletion, col=\"black\", lty = \"dashed\", lwd=1)\n";
			$cmdR .= "abline(h=$seuil_duplication, col=\"black\", lty = \"dashed\", lwd=1)\n";
			}

		##target sample line (green)
		my$r1=0;
		while ($r1<$Nbr_Reg) {
			my$r2=$r1;
			while ($r2<$Nbr_Reg) {
				if (exists ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r2]]{$file}{"depth_ratio"}{$normGraf}) { $r2++;}
				else { last; }
				}
			if (($r2-1) > $r1) {
				$cmdR .= "lines( c(";
				for (my$r=$r1;$r<$r2;$r++)
					{ $cmdR .= ($r+1).","; }
				chop $cmdR;
				$cmdR .= "), c(";
				for (my$r=$r1;$r<$r2;$r++)
					{ $cmdR .= ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{$file}{"depth_ratio"}{$normGraf}.","; }
				chop $cmdR;
				$cmdR .= "), type =\"l\", lwd=1, col=\"green\")\n";
				}
			elsif (($r2-1) == $r1) {
				$cmdR .= "lines( c(".($r1+1)."), c(".${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r1]]{$file}{"depth_ratio"}{$normGraf}."), type =\"p\", lwd=1, col=\"green\")\n";
				}
			$r1 = ($r2+1);
			}
		##points for CNVs
		my$points .= "points( c(";
		for (my$r=0;$r<$Nbr_Reg;$r++) {
			if (exists ${$ResultsRef}{$file}{${$regionOrderRef}{$Chrom}[$r]} && exists ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{$file}{"depth_ratio"}{$normGraf})
				{ $points .= ($r+1).","; }
			}
		chop $points;
		$points .= "), c(";
		foreach my$region (@{ ${$regionOrderRef}{$Chrom} }) {
			if (exists ${$ResultsRef}{$file}{$region} && exists ${$RegionsRef}[$region]{$file}{"depth_ratio"}{$normGraf})
				{ $points .= ${$RegionsRef}[$region]{$file}{"depth_ratio"}{$normGraf}.","; }
			}
		chop $points;
		$points .= "), type =\"p\", pch = 16, lwd=2, col=\"red\")\n";
		if ($Nbr_CNV) { $cmdR .= $points; }

		#if ($normGraf eq "std") { $cmdR .= "abline(h=0, col=\"darkgrey\", lty = \"dashed\", lwd=1)\n"; }
		#else { $cmdR .= "abline(h=1, col=\"darkgrey\", lty = \"dashed\", lwd=1)\n"; }
		$cmdR .= "abline(h=1, col=\"darkgrey\", lty = \"dashed\", lwd=1)\n";

		##all in 1 sheet:
		#if ($c==$Nbr_Chr || $n==$nGraf) {
		#	open (CMDR, ">$outdir/${$sampleNameRef}{$file}\_temp.R") || die;
		#	print CMDR "#!/usr/bin/env Rscript\n\n" ;
		#	if ($nGraf==$Nbr_Chr) { print CMDR "pdf(\"".$outdir."/CNV_${$sampleNameRef}{$file}.pdf\", width=11.69, height=".($nGraf*3).")\npar(mfrow=c($nGraf,1))\n"; }
		#	else {
		#		if ($N>1) { print CMDR "pdf(\"".$outdir."/CNV_${$sampleNameRef}{$file}\_$N.pdf\", width=11.69, height=".($nGraf*3).")\npar(mfrow=c($nGraf,1))\n"; }
		#		else { print CMDR "pdf(\"".$outdir."/CNV_${$sampleNameRef}{$file}\_$N.pdf\", width=11.69, height=".($n*3).")\npar(mfrow=c($n,1))\n"; }
		#		}
		#	print CMDR "$cmdR";
		#	print CMDR "title(main=\"sample: ${$sampleNameRef}{$file}";
		#	if (${$PatientsRef}{$file}{"ecarte"}) { print CMDR " (invalid)\", col.main=\"red\""; }
		#	else { print CMDR "\""; }
		#	print CMDR ", outer=TRUE, line=-2, cex.main=2)\n";
		#	print CMDR "dev.off()\nquit(save=\"no\")\n";
		#	close CMDR;
		#	system "Rscript $outdir/${$sampleNameRef}{$file}\_temp.R";
		#	unlink "$outdir/${$sampleNameRef}{$file}\_temp.R";
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
#print CMDR "title(main=\"sample: ${$sampleNameRef}{$file}";
#if (${$PatientsRef}{$file}{"ecarte"}) { print CMDR " (invalid)\", col.main=\"red\""; }
#else { print CMDR "\""; }
#print CMDR ", outer=TRUE, line=-2, cex.main=2)\n";
print CMDR "dev.off()\nquit(save=\"no\")\n";
close CMDR;
system "Rscript $outdir/${$sampleNameRef}{$file}\_temp.R";
unlink "$outdir/${$sampleNameRef}{$file}\_temp.R";

}





####################
sub graphByChr2 {

my($outdir,$CNV_optRef,$patient,$PatientsRef,$RegionsRef,$ChrOrderRef,$regionOrderRef,$ResultsRef)= @_;

print "cmdR, for ${$PatientsRef}{$patient}{ID}\n";

my$norm = ${$CNV_optRef}{"norm"};
my$normGraf = $norm;
if ($norm eq "std") { $normGraf = "moy"; }

my$maxDepthGraph = ${$CNV_optRef}{"maxDepthGraph"};
my$seuil_deletion = ${$CNV_optRef}{"seuil_deletion"};
my$seuil_duplication = ${$CNV_optRef}{"seuil_duplication"};

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

	my$maxYsup=$seuil_duplication; my$maxYinf=$seuil_deletion;	
	foreach my$region (@{ $regionOrderRef->{$Chrom} }) {
		for (my$p=0;$p<scalar(keys%{$PatientsRef});$p++) {
			if (exists ${$RegionsRef}[$region]{$p}{"depth_ratio"}{$normGraf}) {
				if (${$RegionsRef}[$region]{$p}{"depth_ratio"}{$normGraf} > $maxYsup)
					{ $maxYsup = ${$RegionsRef}[$region]{$p}{"depth_ratio"}{$normGraf}; }
				if ($normGraf eq "std") {
					if (${$RegionsRef}[$region]{$p}{"depth_ratio"}{$normGraf} < $maxYinf)
						{ $maxYinf = ${$RegionsRef}[$region]{$p}{"depth_ratio"}{$normGraf}; }
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
		##all in 1 sheet:
#		$cmdR .= "par(fig=c(0,1,".(1-(($c+0.95)/$Nbr_Chr)).",".(1-(($c+0.05)/$Nbr_Chr))."), new=TRUE)
#plot (c(0,0), xlim=c(0,$maxX), ylim=c($maxYinf,$maxYsup), type =\"n\", main=\"chrom: $ChrName\", xlab=\"\", ylab=\"depth_ratio_to_$normGraf\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n";
		##1 sheet / chr
		$cmdR .= "plot (c(0,0), xlim=c(0,$maxX), ylim=c($maxYinf,$maxYsup), type =\"n\", main=\"chrom: $ChrName\", xlab=\"\", ylab=\"depth_ratio_to_$normGraf\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n";
		}
	else {
		##all in 1 sheet:
#		$cmdR .= "par(fig=c(0,1,".(1-(($c+0.95)/$Nbr_Chr)).",".(1-(($c+0.05)/$Nbr_Chr))."), new=TRUE)
#plot (c(0,0), xlim=c(0,$maxX), ylim=c(0,$maxYsup), type =\"n\", main=\"chrom: $ChrName\", xlab=\"\", ylab=\"depth_ratio_to_$normGraf\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n";
		$cmdR .= "plot (c(0,0), xlim=c(0,$maxX), ylim=c(0,$maxYsup), type =\"n\", main=\"chrom: $ChrName\", xlab=\"\", ylab=\"depth_ratio_to_$normGraf\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n";
		}

	#gene vertical separations
	my$currentGene=""; my$tmpTxt=""; my$Nbr_gene=0;
	for (my$r=0;$r<$Nbr_Reg;$r++) {
		if (${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"} ne "NA" && ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"geneID"} ne $currentGene) {
			$tmpTxt .= "abline(v=".($r+0.5).", col=\"blue\", lty = \"dotted\", lwd=1)\n";
			$currentGene = ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"geneID"};
			$Nbr_gene++;
			}
		}
	if ($Nbr_gene < ${$CNV_optRef}{"maxGeneSep"}) { $cmdR .= $tmpTxt; }

	my$Nbr_CNV=0;
	for (my$r=0;$r<$Nbr_Reg;$r++) {
		if (exists ${$ResultsRef}{$patient}{${$regionOrderRef}{$Chrom}[$r]})
			{ $Nbr_CNV++; }
		}

	#x labels
	my@printReg=();
	if ($maxX < ${$CNV_optRef}{"maxGeneLab"}) {
		##in grey if invalid
		for (my$r=0;$r<$Nbr_Reg;$r++) { 
			if (defined ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"Appel"}) { push(@printReg,$r) ; }
			}
		if (@printReg) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
			chop $cmdR;
			$cmdR .="), labels=c(";
			foreach my$r (@printReg) { $cmdR .= "\"".${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"label"}."\","; }
			chop $cmdR;
			$cmdR .= "), col.axis=\"darkgrey\", las=2";
			if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
			else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_Reg))." )\n"; }
			}
		##in black if valid
		@printReg=();
		for (my$r=0;$r<$Nbr_Reg;$r++) { 
			if(!defined ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"Appel"}) { push(@printReg,$r) ; }
			}
		if (@printReg) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
			chop $cmdR;
			$cmdR .="), labels=c(";
			foreach my$r (@printReg) { $cmdR .= "\"".${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"label"}."\","; }
			chop $cmdR;
			$cmdR .= "), col=\"black\", las=2";
			if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
			else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_Reg))." )\n"; }
			}
		}
	else {
		##in grey if invalid; only ticks
		for (my$r=0;$r<$Nbr_Reg;$r++) { 
			if (defined ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"Appel"}) { push(@printReg,$r) ; }
			}
		if (@printReg && scalar@printReg < ${$CNV_optRef}{"maxGeneSep"}) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
			chop $cmdR;
			$cmdR .= "), labels = FALSE, col.ticks=\"darkgrey\")\n";
			}
		##in black if valid
		@printReg=();
		if ($Nbr_gene < ${$CNV_optRef}{"maxGeneLab"}) {	# && ($Nbr_gene+$Nbr_CNV)>=${$CNV_optRef}{"maxGeneLab"}) {
			$currentGene="";
			for (my$r=0;$r<$Nbr_Reg;$r++) {
				if (${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"} ne "NA" && ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"geneID"} ne $currentGene) {
					push(@printReg,$r);
					$currentGene = ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"geneID"};
					}
				}
			}
		#elsif ($Nbr_gene<${$CNV_optRef}{"maxGeneLab"} && ($Nbr_gene+$Nbr_CNV)<${$CNV_optRef}{"maxGeneLab"}) {
		#	for (my$r=0;$r<$Nbr_Reg;$r++) {
		#		if ( (${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"Gene"} ne "NA" && ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"geneID"} ne $currentGene) || (exists ${$ResultsRef}{$patient}{${$regionOrderRef}{$Chrom}[$r]}) ) {
		#			push(@printReg,$r);
		#			$currentGene = ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"geneID"};
		#			}
		#		}
		#	}
		if (@printReg) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r+1).","; }
			chop $cmdR;
			$cmdR .="), labels=c(";
			foreach my$r (@printReg) { $cmdR .= "\"".${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{"label"}."\","; }
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
		if (@printReg && scalar@printReg < ${$CNV_optRef}{"maxGeneSep"}) {
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
					if (exists ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r2]]{$p}{"depth_ratio"}{$normGraf}) { $r2++;}
					else { last; }
					}
				if (($r2-1) > $r1) {
					$cmdR .= "lines( c(";
					for (my$r=$r1;$r<$r2;$r++)
						{ $cmdR .= ($r+1).","; }
					chop $cmdR;
					$cmdR .= "), c(";
					for (my$r=$r1;$r<$r2;$r++)
						{ $cmdR .= ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{$p}{"depth_ratio"}{$normGraf}.","; }
					chop $cmdR;
					$cmdR .= "), type =\"l\", lwd=1, col=\"darkgrey\")\n";
					}
				elsif (($r2-1) == $r1) {
					$cmdR .= "lines( c(".($r1+1)."), c(".${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r1]]{$p}{"depth_ratio"}{$normGraf}."), type =\"p\", lwd=1, col=\"darkgrey\")\n";
					}
				$r1 = ($r2+1);
				}
			}
		}

	##threshold lines (black):
	if ($norm eq "std" && $normGraf eq "moy") {
		foreach my$gender ("normByR_depth","normByR_depth_fem","normByR_depth_males") {
			my$r1=0;
			while ($r1<$Nbr_Reg) {
				my$r2=$r1;
				while ($r2<$Nbr_Reg) {
					if (${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r2]]{$gender}{"moy"}) { $r2++; }
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
							if ( ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{$gender}{"moy"} )
								{ $cmdR .= (1+$fold*(${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{$gender}{"std"} / ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{$gender}{"moy"})).","; }
							else 	{ $cmdR .= "1,"; }
							}
						chop $cmdR;
						$cmdR .= "), type =\"l\", lwd=1, col=\"black\")\n";
						}
					}
				elsif (($r2-1) == $r1) {
					$cmdR .= "lines( c(".($r1+1)."), c(".(1+(${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r1]]{$gender}{"std"} / ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r1]]{$gender}{"moy"}))."), type =\"p\", lwd=1, col=\"black\")\n";
					$cmdR .= "lines( c(".($r1+1)."), c(".(1-(${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r1]]{$gender}{"std"} / ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r1]]{$gender}{"moy"}))."), type =\"p\", lwd=1, col=\"black\")\n";
					}
				$r1 = ($r2+1);
				}
			}
		}
	else {
		$cmdR .= "abline(h=$seuil_deletion, col=\"black\", lty = \"dashed\", lwd=1)\n";
		$cmdR .= "abline(h=$seuil_duplication, col=\"black\", lty = \"dashed\", lwd=1)\n";
		}

	#target sample line (green):
	my$r1=0;
	while ($r1<$Nbr_Reg) {
		my$r2=$r1;
		while ($r2<$Nbr_Reg) {
			if (exists ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r2]]{$patient}{"depth_ratio"}{$normGraf}) { $r2++;}
			else { last; }
			}
		if (($r2-1) > $r1) {
			$cmdR .= "lines( c(";
			for (my$r=$r1;$r<$r2;$r++)
				{ $cmdR .= ($r+1).","; }
			chop $cmdR;
			$cmdR .= "), c(";
			for (my$r=$r1;$r<$r2;$r++)
				{ $cmdR .= ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{$patient}{"depth_ratio"}{$normGraf}.","; }
			chop $cmdR;
			$cmdR .= "), type =\"l\", lwd=1, col=\"green\")\n";
			}
		elsif (($r2-1) == $r1) {
			$cmdR .= "lines( c(".($r1+1)."), c(".${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r1]]{$patient}{"depth_ratio"}{$normGraf}."), type =\"p\", lwd=1, col=\"green\")\n";
			}
		$r1 = ($r2+1);
		}

	#points for CNVs
	my$points .= "points( c(";
	for (my$r=0;$r<$Nbr_Reg;$r++) {
		if (exists ${$ResultsRef}{$patient}{${$regionOrderRef}{$Chrom}[$r]} && exists ${$RegionsRef}[${$regionOrderRef}{$Chrom}[$r]]{$patient}{"depth_ratio"}{$normGraf})
			{ $points .= ($r+1).","; }
		}
	chop $points;
	$points .= "), c(";
	foreach my$region (@{ ${$regionOrderRef}{$Chrom} }) {
		if (exists ${$ResultsRef}{$patient}{$region} && exists ${$RegionsRef}[$region]{$patient}{"depth_ratio"}{$normGraf})
			{ $points .= ${$RegionsRef}[$region]{$patient}{"depth_ratio"}{$normGraf}.","; }
		}
	chop $points;
	$points .= "), type =\"p\", pch = 16, lwd=2, col=\"red\")\n";
	if ($Nbr_CNV) { $cmdR .= $points; }

	#if ($normGraf eq "std") { $cmdR .= "abline(h=0, col=\"darkgrey\", lty = \"dashed\", lwd=1)\n"; }
	#else { $cmdR .= "abline(h=1, col=\"darkgrey\", lty = \"dashed\", lwd=1)\n"; }
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

my($outdir,$CNV_optRef,$file,$FilesRef,$sampleNameRef,$RegionsRef,$ChromOrderRef,$regionOrderRef,$PatientsRef,$Result2Ref,$Result3Ref)= @_;

print "cmdR, for sample ${$sampleNameRef}{$file}\n";

my$norm = ${$CNV_optRef}{"norm"};
my$normGraf = $norm;
if ($norm eq "std") { $normGraf = "moy"; }

my$ext = ${$CNV_optRef}{"graphCNVpadding"};
my$maxDepthGraph = ${$CNV_optRef}{"maxDepthGraph"};
my$seuil_deletion = ${$CNV_optRef}{"seuil_deletion"};
my$seuil_duplication = ${$CNV_optRef}{"seuil_duplication"};

open (CMDR, ">$outdir/${$sampleNameRef}{$file}\_temp.R") || die;
print CMDR "#!/usr/bin/env Rscript\n\n" ;
print CMDR "pdf(\"$outdir/CNV_".${$sampleNameRef}{$file}.".pdf\", width=11.69, height=4,135)\n";

foreach my$CNV (sort{$a<=>$b}keys%{ ${$Result3Ref}{$file} }) {

	my$cmdR = "";

	my$Chrom = ${$RegionsRef}[$CNV]{"Chrom"};
	$Chrom =~ s/^chr//i;

	my$CNVend = ${$Result3Ref}{$file}{$CNV}{"end"};
	my$firstR = $CNV;
	for (my$r=($CNV - $ext);$r<$CNV;$r++) {
		if (exists ${$RegionsRef}[$r]) {
			$firstR = $r; last;
			}
		}
	my$lastR = $CNVend;
	for (my$r=$CNVend;$r<=($CNVend + $ext);$r++) {
		if (exists ${$RegionsRef}[$r]) {
			$lastR = $r;
			}
		else { last; }
		}
	my$Nbr_Reg = $lastR - $firstR + 1;

	my$maxYsup=$seuil_duplication; my$maxYinf=$seuil_deletion;	
	foreach my$region (@{ ${$regionOrderRef}{$Chrom} }) {
		foreach my$f (@{$FilesRef}) {
			if (exists ${$RegionsRef}[$region]{$f}{"depth_ratio"}{$normGraf}) {
				if (${$RegionsRef}[$region]{$f}{"depth_ratio"}{$normGraf} > $maxYsup)
					{ $maxYsup = ${$RegionsRef}[$region]{$f}{"depth_ratio"}{$normGraf}; }
				if ($normGraf eq "std") {
					if (${$RegionsRef}[$region]{$f}{"depth_ratio"}{$normGraf} < $maxYinf)
						{ $maxYinf = ${$RegionsRef}[$region]{$f}{"depth_ratio"}{$normGraf}; }
					}
				}
			}
		}
	if ($maxDepthGraph && $maxYsup > $maxDepthGraph) { $maxYsup = $maxDepthGraph; }
	if ($normGraf eq "std") {
		if ($maxDepthGraph && $maxYinf < (-$maxDepthGraph)) { $maxYinf = (-$maxDepthGraph); }
		}

	$cmdR .= "plot (c(0,0), xlim=c(0,$Nbr_Reg), ylim=c(";
	if ($normGraf eq "std") { $cmdR .= "$maxYinf,$maxYsup"; }
	else { $cmdR .= "0,$maxYsup"; }
	$cmdR .= "), type =\"n\", main=\"$Chrom:".${$RegionsRef}[$CNV]{"start"}."-".${$RegionsRef}[$CNVend]{"end"}."\", xlab=\"\", ylab=\"depth_ratio_to_$normGraf\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n";

	##gene separations
	my$currentGene=""; my$tmpTxt=""; my$Nbr_gene=0;
	for my$r ($firstR..$lastR) {
		if (${$RegionsRef}[$r]{"Gene"}) {
			if (${$RegionsRef}[$r]{"Gene"} ne "NA" && ${$RegionsRef}[$r]{"Gene"} ne $currentGene)  {
				$tmpTxt .= "abline(v=".($r-$firstR+0.5).", col=\"blue\", lty = \"dotted\", lwd=1)\n";
				$currentGene = ${$RegionsRef}[$r]{"Gene"};
				$Nbr_gene++;
				}
			}
		}
	if ($Nbr_gene < ${$CNV_optRef}{"maxGeneSep"}) { $cmdR .= $tmpTxt; }

	##region labels
	my@printReg=();
	if ($Nbr_Reg < ${$CNV_optRef}{"maxGeneLab"}) {
		##in grey if invalid
		for my$r ($firstR..$lastR){ 
			if (defined ${$RegionsRef}[$r]{"Appel"}) { push(@printReg,$r) ; }
			}
		if (@printReg) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r-$firstR+1).","; }
			chop $cmdR;
			$cmdR .="), labels=c(";
			foreach my$r (@printReg) { $cmdR .= "\"".${$RegionsRef}[$r]{"label"}."\","; }
			chop $cmdR;
			$cmdR .= "), col.axis=\"darkgrey\", las=2";
			if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
			else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_Reg))." )\n"; }
			}
		##in black if valid
		@printReg=();
		for my$r ($firstR..$lastR) { 
			if(!defined ${$RegionsRef}[$r]{"Appel"}) { push(@printReg,$r) ; }
			}
		if (@printReg) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r-$firstR+1).","; }
			chop $cmdR;
			$cmdR .="), labels=c(";
			foreach my$r (@printReg) { $cmdR .= "\"".${$RegionsRef}[$r]{"label"}."\","; }
			chop $cmdR;
			$cmdR .= "), col.axis=\"black\", las=2";
			if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
			else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_Reg))." )\n"; }
			}
		}
	else {
		##in grey if invalid; only ticks
		for my$r ($firstR..$lastR) { 
			if (defined ${$RegionsRef}[$r]{"Appel"}) { push(@printReg,$r) ; }
			}
		if (@printReg && scalar@printReg < ${$CNV_optRef}{"maxGeneSep"}) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r-$firstR+1).","; }
			chop $cmdR;
			$cmdR .= "), labels = FALSE, col.ticks=\"darkgrey\")\n";
			}
		##in black if valid
		@printReg=();
		if ($Nbr_gene < ${$CNV_optRef}{"maxGeneLab"}) {
			$currentGene="";
			for my$r ($firstR..$lastR) {
				if (${$RegionsRef}[$r]{"Gene"}) {
					if (${$RegionsRef}[$r]{"Gene"} ne "NA" && ${$RegionsRef}[$r]{"Gene"} ne $currentGene)  { 
						push(@printReg,$r); 
						$currentGene = ${$RegionsRef}[$r]{"Gene"};
						}
					}
				}
			}
		if (@printReg) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r-$firstR+1).","; }
			chop $cmdR;
			$cmdR .="), labels=c(";
			foreach my$r (@printReg) { $cmdR .= "\"".${$RegionsRef}[$r]{"Gene"}."\","; }
			chop $cmdR;
			$cmdR .= "), col.axis=\"black\", las=2";
			if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
			else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_gene))." )\n"; }
			}
		##in red if CNV; only ticks
		@printReg=();
		for (my$r=0;$r<$Nbr_Reg;$r++) {
			if (exists ${$Result2Ref}{$file}{$r}) { push(@printReg,$r) ; }
			}
		if (@printReg && scalar@printReg < ${$CNV_optRef}{"maxGeneSep"}) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r-$firstR+1).","; }
			chop $cmdR;
			$cmdR .= "), labels = FALSE, col.ticks=\"red\")\n";
			}


		}

	##all not target sample lines (grey):
	for (my$f=0;$f<scalar@{$FilesRef};$f++) {
		unless (${$FilesRef}[$f] eq $file) {
			my$r1 = $firstR;
			while ($r1 <= $lastR) {
				my$r2=$r1;
				while ($r2 <= $lastR) {
					if (exists ${$RegionsRef}[$r2]{${$FilesRef}[$f]}{"depth_ratio"}{$normGraf}) { $r2++; }
					else { last; }
					}
				if (($r2-1) > $r1) {
					$cmdR .= "lines( c(";
					for (my$r=$r1;$r<$r2;$r++)
						{ $cmdR .= ($r-$firstR+1).","; }
					chop $cmdR;
					$cmdR .= "), c(";
					for (my$r=$r1;$r<$r2;$r++)
						{ $cmdR .= ${$RegionsRef}[$r]{${$FilesRef}[$f]}{"depth_ratio"}{$normGraf}.","; }
					chop $cmdR;
					$cmdR .= "), type =\"l\", lwd=1, col=\"darkgrey\")\n";
					}
				elsif (($r2-1) == $r1) {
					$cmdR .= "lines( c(".($r1-$firstR+1)."), c(".${$RegionsRef}[$r1]{${$FilesRef}[$f]}{"depth_ratio"}{$normGraf}."), type =\"p\", lwd=1, col=\"darkgrey\")\n";
					}
				$r1 = ($r2+1);
				}
			}
		}

	##threshold lines (black):
	if ($norm eq "std" && $normGraf eq "moy") {
		foreach my$gender ("normByR_depth","normByR_depth_fem","normByR_depth_males") {
			my$r1 = $firstR;
			while ($r1 <= $lastR) {
				my$r2=$r1;
				while ($r2 <= $lastR) {
					if (${$RegionsRef}[$r2]{$gender}{"moy"}) { $r2++;}
					else { last; }
					}
				if (($r2-1) > $r1) {
					foreach my$fold ($seuil_duplication,$seuil_deletion) {
						$cmdR .= "lines( c(";
						for (my$r=$r1;$r<$r2;$r++)
							{ $cmdR .= ($r-$firstR+1).","; }
						chop $cmdR;
						$cmdR .= "), c(";
						for (my$r=$r1;$r<$r2;$r++) {
							if ( ${$RegionsRef}[$r]{$gender}{"moy"} ) {
								$cmdR .= (1+$fold*(${$RegionsRef}[$r]{$gender}{"std"} / ${$RegionsRef}[$r]{$gender}{"moy"})).",";
								}
							else 	{ $cmdR .= "1,"; }
							}
						chop $cmdR;
						$cmdR .= "), type =\"l\", lwd=1, col=\"black\")\n";
						}
					}
				elsif (($r2-1) == $r1) {
					$cmdR .= "lines( c(".($r1-$firstR+1)."), c(".(1+(${$RegionsRef}[$r1]{$gender}{"std"} / ${$RegionsRef}[$r1]{$gender}{"moy"}))."), type =\"p\", lwd=1, col=\"black\")\n";
					$cmdR .= "lines( c(".($r1-$firstR+1)."), c(".(1-(${$RegionsRef}[$r1]{$gender}{"std"} / ${$RegionsRef}[$r1]{$gender}{"moy"}))."), type =\"p\", lwd=1, col=\"black\")\n";
					}
				$r1 = ($r2+1);
				}
			}
		}
	else {
		$cmdR .= "abline(h=$seuil_deletion, col=\"black\", lty = \"dashed\", lwd=1)\n";
		$cmdR .= "abline(h=$seuil_duplication, col=\"black\", lty = \"dashed\", lwd=1)\n";
		}

	##target sample line (green)
	my$r1 = $firstR;
	while ($r1 <= $lastR) {
		my$r2=$r1;
		while ($r2 <= $lastR) {
			if (exists ${$RegionsRef}[$r2]{$file}{"depth_ratio"}{$normGraf}) { $r2++;}
			else { last; }
			}
		if (($r2-1) > $r1) {
			$cmdR .= "lines( c(";
			for (my$r=$r1;$r<$r2;$r++)
				{ $cmdR .= ($r-$firstR+1).","; }
			chop $cmdR;
			$cmdR .= "), c(";
			for (my$r=$r1;$r<$r2;$r++)
				{ $cmdR .= ${$RegionsRef}[$r]{$file}{"depth_ratio"}{$normGraf}.","; }
			chop $cmdR;
			$cmdR .= "), type =\"l\", lwd=1, col=\"green\")\n";
			}
		elsif (($r2-1) == $r1) {
			$cmdR .= "lines( c(".($r1-$firstR+1)."), c(".${$RegionsRef}[$r1]{$file}{"depth_ratio"}{$normGraf}."), type =\"p\", lwd=1, col=\"green\")\n";
			}
		$r1 = ($r2+1);
		}

	##points for CNVs
	my$points .= "points( c(";
	for my$r ($firstR..$lastR) {
		if (exists ${$Result2Ref}{$file}{$r} && exists ${$RegionsRef}[$r]{$file}{"depth_ratio"}{$normGraf})
			{ $points .= ($r-$firstR+1).","; }
		}
	chop $points;
	$points .= "), c(";
	for my$r ($firstR..$lastR){
		if (exists ${$Result2Ref}{$file}{$r} && exists ${$RegionsRef}[$r]{$file}{"depth_ratio"}{$normGraf})
			{ $points .= ${$RegionsRef}[$r]{$file}{"depth_ratio"}{$normGraf}.","; }
		}
	chop $points;
	$points .= "), type =\"p\", pch = 16, lwd=2, col=\"red\")\n";
	$cmdR .= $points;

	$cmdR .= "abline(h=1, col=\"darkgrey\", lty = \"dashed\", lwd=1)\n";

	print CMDR "$cmdR";
	
	}

##1 sheet / chr :
#print CMDR "title(main=\"sample: ${$sampleNameRef}{$file}";
#if (${$PatientsRef}{$file}{"ecarte"}) { print CMDR " (invalid)\", col.main=\"red\""; }
#else { print CMDR "\""; }
#print CMDR ", outer=TRUE, line=-2, cex.main=2)\n";
print CMDR "dev.off()\nquit(save=\"no\")\n";
close CMDR;
system "Rscript $outdir/${$sampleNameRef}{$file}\_temp.R";
unlink "$outdir/${$sampleNameRef}{$file}\_temp.R";

}


####################
sub graphByCNV2 {

my($outdir,$CNV_optRef,$patient,$PatientsRef,$RegionsRef,$ChrOrderRef,$regionOrderRef,$Result2Ref,$Result3Ref)= @_;

my$norm = ${$CNV_optRef}{"norm"};
my$normGraf = $norm;
if ($norm eq "std") { $normGraf = "moy"; }

my$ext = ${$CNV_optRef}{"graphCNVpadding"};
my$maxDepthGraph = ${$CNV_optRef}{"maxDepthGraph"};
my$seuil_deletion = ${$CNV_optRef}{"seuil_deletion"};
my$seuil_duplication = ${$CNV_optRef}{"seuil_duplication"};

print "cmdR, for ${$PatientsRef}{$patient}{ID}\n";
open (CMDR, ">$outdir/${$PatientsRef}{$patient}{ID}\_temp.R") || die;
print CMDR "#!/usr/bin/env Rscript\n\n" ;
print CMDR "pdf(\"".$outdir."/CNV_".${$PatientsRef}{$patient}{"ID"}.".pdf\", width=11.69, height=4.135)\n";

foreach my$CNV (sort{$a<=>$b}keys%{ ${$Result3Ref}{$patient} }) {

	my$cmdR = "";

	my$Chrom = ${$RegionsRef}[$CNV]{"Chrom"};
	my$ChrName = $Chrom; $ChrName =~ s/^chr//;

	my$CNVend = ${$Result3Ref}{$patient}{$CNV}{"end"};
	my$firstR = $CNV;
	for (my$r=($CNV - $ext);$r<$CNV;$r++) {
		if (exists ${$RegionsRef}[$r]) {
			$firstR = $r; last;
			}
		}
	my$lastR = $CNVend;
	for (my$r=$CNVend;$r<=($CNVend + $ext);$r++) {
		if (exists ${$RegionsRef}[$r]) {
			$lastR = $r;
			}
		else { last; }
		}
	my$Nbr_Reg = $lastR - $firstR + 1;

	my$maxYsup=$seuil_duplication; my$maxYinf=$seuil_deletion;	
	for my$r ($firstR..$lastR) {
		for (my$p=0;$p<scalar(keys%{$PatientsRef});$p++) {
			if (exists ${$RegionsRef}[$r]{$p}{"depth_ratio"}{$normGraf}) {
				if (${$RegionsRef}[$r]{$p}{"depth_ratio"}{$normGraf} > $maxYsup)
					{ $maxYsup = ${$RegionsRef}[$r]{$p}{"depth_ratio"}{$normGraf}; }
				if ($normGraf eq "std") {
					if (${$RegionsRef}[$r]{$p}{"depth_ratio"}{$normGraf} < $maxYinf)
						{ $maxYinf = ${$RegionsRef}[$r]{$p}{"depth_ratio"}{$normGraf}; }
					}
				}
			}
		}
	if ($maxDepthGraph && $maxYsup > $maxDepthGraph) { $maxYsup = $maxDepthGraph; }
	if ($normGraf eq "std") {
		if ($maxDepthGraph && $maxYinf < (-$maxDepthGraph)) { $maxYinf = (-$maxDepthGraph); }
		}

	if ($normGraf eq "std") {
		$cmdR .= "plot (c(0,0), xlim=c(0,$Nbr_Reg), ylim=c($maxYinf,$maxYsup), type =\"n\", main=\"$Chrom:".${$RegionsRef}[$CNV]{"start"}."-".${$RegionsRef}[$CNVend]{"end"}."\", xlab=\"\", ylab=\"depth_ratio_to_$normGraf\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n";
		}

	else {
		$cmdR .= "plot (c(0,0), xlim=c(0,$Nbr_Reg), ylim=c(0,$maxYsup), type =\"n\", main=\"$Chrom:".${$RegionsRef}[$CNV]{"start"}."-".${$RegionsRef}[$CNVend]{"end"}."\", xlab=\"\", ylab=\"depth_ratio_to_$normGraf\", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, xaxt=\"n\")\n";
		}

	#gene vertical separations
	my$currentGene=""; my$tmpTxt=""; my$Nbr_gene=0;
	for my$r ($firstR..$lastR) {
		if (${$RegionsRef}[$r]{"Gene"} ne "NA" && ${$RegionsRef}[$r]{"geneID"} ne $currentGene) {
			$tmpTxt .= "abline(v=".($r-$firstR+0.5).", col=\"blue\", lty = \"dotted\", lwd=1)\n";
			$currentGene = ${$RegionsRef}[$r]{"geneID"};
			$Nbr_gene++;
			}
		}
	if ($Nbr_gene < ${$CNV_optRef}{"maxGeneSep"}) { $cmdR .= $tmpTxt; }


	#x labels
	my@printReg=();
	if ($Nbr_Reg < ${$CNV_optRef}{"maxGeneLab"}) {
		##in grey if invalid
		for my$r ($firstR..$lastR) { 
			if (defined ${$RegionsRef}[$r]{"Appel"}) { push(@printReg,$r) ; }
			}
		if (@printReg) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r-$firstR+1).","; }
			chop $cmdR;
			$cmdR .="), labels=c(";
			foreach my$r (@printReg) { $cmdR .= "\"".${$RegionsRef}[$r]{"label"}."\","; }
			chop $cmdR;
			$cmdR .= "), col.axis=\"darkgrey\", las=2";
			if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
			else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_Reg))." )\n"; }
			}
		##in black if valid
		@printReg=();
		for my$r ($firstR..$lastR) { 
			if(!defined ${$RegionsRef}[$r]{"Appel"}) { push(@printReg,$r) ; }
			}
		if (@printReg) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r-$firstR+1).","; }
			chop $cmdR;
			$cmdR .="), labels=c(";
			foreach my$r (@printReg) { $cmdR .= "\"".${$RegionsRef}[$r]{"label"}."\","; }
			chop $cmdR;
			$cmdR .= "), col=\"black\", las=2";
			if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
			else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_Reg))." )\n"; }
			}
		}
	else {
		##in grey if invalid; only ticks
		for my$r ($firstR..$lastR) { 
			if (defined ${$RegionsRef}[$r]{"Appel"}) { push(@printReg,$r) ; }
			}
		if (@printReg && scalar@printReg < ${$CNV_optRef}{"maxGeneSep"}) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r-$firstR+1).","; }
			chop $cmdR;
			$cmdR .= "), labels = FALSE, col.ticks=\"darkgrey\")\n";
			}
		##in black if valid
		@printReg=();
		if ($Nbr_gene < ${$CNV_optRef}{"maxGeneLab"}) {	# && ($Nbr_gene+$Nbr_CNV)>=${$CNV_optRef}{"maxGeneLab"}) {
			$currentGene="";
			for my$r ($firstR..$lastR) {
				if (${$RegionsRef}[$r]{"Gene"} ne "NA" && ${$RegionsRef}[$r]{"geneID"} ne $currentGene) {
					push(@printReg,$r);
					$currentGene = ${$RegionsRef}[$r]{"geneID"};
					}
				}
			}
		if (@printReg) {
			$cmdR .= "axis(1, at=c(";
			foreach my$r (@printReg) { $cmdR .= ($r-$firstR+1).","; }
			chop $cmdR;
			$cmdR .="), labels=c(";
			foreach my$r (@printReg) { $cmdR .= "\"".${$RegionsRef}[$r]{"label"}."\","; }
			chop $cmdR;
			$cmdR .= "), col.axis=\"black\", las=2";
			if ($Nbr_gene<=5) { $cmdR .= ", cex.axis=1 )\n"; }
			else  { $cmdR .= ", cex.axis=".(log(5)/log($Nbr_gene))." )\n"; }
			}
		##in red if CNV; only ticks
		@printReg=();
		for my$r ($firstR..$lastR)  {
			if (exists ${$Result2Ref}{$patient}{$r}) { push(@printReg,$r) ; }
			}
		if (@printReg && scalar@printReg < ${$CNV_optRef}{"maxGeneSep"}) {
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
					if (exists ${$RegionsRef}[$r2]{$p}{"depth_ratio"}{$normGraf}) { $r2++;}
					else { last; }
					}
				if (($r2-1) > $r1) {
					$cmdR .= "lines( c(";
					for (my$r=$r1;$r<$r2;$r++)
						{ $cmdR .= ($r-$firstR+1).","; }
					chop $cmdR;
					$cmdR .= "), c(";
					for (my$r=$r1;$r<$r2;$r++)
						{ $cmdR .= ${$RegionsRef}[$r]{$p}{"depth_ratio"}{$normGraf}.","; }
					chop $cmdR;
					$cmdR .= "), type =\"l\", lwd=1, col=\"darkgrey\")\n";
					}
				elsif (($r2-1) == $r1) {
					$cmdR .= "lines( c(".($r1-$firstR+1)."), c(".${$RegionsRef}[$r1]{$p}{"depth_ratio"}{$normGraf}."), type =\"p\", lwd=1, col=\"darkgrey\")\n";
					}
				$r1 = ($r2+1);
				}
			}
		}

	##threshold lines (black):
	if ($norm eq "std" && $normGraf eq "moy") {
		foreach my$gender ("normByR_depth","normByR_depth_fem","normByR_depth_males") {
			my$r1 = $firstR;
			while ($r1 <= $lastR) {
				my$r2=$r1;
				while ($r2 <= $lastR) {
					if (${$RegionsRef}[$r2]{$gender}{"moy"}) { $r2++; }
					else { last; }
					}
				if (($r2-1) > $r1) {
					foreach my$fold ($seuil_duplication,$seuil_deletion) {
						$cmdR .= "lines( c(";
						for (my$r=$r1;$r<$r2;$r++)
							{ $cmdR .= ($r-$firstR+1).","; }
						chop $cmdR;
						$cmdR .= "), c(";
						for (my$r=$r1;$r<$r2;$r++) {
							if ( ${$RegionsRef}[$r]{$gender}{"moy"} )
								{ $cmdR .= (1+$fold*(${$RegionsRef}[$r]{$gender}{"std"} / ${$RegionsRef}[$r]{$gender}{"moy"})).","; }
							else 	{ $cmdR .= "1,"; }
							}
						chop $cmdR;
						$cmdR .= "), type =\"l\", lwd=1, col=\"black\")\n";
						}
					}
				elsif (($r2-1) == $r1) {
					$cmdR .= "lines( c(".($r1-$firstR+1)."), c(".(1+(${$RegionsRef}[$r1]{$gender}{"std"} / ${$RegionsRef}[$r1]{$gender}{"moy"}))."), type =\"p\", lwd=1, col=\"black\")\n";
					$cmdR .= "lines( c(".($r1-$firstR+1)."), c(".(1-(${$RegionsRef}[$r1]{$gender}{"std"} / ${$RegionsRef}[$r1]{$gender}{"moy"}))."), type =\"p\", lwd=1, col=\"black\")\n";
					}
				$r1 = ($r2+1);
				}
			}
		}
	else {
		$cmdR .= "abline(h=$seuil_deletion, col=\"black\", lty = \"dashed\", lwd=1)\n";
		$cmdR .= "abline(h=$seuil_duplication, col=\"black\", lty = \"dashed\", lwd=1)\n";
		}

	#target sample line (green):
	my$r1 = $firstR;
	while ($r1 <= $lastR) {
		my$r2=$r1;
		while ($r2 <= $lastR) {
			if (exists ${$RegionsRef}[$r2]{$patient}{"depth_ratio"}{$normGraf}) { $r2++;}
			else { last; }
			}
		if (($r2-1) > $r1) {
			$cmdR .= "lines( c(";
			for (my$r=$r1;$r<$r2;$r++)
				{ $cmdR .= ($r-$firstR+1).","; }
			chop $cmdR;
			$cmdR .= "), c(";
			for (my$r=$r1;$r<$r2;$r++)
				{ $cmdR .= ${$RegionsRef}[$r]{$patient}{"depth_ratio"}{$normGraf}.","; }
			chop $cmdR;
			$cmdR .= "), type =\"l\", lwd=1, col=\"green\")\n";
			}
		elsif (($r2-1) == $r1) {
			$cmdR .= "lines( c(".($r1-$firstR+1)."), c(".${$RegionsRef}[$r1]{$patient}{"depth_ratio"}{$normGraf}."), type =\"p\", lwd=1, col=\"green\")\n";
			}
		$r1 = ($r2+1);
		}

	#points for CNVs
	my$points .= "points( c(";
	for my$r ($firstR..$lastR) {
		if (exists ${$Result2Ref}{$patient}{$r} && exists ${$RegionsRef}[$r]{$patient}{"depth_ratio"}{$normGraf})
			{ $points .= ($r-$firstR+1).","; }
		}
	chop $points;
	$points .= "), c(";
	for my$r ($firstR..$lastR) {
		if (exists ${$Result2Ref}{$patient}{$r} && exists ${$RegionsRef}[$r]{$patient}{"depth_ratio"}{$normGraf})
			{ $points .= ${$RegionsRef}[$r]{$patient}{"depth_ratio"}{$normGraf}.","; }
		}
	chop $points;
	$points .= "), type =\"p\", pch = 16, lwd=2, col=\"red\")\n";
	$cmdR .= $points;

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


1;

