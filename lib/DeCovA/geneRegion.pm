package DeCovA::geneRegion;


use strict;
use warnings;



#########################
## NM2GeneRegion($Genes_r,$geneNM_r,$NMchr_r,$Regions{"raw"},$NM_Ex{"raw"});
# @{ $geneNM{$gene} } = [$NM1,...]
# %NMchr : key = NM, value = chr
# $Regions{chr}{NM}{start of region} = end of region
# %NM_Ex{NM}{start of region}{start of exon} = end of exon

sub NM2GeneRegion {

my($Genes,$geneNM,$NMchr,$RegionR,$NM_ExR) = @_;
my%RegionG;
foreach my$gene (keys%{$Genes}) {
	my@NMs =  @{ ${$geneNM}{$gene} };
	my$chr = ${$NMchr}{$NMs[0]};
#	print "change intervals for ".$gene." according to:\n\t".$NMs[0];
	%{ $RegionG{$chr}{$gene} } = %{ ${$RegionR}{$chr}{$NMs[0]} };
	my%geneNM_Ex;
	%{ $geneNM_Ex{$NMs[0]} } = %{ ${$NM_ExR}{$NMs[0]} };
	for(my$i=1;$i<scalar@NMs;$i++) {
#		print "\t".$NMs[$i];
		changeRegion1($NMs[$i],\%{ $RegionG{$chr}{$gene} },\%geneNM_Ex,\%{ ${$RegionR}{$chr}{$NMs[$i]} },\%{ ${$NM_ExR}{$NMs[$i]} });
		}
#	print "\n";
	foreach (@NMs) {
		%{ ${$NM_ExR}{$_} } = %{ $geneNM_Ex{$_} };
		}
	}
%{$RegionR} = %RegionG;
}


#########################
#make 1 common region from all transcripts of 1 gene
## changeRegion1($NMs[$i],$gene,$chr,\%{ $RegionG{$chr}{$gene} },\%geneNM_Ex,\%{ ${$RegionR}{$chr}{$NMs[$i]} },\%{ ${$NM_ExR}{$NMs[$i]} });

sub changeRegion1 {

my($NM,$RegionG1,$geneNM_Ex1,$RegionN,$NM_Ex_r) = @_;

my@StartG = sort{$a<=>$b}(keys%{$RegionG1});
my@StartN = sort{$a<=>$b}(keys%{$RegionN});
my(%RegionG2,%geneNM_Ex2);
my$c=0;	#$c: count of @StartG
my$i=0;	#$i: count of @StartN
my$endN = ${$RegionN}{$StartN[0]};

#if regionN ends before start of regionG
#while ( ($i < scalar@StartN) && ($endN < $StartG[0]) ) {
#	$RegionG2{$StartN[$i]} = $endN;
#	foreach my$startEx (keys%{ ${$NM_Ex_r}{$StartN[$i]} })
#		{ $geneNM_Ex2{$NM}{$StartN[$i]}{$startEx} = ${$NM_Ex_r}{$StartN[$i]}{$startEx}; }
#	$i++;
#	$endN = ${$RegionN}{$StartN[$i]};
#	}
while ( $i < scalar@StartN ) {
	$endN = ${$RegionN}{$StartN[$i]};				
	if ( $endN < $StartG[$c] ) {
		$RegionG2{$StartN[$i]} = $endN;
		foreach my$startEx (keys%{ ${$NM_Ex_r}{$StartN[$i]} })
			{ $geneNM_Ex2{$NM}{$StartN[$i]}{$startEx} = ${$NM_Ex_r}{$StartN[$i]}{$startEx}; }
		$i++; next;
		}
	while ( ($c < (scalar@StartG -1)) && ($StartN[$i] > ${$RegionG1}{$StartG[$c]}) ) {
		$RegionG2{$StartG[$c]} = ${$RegionG1}{$StartG[$c]};
		foreach my$nm (keys%{$geneNM_Ex1}) {
			foreach my$startEx (keys%{ ${$geneNM_Ex1}{$nm}{$StartG[$c]} })
				{ $geneNM_Ex2{$nm}{$StartG[$c]}{$startEx} = ${$geneNM_Ex1}{$nm}{$StartG[$c]}{$startEx}; }
			}
		$c++; 
		}
	while ( ($c < (scalar@StartG -1)) && ($StartN[$i] <= ${$RegionG1}{$StartG[$c]}) && ($endN >= $StartG[$c]) ) {
		if ($StartN[$i] < $StartG[$c]) {
			if ($endN > ${$RegionG1}{$StartG[$c]})
				{ $RegionG2{$StartN[$i]} = $endN; }
			else
				{ $RegionG2{$StartN[$i]} = ${$RegionG1}{$StartG[$c]}; }
			foreach my$nm (keys%{$geneNM_Ex1}) { 
				foreach my$startEx (keys%{ ${$geneNM_Ex1}{$nm}{$StartG[$c]} })
					{ $geneNM_Ex2{$nm}{$StartN[$i]}{$startEx} = ${$geneNM_Ex1}{$nm}{$StartG[$c]}{$startEx}; }
				}
			foreach my$startEx (keys%{ ${$NM_Ex_r}{$StartN[$i]} })
				{ $geneNM_Ex2{$NM}{$StartN[$i]}{$startEx} = ${$NM_Ex_r}{$StartN[$i]}{$startEx}; }
			}
		else	{
			if ($endN > ${$RegionG1}{$StartG[$c]})
				{ $RegionG2{$StartG[$c]} = $endN; }
			else
				{ $RegionG2{$StartG[$c]} = ${$RegionG1}{$StartG[$c]}; }
			foreach my$nm (keys%{$geneNM_Ex1}) { 
				foreach my$startEx (keys%{ ${$geneNM_Ex1}{$nm}{$StartG[$c]} })
					{ $geneNM_Ex2{$nm}{$StartG[$c]}{$startEx} = ${$geneNM_Ex1}{$nm}{$StartG[$c]}{$startEx}; }
				}
			foreach my$startEx (keys%{ ${$NM_Ex_r}{$StartN[$i]} })
				{ $geneNM_Ex2{$NM}{$StartG[$c]}{$startEx} = ${$NM_Ex_r}{$StartN[$i]}{$startEx}; }
			}
		$c++;
		}
	#current $c
	if ( ($endN < $StartG[$c]) && ($StartN[$i]>${$RegionG1}{$StartG[$c-1]}) ) {
		$RegionG2{$StartN[$i]} = $endN;			
		foreach my$startEx (keys%{ ${$NM_Ex_r}{$StartN[$i]} })
			{ $geneNM_Ex2{$NM}{$StartN[$i]}{$startEx} = ${$NM_Ex_r}{$StartN[$i]}{$startEx}; }
		}
	if ( $StartN[$i] > ${$RegionG1}{$StartG[$c]} ) {
		$RegionG2{$StartG[$c]} = ${$RegionG1}{$StartG[$c]};
		foreach my$nm (keys%{$geneNM_Ex1}) {
			foreach my$startEx (keys%{ ${$geneNM_Ex1}{$nm}{$StartG[$c]} })
				{ $geneNM_Ex2{$nm}{$StartG[$c]}{$startEx} = ${$geneNM_Ex1}{$nm}{$StartG[$c]}{$startEx}; }
			} 
		}
	if ( ($StartN[$i] <= ${$RegionG1}{$StartG[$c]}) && ($endN >= $StartG[$c]) ) {
		if ($StartN[$i] < $StartG[$c]) {
			if ($endN > ${$RegionG1}{$StartG[$c]})
				{ $RegionG2{$StartN[$i]} = $endN; }
			else
				{ $RegionG2{$StartN[$i]} = ${$RegionG1}{$StartG[$c]}; }
			foreach my$nm (keys%{$geneNM_Ex1}) { 
				foreach my$startEx (keys%{ ${$geneNM_Ex1}{$nm}{$StartG[$c]} })
					{ $geneNM_Ex2{$nm}{$StartN[$i]}{$startEx} = ${$geneNM_Ex1}{$nm}{$StartG[$c]}{$startEx}; }
				}
			foreach my$startEx (keys%{ ${$NM_Ex_r}{$StartN[$i]} })
				{ $geneNM_Ex2{$NM}{$StartN[$i]}{$startEx} = ${$NM_Ex_r}{$StartN[$i]}{$startEx}; }
			}
		else {
			if ($endN > ${$RegionG1}{$StartG[$c]})
				{ $RegionG2{$StartG[$c]} = $endN; }
			else
				{ $RegionG2{$StartG[$c]} = ${$RegionG1}{$StartG[$c]}; }
			foreach my$nm (keys%{$geneNM_Ex1}) { 
				foreach my$startEx (keys%{ ${$geneNM_Ex1}{$nm}{$StartG[$c]} })
					{ $geneNM_Ex2{$nm}{$StartG[$c]}{$startEx} = ${$geneNM_Ex1}{$nm}{$StartG[$c]}{$startEx}; }
				}
			foreach my$startEx (keys%{ ${$NM_Ex_r}{$StartN[$i]} })
				{ $geneNM_Ex2{$NM}{$StartG[$c]}{$startEx} = ${$NM_Ex_r}{$StartN[$i]}{$startEx}; }
			}
		}
	#if regionN keep going on after end of regionG
	if ( $StartN[$i] > ${$RegionG1}{$StartG[-1]} ) {
		$RegionG2{$StartN[$i]} = $endN;		
		foreach my$startEx (keys%{ ${$NM_Ex_r}{$StartN[$i]} })
			{ $geneNM_Ex2{$NM}{$StartN[$i]}{$startEx} = ${$NM_Ex_r}{$StartN[$i]}{$startEx}; }
		}
	$i++;
	}
#if regionG keep going on after end of regionN
while ( ($c < scalar@StartG) && ($StartG[$c] > $endN) ) {
	$RegionG2{$StartG[$c]} = ${$RegionG1}{$StartG[$c]};
	foreach my$nm (keys%{$geneNM_Ex1}) { 
		foreach my$startEx (keys%{ ${$geneNM_Ex1}{$nm}{$StartG[$c]} })
			{ $geneNM_Ex2{$nm}{$StartG[$c]}{$startEx} = ${$geneNM_Ex1}{$nm}{$StartG[$c]}{$startEx}; }
		}
	$c++; 
	}

#merge intervals
@StartG = sort{$a<=>$b}(keys%RegionG2);
my(%RegionG3,%geneNM_Ex3);
my$startReg = $StartG[0];
my$endReg = $RegionG2{$startReg};
foreach my$nm (keys%geneNM_Ex2) { 
	foreach my$startEx (keys%{ $geneNM_Ex2{$nm}{$startReg} })
		{ $geneNM_Ex3{$nm}{$startReg}{$startEx} = $geneNM_Ex2{$nm}{$startReg}{$startEx}; }
	}
for (my$i=1;$i<scalar(@StartG);$i++) {
	if ($StartG[$i] <= $endReg) { 
		if ($RegionG2{$StartG[$i]} > $endReg)
			{ $endReg = $RegionG2{$StartG[$i]}; }
		foreach my$nm (keys%geneNM_Ex2) {
			foreach my$startEx (keys%{ $geneNM_Ex2{$nm}{$StartG[$i]} })
				{ $geneNM_Ex3{$nm}{$startReg}{$startEx} = $geneNM_Ex2{$nm}{$StartG[$i]}{$startEx}; }
			}
		}
	else	{ 
		$RegionG3{$startReg} = $endReg;
		$startReg = $StartG[$i];
		$endReg = $RegionG2{$startReg};
		foreach my$nm (keys%geneNM_Ex2) {
			foreach my$startEx (keys%{ $geneNM_Ex2{$nm}{$startReg} })
				{ $geneNM_Ex3{$nm}{$startReg}{$startEx} = $geneNM_Ex2{$nm}{$startReg}{$startEx}; }
			}
		}
	}
$RegionG3{$startReg} = $endReg;

#return(\%RegionG3,\%geneNM_Ex3);
%{$RegionG1} = %RegionG3;
%{$geneNM_Ex1} = %geneNM_Ex3;

}


#########################
#extract bed from file, and intersect with %Regions (regions of selected genes) 
#create 1 hash :
#@hashSub = changeRegion2G(\@NMs,\%{ $Regions{"raw"}{$chr}{$gene} },\%{ $NM_Ex{"raw"} },\%{ $Bed{$chr} });
#%{ $RegBed{"raw"}{$gene} } = %{$hashSub[0]};	#$RegBed{$gene}{$startReg}{$startBed} = $endBed;

sub changeRegionG {

my($gene2NM,$Regions_r,$NM_Ex_r,$Bed_r) = @_;
# @{ ${$gene2NM}{$gene} } = [$NM1,...]
# ${$Regions_r}{chr}{NM}{start of region} = end of region
# ${$NM_Ex_r}{NM}{start of region}{start of exon} = end of exon
# ${$Bed_r}{$chr}{$start} = $end;
my%RegBed;	#$RegBed{$gene}{$startReg}{$startBed} = $endBed;

foreach my$chr (keys%{$Regions_r}) {

	my@bedStarts = sort{$a<=>$b}keys%{ ${$Bed_r}{$chr} };
	my$i1=0;		# idx within @bedStarts

	my%Genes = ();
	foreach my$gene (keys%{ ${$Regions_r}{$chr} }) {
		@{ $Genes{$gene}{"starts"} } = sort{$a<=>$b}(keys%{ ${$Regions_r}{$chr}{$gene} });
		}
	my@Genes = sort { $Genes{$a}{"starts"}[0]<=>$Genes{$b}{"starts"}[0] } (keys%Genes);	##sort by 1st start pos

	foreach my$gene (@Genes) {

		my(%interval2,%NM_Ex2,%RegBed1);

		my@regStarts = @{ $Genes{$gene}{"starts"} };
		my$c = 0;	#$c: idx within @regStarts

		while ( ($i1 < $#bedStarts) && (${$Bed_r}{$chr}{$bedStarts[$i1]} < $regStarts[0]) ) { $i1++; }	##reg progression to reach current bed interval

		my$i2 = $i1;
		my$endBed;
		while ( ($i2 < scalar@bedStarts) && ($bedStarts[$i2] <= ${$Regions_r}{$chr}{$gene}{$regStarts[-1]}) ) {
			$endBed = ${$Bed_r}{$chr}{$bedStarts[$i2]};
			if ( $endBed < $regStarts[$c] ) {
				if ($bedStarts[$i2] > ${$Regions_r}{$chr}{$gene}{$regStarts[$c-1]}) { 
					$interval2{$bedStarts[$i2]} = $endBed; 
					$RegBed1{$bedStarts[$i2]}{$bedStarts[$i2]} = $endBed;
					$i2++;
					next;
					}
				else { $c--; }
				}
			while ( ($c < $#regStarts) && ($bedStarts[$i2] > ${$Regions_r}{$chr}{$gene}{$regStarts[$c]}) ) {	##reg progression to reach current bed interval
				$interval2{$regStarts[$c]} = ${$Regions_r}{$chr}{$gene}{$regStarts[$c]};
				foreach my$NM (@{ ${$gene2NM}{$gene} }) {
					foreach my$startEx (keys%{ ${$NM_Ex_r}{$NM}{$regStarts[$c]} }) {
						$NM_Ex2{$NM}{$regStarts[$c]}{$startEx} = ${$NM_Ex_r}{$NM}{$regStarts[$c]}{$startEx};
						}
					}
				$c++; 
				}
			if ( $bedStarts[$i2] > ${$Regions_r}{$chr}{$gene}{$regStarts[$c]} ) { #for current $c
				$interval2{$regStarts[$c]} = ${$Regions_r}{$chr}{$gene}{$regStarts[$c]};
				foreach my$NM (@{ ${$gene2NM}{$gene} }) {
					foreach my$startEx (keys%{ ${$NM_Ex_r}{$NM}{$regStarts[$c]} }) {
						$NM_Ex2{$NM}{$regStarts[$c]}{$startEx} = ${$NM_Ex_r}{$NM}{$regStarts[$c]}{$startEx};
						}
					} 
				}
			while ( ($c < scalar@regStarts) && ($bedStarts[$i2] <= ${$Regions_r}{$chr}{$gene}{$regStarts[$c]}) && ($endBed >= $regStarts[$c]) ) {
				if ($bedStarts[$i2] < $regStarts[$c]) {
					if ($endBed > ${$Regions_r}{$chr}{$gene}{$regStarts[$c]}) {
						$interval2{$bedStarts[$i2]} = $endBed;
						}
					else {
						$interval2{$bedStarts[$i2]} = ${$Regions_r}{$chr}{$gene}{$regStarts[$c]};
						}
					$RegBed1{$bedStarts[$i2]}{$bedStarts[$i2]} = $endBed;
					foreach my$NM (@{ ${$gene2NM}{$gene} }) {
						foreach my$startEx (keys%{ ${$NM_Ex_r}{$NM}{$regStarts[$c]} }) {
							$NM_Ex2{$NM}{$bedStarts[$i2]}{$startEx} = ${$NM_Ex_r}{$NM}{$regStarts[$c]}{$startEx};
							}
						}
					}
				else {
					if ($endBed > ${$Regions_r}{$chr}{$gene}{$regStarts[$c]}) {
						$interval2{$regStarts[$c]} = $endBed;
						}
					else {
						$interval2{$regStarts[$c]} = ${$Regions_r}{$chr}{$gene}{$regStarts[$c]};
						}
					$RegBed1{$regStarts[$c]}{$bedStarts[$i2]} = $endBed;
					foreach my$NM (@{ ${$gene2NM}{$gene} }) {
						foreach my$startEx (keys%{ ${$NM_Ex_r}{$NM}{$regStarts[$c]} }) {
							$NM_Ex2{$NM}{$regStarts[$c]}{$startEx} = ${$NM_Ex_r}{$NM}{$regStarts[$c]}{$startEx};
							}
						}
					}
				$c++;
				}
			$i2++;
			}
		#if region keep going on after end of last bed interval
		while ( ($c < scalar@regStarts) && ($regStarts[$c] > $endBed) ) {
			$interval2{$regStarts[$c]} = ${$Regions_r}{$chr}{$gene}{$regStarts[$c]};
			foreach my$NM (@{ ${$gene2NM}{$gene} }) {
				foreach my$startEx (keys%{ ${$NM_Ex_r}{$NM}{$regStarts[$c]} }) {
					$NM_Ex2{$NM}{$regStarts[$c]}{$startEx} = ${$NM_Ex_r}{$NM}{$regStarts[$c]}{$startEx};
					}
				}
			$c++; 
			}

		#merge intervals
		my(%interval3,%NM_Ex3,%RegBed2);
		my@regStart2 = sort{$a<=>$b}(keys%interval2);
		my$startReg = $regStart2[0];
		my$endReg = $interval2{$startReg};
		foreach my$NM (keys%NM_Ex2) {
			foreach my$startEx (keys%{ $NM_Ex2{$NM}{$startReg} }) {
				$NM_Ex3{$NM}{$startReg}{$startEx} = $NM_Ex2{$NM}{$startReg}{$startEx};
				}
			}
		foreach my$startBed (keys%{ $RegBed1{$startReg} }) {
			$RegBed2{$startReg}{$startBed} = $RegBed1{$startReg}{$startBed};
			}
		for (my$i=1;$i<scalar(@regStart2);$i++) {
			if ($regStart2[$i] <= $endReg) { 
				if ($interval2{$regStart2[$i]} > $endReg) {
					$endReg = $interval2{$regStart2[$i]};
					}
				foreach my$NM (keys%NM_Ex2) {
					foreach my$startEx (keys%{ $NM_Ex2{$NM}{$regStart2[$i]} }) {
						$NM_Ex3{$NM}{$startReg}{$startEx} = $NM_Ex2{$NM}{$regStart2[$i]}{$startEx};
						}
					}
				foreach my$startBed (keys%{ $RegBed1{$regStart2[$i]} }) {
					$RegBed2{$startReg}{$startBed} = $RegBed1{$regStart2[$i]}{$startBed};
					}
				}
			else { 
				$interval3{$startReg} = $endReg;
				$startReg = $regStart2[$i];
				$endReg = $interval2{$startReg};
				foreach my$NM (keys%NM_Ex2) {
					foreach my$startEx (keys%{ $NM_Ex2{$NM}{$startReg} }) {
						$NM_Ex3{$NM}{$startReg}{$startEx} = $NM_Ex2{$NM}{$startReg}{$startEx};
						}
					}
				foreach my$startBed (keys%{ $RegBed1{$startReg} }) {
					$RegBed2{$startReg}{$startBed} = $RegBed1{$startReg}{$startBed};
					}
				}
			}
		$interval3{$startReg} = $endReg;

		%{ ${$Regions_r}{$chr}{$gene} }  = %interval3;
		foreach my$NM (keys%NM_Ex3) {
			%{ ${$NM_Ex_r}{$NM} } = %{ $NM_Ex3{$NM} };
			}
		%{ $RegBed{$gene} } = %RegBed2;

		}
	}

return(\%RegBed);
}

#########################
#extract bed from file, and intersect with %Regions (regions of selected genes) 
#create 1 hash :
#@hashSub = changeRegion2G(\@NMs,\%{ $Regions{"raw"}{$chr}{$gene} },\%{ $NM_Ex{"raw"} },\%{ $Bed{$chr} });
#%{ $RegBed{"raw"}{$gene} } = %{$hashSub[0]};	#$RegBed{$gene}{$startReg}{$startBed} = $endBed;

sub changeRegionG2 {

my($gene2NM,$Regions_r,$NM_Ex_r,$Bed_r) = @_;
# @{ ${$gene2NM}{$gene} } = [$NM1,...]
# ${$Regions_r}{chr}{NM}{start of region} = end of region
# ${$NM_Ex_r}{NM}{start of region}{start of exon} = end of exon
# ${$Bed_r}{$chr}{$start} = $end;
my%RegBed;	#$RegBed{$gene}{$startReg}{$startBed} = $endBed;

foreach my$chr (keys%{$Regions_r}) {

	my@bedStarts = sort{$a<=>$b}keys%{ ${$Bed_r}{$chr} };
	my$i1=0;		# idx within @bedStarts

	my%Genes = ();
	foreach my$gene (keys%{ ${$Regions_r}{$chr} }) {
		@{ $Genes{$gene}{"starts"} } = sort{$a<=>$b}(keys%{ ${$Regions_r}{$chr}{$gene} });
		}
	my@Genes = sort { $Genes{$a}{"starts"}[0]<=>$Genes{$b}{"starts"}[0] } (keys%Genes);	##sort by 1st start pos

	foreach my$gene (@Genes) {

		my(%interval2,%NM_Ex2,%RegBed1);

		my@regStarts = @{ $Genes{$gene}{"starts"} };
		my$c1 = 0;	#$c: idx within @regStarts

		I1LOOP: while ( ($i1 < $#bedStarts) && (${$Bed_r}{$chr}{$bedStarts[$i1]} < $regStarts[0]) ) { $i1++; }			##bed progression to reach 1st reg interval

		my$i2 = $i1;
		my$endBed;
		I2LOOP: while ( ($i2 < scalar@bedStarts) && ($bedStarts[$i2] <= ${$Regions_r}{$chr}{$gene}{$regStarts[-1]}) ) {

			$endBed = ${$Bed_r}{$chr}{$bedStarts[$i2]};

			C1LOOP: while ( ($c1 < $#regStarts) && ($bedStarts[$i2] > ${$Regions_r}{$chr}{$gene}{$regStarts[$c1]}) ) {	##reg progression to reach current bed interval
				if ( $endBed < $regStarts[$c1+1] ) {			## bed interval between 2 reg intervals
					$interval2{$bedStarts[$i2]} = $endBed; 
					$RegBed1{$bedStarts[$i2]}{$bedStarts[$i2]} = $endBed;
					$i2++;
					next I2LOOP;
					}
				else {
					unless (exists $interval2{$regStarts[$c1]}) {
						$interval2{$regStarts[$c1]} = ${$Regions_r}{$chr}{$gene}{$regStarts[$c1]};
						foreach my$NM (@{ ${$gene2NM}{$gene} }) {
							foreach my$startEx (keys%{ ${$NM_Ex_r}{$NM}{$regStarts[$c1]} }) {
								$NM_Ex2{$NM}{$regStarts[$c1]}{$startEx} = ${$NM_Ex_r}{$NM}{$regStarts[$c1]}{$startEx};
								}
							}
						}
					}
				$c1++; 
				}
			if ( $bedStarts[$i2] > ${$Regions_r}{$chr}{$gene}{$regStarts[$c1]} ) { #for last $c
				unless (exists $interval2{$regStarts[$c1]}) {
					$interval2{$regStarts[$c1]} = ${$Regions_r}{$chr}{$gene}{$regStarts[$c1]};
					foreach my$NM (@{ ${$gene2NM}{$gene} }) {
						foreach my$startEx (keys%{ ${$NM_Ex_r}{$NM}{$regStarts[$c1]} }) {
							$NM_Ex2{$NM}{$regStarts[$c1]}{$startEx} = ${$NM_Ex_r}{$NM}{$regStarts[$c1]}{$startEx};
							}
						}
					}
				}

			my$c2 = $c1;
			C2LOOP: while ( ($c2 < scalar@regStarts) && ($bedStarts[$i2] <= ${$Regions_r}{$chr}{$gene}{$regStarts[$c2]}) && ($endBed >= $regStarts[$c2]) ) {
				if ($bedStarts[$i2] < $regStarts[$c2]) {
					if ($endBed > ${$Regions_r}{$chr}{$gene}{$regStarts[$c2]}) {
						$interval2{$bedStarts[$i2]} = $endBed;
						}
					else {
						$interval2{$bedStarts[$i2]} = ${$Regions_r}{$chr}{$gene}{$regStarts[$c2]};
						}
					$RegBed1{$bedStarts[$i2]}{$bedStarts[$i2]} = $endBed;
					foreach my$NM (@{ ${$gene2NM}{$gene} }) {
						foreach my$startEx (keys%{ ${$NM_Ex_r}{$NM}{$regStarts[$c2]} }) {
							$NM_Ex2{$NM}{$bedStarts[$i2]}{$startEx} = ${$NM_Ex_r}{$NM}{$regStarts[$c2]}{$startEx};
							}
						}
					}
				else {
					if ($endBed > ${$Regions_r}{$chr}{$gene}{$regStarts[$c2]}) {
						$interval2{$regStarts[$c2]} = $endBed;
						}
					else {
						$interval2{$regStarts[$c2]} = ${$Regions_r}{$chr}{$gene}{$regStarts[$c2]};
						}
					$RegBed1{$regStarts[$c2]}{$bedStarts[$i2]} = $endBed;
					foreach my$NM (@{ ${$gene2NM}{$gene} }) {
						foreach my$startEx (keys%{ ${$NM_Ex_r}{$NM}{$regStarts[$c2]} }) {
							$NM_Ex2{$NM}{$regStarts[$c2]}{$startEx} = ${$NM_Ex_r}{$NM}{$regStarts[$c2]}{$startEx};
							}
						}
					}
				$c2++;
				}
			$i2++;
			}

		#if region keep going on after end of last bed interval
		while ($c1 < scalar@regStarts) {
			if ($regStarts[$c1] > $endBed) {
				$interval2{$regStarts[$c1]} = ${$Regions_r}{$chr}{$gene}{$regStarts[$c1]};
				foreach my$NM (@{ ${$gene2NM}{$gene} }) {
					foreach my$startEx (keys%{ ${$NM_Ex_r}{$NM}{$regStarts[$c1]} }) {
						$NM_Ex2{$NM}{$regStarts[$c1]}{$startEx} = ${$NM_Ex_r}{$NM}{$regStarts[$c1]}{$startEx};
						}
					}
				}
			$c1++; 
			}

		#merge intervals
		my(%interval3,%NM_Ex3,%RegBed2);
		my@regStart2 = sort{$a<=>$b}(keys%interval2);
		my$startReg = $regStart2[0];
		my$endReg = $interval2{$startReg};
		foreach my$NM (keys%NM_Ex2) {
			foreach my$startEx (keys%{ $NM_Ex2{$NM}{$startReg} }) {
				$NM_Ex3{$NM}{$startReg}{$startEx} = $NM_Ex2{$NM}{$startReg}{$startEx};
				}
			}
		foreach my$startBed (keys%{ $RegBed1{$startReg} }) {
			$RegBed2{$startReg}{$startBed} = $RegBed1{$startReg}{$startBed};
			}
		for (my$i=1;$i<scalar(@regStart2);$i++) {
			if ($regStart2[$i] <= $endReg) { 
				if ($interval2{$regStart2[$i]} > $endReg) {
					$endReg = $interval2{$regStart2[$i]};
					}
				foreach my$NM (keys%NM_Ex2) {
					foreach my$startEx (keys%{ $NM_Ex2{$NM}{$regStart2[$i]} }) {
						$NM_Ex3{$NM}{$startReg}{$startEx} = $NM_Ex2{$NM}{$regStart2[$i]}{$startEx};
						}
					}
				foreach my$startBed (keys%{ $RegBed1{$regStart2[$i]} }) {
					$RegBed2{$startReg}{$startBed} = $RegBed1{$regStart2[$i]}{$startBed};
					}
				}
			else { 
				$interval3{$startReg} = $endReg;
				$startReg = $regStart2[$i];
				$endReg = $interval2{$startReg};
				foreach my$NM (keys%NM_Ex2) {
					foreach my$startEx (keys%{ $NM_Ex2{$NM}{$startReg} }) {
						$NM_Ex3{$NM}{$startReg}{$startEx} = $NM_Ex2{$NM}{$startReg}{$startEx};
						}
					}
				foreach my$startBed (keys%{ $RegBed1{$startReg} }) {
					$RegBed2{$startReg}{$startBed} = $RegBed1{$startReg}{$startBed};
					}
				}
			}
		$interval3{$startReg} = $endReg;

		%{ ${$Regions_r}{$chr}{$gene} }  = %interval3;
		foreach my$NM (keys%NM_Ex3) {
			%{ ${$NM_Ex_r}{$NM} } = %{ $NM_Ex3{$NM} };
			}
		%{ $RegBed{$gene} } = %RegBed2;

		}
	}

return(\%RegBed);
}


#########################
#extract bed from file, and intersect with %Regions (regions of selected genes) 
#create 1 hash :
#@hashSub = changeRegion2N(\%{ $Regions{"raw"}{$NMchr{$NM}}{$NM} },\%{ $NM_Ex{"raw"}{$NM} },\%{ $Bed{$NMchr{$NM}} });
#%{ $RegBed{"raw"}{$NM} } = %{$hashSub[0]};	#$RegBed{$NM}{$startReg}{$startBed} = $endBed;

sub changeRegionN2 {

my($Regions_r,$NM_Ex_r,$Bed_r) = @_;
# ${$Regions_r}{chr}{NM}{start of region} = end of region
# ${$NM_Ex_r}{NM}{start of region}{start of exon} = end of exon
# ${$Bed_r}{$chr}{$start} = $end;
my%RegBed;	#$RegBed{$NM}{$startReg}{$startBed} = $endBed;

foreach my$chr (keys%{$Regions_r}) {

	my@bedStarts = sort{$a<=>$b}keys%{ ${$Bed_r}{$chr} };
	my$i1=0;		# idx within @bedStarts

	my%NMs = ();
	foreach my$nm (keys%{ ${$Regions_r}{$chr} }) {
		@{ $NMs{$nm}{"starts"} } = sort{$a<=>$b}(keys%{ ${$Regions_r}{$chr}{$nm} });
		}
	my@NMs = sort { $NMs{$a}{"starts"}[0]<=>$NMs{$b}{"starts"}[0] } (keys%NMs);	##sort by 1st start pos

	foreach my$nm (@NMs) {

		my(%interval2,%NM_Ex2,%RegBed1);

		my@regStarts = @{ $NMs{$nm}{"starts"} };
		my$c1 = 0;	#$c: idx within @regStarts

		I1LOOP: while ( ($i1 < $#bedStarts) && (${$Bed_r}{$chr}{$bedStarts[$i1]} < $regStarts[0]) ) { $i1++; }			##bed progression to reach 1st reg interval

		my$i2 = $i1;
		my$endBed;

		I2LOOP: while ( ($i2 < scalar@bedStarts) && ($bedStarts[$i2] <= ${$Regions_r}{$chr}{$nm}{$regStarts[-1]}) ) {

			$endBed = ${$Bed_r}{$chr}{$bedStarts[$i2]};

			C1LOOP: while ( ($c1 < $#regStarts) && ($bedStarts[$i2] > ${$Regions_r}{$chr}{$nm}{$regStarts[$c1]}) ) {	##reg progression to reach current bed interval
				if ( $endBed < $regStarts[$c1+1] ) {			## bed interval between 2 reg intervals
					$interval2{$bedStarts[$i2]} = $endBed; 
					$RegBed1{$bedStarts[$i2]}{$bedStarts[$i2]} = $endBed;
					$i2++;
					next I2LOOP;
					}
				else {
					unless (exists $interval2{$regStarts[$c1]}) {
						$interval2{$regStarts[$c1]} = ${$Regions_r}{$chr}{$nm}{$regStarts[$c1]};
						foreach my$startEx (keys%{ ${$NM_Ex_r}{$nm}{$regStarts[$c1]} }) {
							$NM_Ex2{$regStarts[$c1]}{$startEx} = ${$NM_Ex_r}{$nm}{$regStarts[$c1]}{$startEx};
							}
						}
					}
				$c1++; 
				}
			if ( $bedStarts[$i2] > ${$Regions_r}{$chr}{$nm}{$regStarts[$c1]} ) { #for last $c
				unless (exists $interval2{$regStarts[$c1]}) {
					$interval2{$regStarts[$c1]} = ${$Regions_r}{$chr}{$nm}{$regStarts[$c1]};
					foreach my$startEx (keys%{ ${$NM_Ex_r}{$nm}{$regStarts[$c1]} }) {
						$NM_Ex2{$regStarts[$c1]}{$startEx} = ${$NM_Ex_r}{$nm}{$regStarts[$c1]}{$startEx};
						}
					}
				}

			my$c2 = $c1;
			C2LOOP: while ( ($c2 < scalar@regStarts) && ($bedStarts[$i2] <= ${$Regions_r}{$chr}{$nm}{$regStarts[$c2]}) && ($endBed >= $regStarts[$c2]) ) {
				if ($bedStarts[$i2] < $regStarts[$c2]) {
					if ($endBed > ${$Regions_r}{$chr}{$nm}{$regStarts[$c2]}) {
						$interval2{$bedStarts[$i2]} = $endBed;
						}
					else {
						$interval2{$bedStarts[$i2]} = ${$Regions_r}{$chr}{$nm}{$regStarts[$c2]};
						}
					$RegBed1{$bedStarts[$i2]}{$bedStarts[$i2]} = $endBed;
					foreach my$startEx (keys%{ ${$NM_Ex_r}{$nm}{$regStarts[$c2]} }) {
						$NM_Ex2{$bedStarts[$i2]}{$startEx} = ${$NM_Ex_r}{$nm}{$regStarts[$c2]}{$startEx};
						}
					}
				else {
					if ($endBed > ${$Regions_r}{$chr}{$nm}{$regStarts[$c2]}) {
						$interval2{$regStarts[$c2]} = $endBed;
						}
					else {
						$interval2{$regStarts[$c2]} = ${$Regions_r}{$chr}{$nm}{$regStarts[$c2]};
						}
					$RegBed1{$regStarts[$c2]}{$bedStarts[$i2]} = $endBed;
					foreach my$startEx (keys%{ ${$NM_Ex_r}{$nm}{$regStarts[$c2]} }) {
						$NM_Ex2{$regStarts[$c2]}{$startEx} = ${$NM_Ex_r}{$nm}{$regStarts[$c2]}{$startEx};
						}
					}
				$c2++;
				}
			$i2++;
			}
		#if region keep going on after end of bed
		while ($c1 < scalar@regStarts) {
			if ($regStarts[$c1] > $endBed) {
				$interval2{$regStarts[$c1]} = ${$Regions_r}{$chr}{$nm}{$regStarts[$c1]};
				foreach my$startEx (keys%{ ${$NM_Ex_r}{$nm}{$regStarts[$c1]} }) {
					$NM_Ex2{$regStarts[$c1]}{$startEx} = ${$NM_Ex_r}{$nm}{$regStarts[$c1]}{$startEx};
					}
				}
			$c1++; 
			}

		#merge intervals
		my(%interval3,%NM_Ex3,%RegBed2);
		my@regStart2 = sort{$a<=>$b}(keys%interval2);
		my$startReg = $regStart2[0];
		my$endReg = $interval2{$startReg};
		foreach my$startEx (keys%{ $NM_Ex2{$startReg} }) {
			$NM_Ex3{$startReg}{$startEx} = $NM_Ex2{$startReg}{$startEx};
			}
		foreach my$startBed (keys%{ $RegBed1{$startReg} }) {
			$RegBed2{$startReg}{$startBed} = $RegBed1{$startReg}{$startBed};
			}
		for (my$i=1;$i<scalar(@regStart2);$i++) {
			if ($regStart2[$i] <= $endReg) { 
				if ($interval2{$regStart2[$i]} > $endReg) {
					$endReg = $interval2{$regStart2[$i]};
					}
				foreach my$startEx (keys%{ $NM_Ex2{$regStart2[$i]} }) {
					$NM_Ex3{$startReg}{$startEx} = $NM_Ex2{$regStart2[$i]}{$startEx};
					}
				foreach my$startBed (keys%{ $RegBed1{$regStart2[$i]} }) {
					$RegBed2{$startReg}{$startBed} = $RegBed1{$regStart2[$i]}{$startBed};
					}
				}
			else { 
				$interval3{$startReg} = $endReg;
				$startReg = $regStart2[$i];
				$endReg = $interval2{$startReg};
				foreach my$startEx (keys%{ $NM_Ex2{$startReg} }) {
					$NM_Ex3{$startReg}{$startEx} = $NM_Ex2{$startReg}{$startEx};
					}
				foreach my$startBed (keys%{ $RegBed1{$startReg} }) {
					$RegBed2{$startReg}{$startBed} = $RegBed1{$startReg}{$startBed};
					}
				}
			}
		$interval3{$startReg} = $endReg;

		%{ ${$Regions_r}{$chr}{$nm} }  = %interval3;
		%{ ${$NM_Ex_r}{$nm} } = %NM_Ex3;
		%{ $RegBed{$nm} } = %RegBed2;

		}
	}

return(\%RegBed);
}

#########################
#extract bed from file, and intersect with %Regions (regions of selected genes) 
#create 1 hash :
#@hashSub = changeRegion2N(\%{ $Regions{"raw"}{$NMchr{$NM}}{$NM} },\%{ $NM_Ex{"raw"}{$NM} },\%{ $Bed{$NMchr{$NM}} });
#%{ $RegBed{"raw"}{$NM} } = %{$hashSub[0]};	#$RegBed{$NM}{$startReg}{$startBed} = $endBed;

sub changeRegionN {

my($Regions_r,$NM_Ex_r,$Bed_r) = @_;
# ${$Regions_r}{chr}{NM}{start of region} = end of region
# ${$NM_Ex_r}{NM}{start of region}{start of exon} = end of exon
# ${$Bed_r}{$chr}{$start} = $end;
my%RegBed;	#$RegBed{$NM}{$startReg}{$startBed} = $endBed;

foreach my$chr (keys%{$Regions_r}) {

	my@bedStarts = sort{$a<=>$b}keys%{ ${$Bed_r}{$chr} };
	my$i1=0;		# idx within @bedStarts

	my%NMs = ();
	foreach my$nm (keys%{ ${$Regions_r}{$chr} }) {
		@{ $NMs{$nm}{"starts"} } = sort{$a<=>$b}(keys%{ ${$Regions_r}{$chr}{$nm} });
		}
	my@NMs = sort { $NMs{$a}{"starts"}[0]<=>$NMs{$b}{"starts"}[0] } (keys%NMs);	##sort by 1st start pos

	foreach my$nm (@NMs) {

		my(%interval2,%NM_Ex2,%RegBed1);

		my@regStarts = @{ $NMs{$nm}{"starts"} };
		my$c = 0;	#$c: idx within @regStarts

		while ( ($i1 < $#bedStarts) && (${$Bed_r}{$chr}{$bedStarts[$i1]} < $regStarts[0]) ) { $i1++; }

		my$i2 = $i1;
		my$endBed;

		while ( ($i2 < scalar@bedStarts) && ($bedStarts[$i2] <= ${$Regions_r}{$chr}{$nm}{$regStarts[-1]}) ) {
			$endBed = ${$Bed_r}{$chr}{$bedStarts[$i2]};
			if ( $endBed < $regStarts[$c] ) {
				if ($bedStarts[$i2] > ${$Regions_r}{$chr}{$nm}{$regStarts[$c-1]}) { 
					$interval2{$bedStarts[$i2]} = $endBed; 
					$RegBed1{$bedStarts[$i2]}{$bedStarts[$i2]} = $endBed;
					$i2++;
					next;
					}
				else { $c--; }
				}
			while ( ($c < $#regStarts) && ($bedStarts[$i2] > ${$Regions_r}{$chr}{$nm}{$regStarts[$c]}) ) { 
				$interval2{$regStarts[$c]} = ${$Regions_r}{$chr}{$nm}{$regStarts[$c]};
				foreach my$startEx (keys%{ ${$NM_Ex_r}{$nm}{$regStarts[$c]} }) {
					$NM_Ex2{$regStarts[$c]}{$startEx} = ${$NM_Ex_r}{$nm}{$regStarts[$c]}{$startEx};
					}
				$c++; 
				}
			if ( $bedStarts[$i2] > ${$Regions_r}{$chr}{$nm}{$regStarts[$c]} ) {	#for current $c 
				$interval2{$regStarts[$c]} = ${$Regions_r}{$chr}{$nm}{$regStarts[$c]};
				foreach my$startEx (keys%{ ${$NM_Ex_r}{$nm}{$regStarts[$c]} }) {
					$NM_Ex2{$regStarts[$c]}{$startEx} = ${$NM_Ex_r}{$nm}{$regStarts[$c]}{$startEx};
					}
				}

			while ( ($c < scalar@regStarts) && ($bedStarts[$i2] <= ${$Regions_r}{$chr}{$nm}{$regStarts[$c]}) && ($endBed >= $regStarts[$c]) ) {
				if ($bedStarts[$i2] < $regStarts[$c]) {
					if ($endBed > ${$Regions_r}{$chr}{$nm}{$regStarts[$c]}) {
						$interval2{$bedStarts[$i2]} = $endBed;
						}
					else {
						$interval2{$bedStarts[$i2]} = ${$Regions_r}{$chr}{$nm}{$regStarts[$c]};
						}
					$RegBed1{$bedStarts[$i2]}{$bedStarts[$i2]} = $endBed;
					foreach my$startEx (keys%{ ${$NM_Ex_r}{$nm}{$regStarts[$c]} }) {
						$NM_Ex2{$bedStarts[$i2]}{$startEx} = ${$NM_Ex_r}{$nm}{$regStarts[$c]}{$startEx};
						}
					}
				else {
					if ($endBed > ${$Regions_r}{$chr}{$nm}{$regStarts[$c]}) {
						$interval2{$regStarts[$c]} = $endBed;
						}
					else {
						$interval2{$regStarts[$c]} = ${$Regions_r}{$chr}{$nm}{$regStarts[$c]};
						}
					$RegBed1{$regStarts[$c]}{$bedStarts[$i2]} = $endBed;
					foreach my$startEx (keys%{ ${$NM_Ex_r}{$nm}{$regStarts[$c]} }) {
						$NM_Ex2{$regStarts[$c]}{$startEx} = ${$NM_Ex_r}{$nm}{$regStarts[$c]}{$startEx};
						}
					}
				$c++;
				}
			}
		#if region keep going on after end of bed
		while ( ($c < scalar@regStarts) && ($regStarts[$c] > $endBed) ) {
			$interval2{$regStarts[$c]} = ${$Regions_r}{$chr}{$nm}{$regStarts[$c]};
			foreach my$startEx (keys%{ ${$NM_Ex_r}{$nm}{$regStarts[$c]} }) {
				$NM_Ex2{$regStarts[$c]}{$startEx} = ${$NM_Ex_r}{$nm}{$regStarts[$c]}{$startEx};
				}
			$c++; 
			}

		#merge intervals
		my(%interval3,%NM_Ex3,%RegBed2);
		my@regStart2 = sort{$a<=>$b}(keys%interval2);
		my$startReg = $regStart2[0];
		my$endReg = $interval2{$startReg};
		foreach my$startEx (keys%{ $NM_Ex2{$startReg} }) {
			$NM_Ex3{$startReg}{$startEx} = $NM_Ex2{$startReg}{$startEx};
			}
		foreach my$startBed (keys%{ $RegBed1{$startReg} }) {
			$RegBed2{$startReg}{$startBed} = $RegBed1{$startReg}{$startBed};
			}
		for (my$i=1;$i<scalar(@regStart2);$i++) {
			if ($regStart2[$i] <= $endReg) { 
				if ($interval2{$regStart2[$i]} > $endReg) {
					$endReg = $interval2{$regStart2[$i]};
					}
				foreach my$startEx (keys%{ $NM_Ex2{$regStart2[$i]} }) {
					$NM_Ex3{$startReg}{$startEx} = $NM_Ex2{$regStart2[$i]}{$startEx};
					}
				foreach my$startBed (keys%{ $RegBed1{$regStart2[$i]} }) {
					$RegBed2{$startReg}{$startBed} = $RegBed1{$regStart2[$i]}{$startBed};
					}
				}
			else { 
				$interval3{$startReg} = $endReg;
				$startReg = $regStart2[$i];
				$endReg = $interval2{$startReg};
				foreach my$startEx (keys%{ $NM_Ex2{$startReg} }) {
					$NM_Ex3{$startReg}{$startEx} = $NM_Ex2{$startReg}{$startEx};
					}
				foreach my$startBed (keys%{ $RegBed1{$startReg} }) {
					$RegBed2{$startReg}{$startBed} = $RegBed1{$startReg}{$startBed};
					}
				}
			}
		$interval3{$startReg} = $endReg;

		%{ ${$Regions_r}{$chr}{$nm} }  = %interval3;
		%{ ${$NM_Ex_r}{$nm} } = %NM_Ex3;
		%{ $RegBed{$nm} } = %RegBed2;

		}
	}

return(\%RegBed);
}


#########################
#extract bed from file, and intersect with %Regions (regions of selected genes) 
#no -L option: shorten bed eventually
#create 1 hash :
#%{ $RegBed{"raw"} } = linkBed(\%{ $Regions{"raw"} },\%Bed);

sub linkBed {

my($interval,$Bed_r) = @_;
			
my%RegBed;			# $RegBed{$NM}{$startRegion}{$startBed} = $endBed;

foreach my$chr (keys%{$interval}) {

	my@bedStarts = sort{$a<=>$b}keys%{ ${$Bed_r}{$chr} };
	my$i1 = 0;		# idx of @bedStarts

	my%NMs = ();
	foreach my$NM (keys%{ ${$interval}{$chr} }) {
		@{ $NMs{$NM}{"starts"} } = sort{$a<=>$b}(keys%{ ${$interval}{$chr}{$NM} });
		}
	my@NMs = sort { $NMs{$a}{"starts"}[0]<=>$NMs{$b}{"starts"}[0] } (keys%NMs);	##sort en fonction du pemrier start
	
	foreach my$NM (@NMs) {

		while ( ($i1 < $#bedStarts) && (${$Bed_r}{$chr}{$bedStarts[$i1]} < $NMs{$NM}{"starts"}[0]) ) { $i1++; }

		my$c = 0;		# idx of @{ $NMs{$NM}{"starts"} }
		my$i2 = $i1;			# idx of @bedStarts within NM loop

		while ($i2 < scalar@bedStarts) {

			my$bedEnd = ${$Bed_r}{$chr}{$bedStarts[$i2]};

			while ( ($c < (scalar@{ $NMs{$NM}{"starts"} } - 1)) && ($bedStarts[$i2] > ${$interval}{$chr}{$NM}{$NMs{$NM}{"starts"}[$c]}) ) { $c++; }

			my$c2 = $c;			# idx of @{ $NMs{$NM}{"starts"} within i2 loop
			while ( ($c2 < scalar@{ $NMs{$NM}{"starts"} }) && ($bedStarts[$i2] <= ${$interval}{$chr}{$NM}{$NMs{$NM}{"starts"}[$c2]}) && ($bedEnd >= $NMs{$NM}{"starts"}[$c2]) ) { 
				#eventually restrict %Bed
				if ($bedStarts[$i2] < $NMs{$NM}{"starts"}[$c2]) {
					if ($bedEnd > ${$interval}{$chr}{$NM}{$NMs{$NM}{"starts"}[$c2]}) { 
						$RegBed{$NM}{$NMs{$NM}{"starts"}[$c2]}{$NMs{$NM}{"starts"}[$c2]} = ${$interval}{$chr}{$NM}{$NMs{$NM}{"starts"}[$c2]};
						}
					else { 
						$RegBed{$NM}{$NMs{$NM}{"starts"}[$c2]}{$NMs{$NM}{"starts"}[$c2]} = $bedEnd;
						}
					}
				else {
					if ($bedEnd > ${$interval}{$chr}{$NM}{$NMs{$NM}{"starts"}[$c2]}) { 
						$RegBed{$NM}{$NMs{$NM}{"starts"}[$c2]}{$bedStarts[$i2]} = ${$interval}{$chr}{$NM}{$NMs{$NM}{"starts"}[$c2]};
						}
					else {
						$RegBed{$NM}{$NMs{$NM}{"starts"}[$c2]}{$bedStarts[$i2]} = $bedEnd;
						}
					}
				$c2++;
				}
			$i2++;
			}
		}
	}

return(\%RegBed);
}


#########################
#eliminates $NM with no corresponding cov bed
#if ($All) { notAnalysedG($geneNM_r,\%Genes,\%{ $Regions{"raw"} },\%{ $NM_Ex{"raw"} },\%{ $RegBed{"raw"} }); }

sub notAnalysedG {

my($geneNM_r,$Gene_r,$Region_r,$NM_Ex_r,$RegBed_r) = @_;
=pod
my(@Not,%Genes2,%Regions2,%NM_Ex2,%RegBed2);

foreach my$gene (keys%{$Gene_r}) {
	my@startBed = keys%{ ${$RegBed_r}{$gene} };
	unless (@startBed) {
#		print "$gene not in analysed bed\n";
		push(@Not ,$gene);
		}
	}
foreach my$chr (sort(keys%{$Region_r})) {
	foreach my$gene (sort(keys%{ ${$Region_r}{$chr} })) {
		my$ok=1;
		foreach my$not(@Not) {
			if ($gene eq $not)
				{ $ok=0; last; }
			}
		if ($ok) {
			$Genes2{$gene} = 1;
			%{ $Regions2{$chr}{$gene} } = %{ ${$Region_r}{$chr}{$gene} };
			foreach (@{ ${$geneNM_r}{$gene}})
				{ %{ $NM_Ex2{$_} } = %{ ${$NM_Ex_r}{$_} }; }
			%{ $RegBed2{$gene} } = %{ ${$RegBed_r}{$gene} };
			}	
		}
	}

%{$Gene_r} = %Genes2;
%{$Region_r} = %Regions2;
%{$NM_Ex_r} = %NM_Ex2;
%{$RegBed_r} = %RegBed2;
=cut

foreach my$chr (sort(keys%{$Region_r})) {
	foreach my$gene (sort(keys%{ ${$Region_r}{$chr} })) {
		if (scalar(keys%{ ${$RegBed_r}{$gene} }) == 0) {
			delete(${$Gene_r}{$gene});
			delete(${$Region_r}{$chr}{$gene});
			foreach (@{ ${$geneNM_r}{$gene}}) { delete(${$NM_Ex_r}{$_}); }
			delete(${$RegBed_r}{$gene});
			}
		}
	}

}


#############
#eliminates $NM with no corresponding cov bed
#else  { notAnalysedN(\%Genes,\%NMgene,\%{ $Regions{"raw"} },\%{ $NM_Ex{"raw"} },\%{ $RegBed{"raw"} }); }

sub notAnalysedN {

my($Gene_r,$NMgene_r,$Region_r,$NM_Ex_r,$RegBed_r) = @_;
=pod
my(@Not,%Genes2,%Regions2,%NM_Ex2,%RegBed2);

foreach my$NM (keys%{$NM_Ex_r}) {
	my@startBed = keys%{ ${$RegBed_r}{$NM} };
	unless (@startBed) {
#		print "$NM not in analysed bed\n";
		push(@Not ,$NM);
		}
	}
foreach my$chr (sort(keys%{$Region_r})) {
	foreach my$NM (sort(keys%{ ${$Region_r}{$chr} })) {
		my$ok=1;
		foreach my$not(@Not) {
			if ($NM eq $not)
				{ $ok=0; last; }
			}
		if ($ok) {
			$Genes2{${$NMgene_r}{$NM}} = 1;
			%{ $Regions2{$chr}{$NM} } = %{ ${$Region_r}{$chr}{$NM} };
			%{ $NM_Ex2{$NM} } = %{ ${$NM_Ex_r}{$NM} };
			%{ $RegBed2{$NM} } = %{ ${$RegBed_r}{$NM} };
			}	
		}
	}

%{$Gene_r} = %Genes2;
%{$Region_r} = %Regions2;
%{$NM_Ex_r} = %NM_Ex2;
%{$RegBed_r} = %RegBed2;
=cut

foreach my$chr (sort(keys%{$Region_r})) {
	foreach my$NM (sort(keys%{ ${$Region_r}{$chr} })) {
		if (scalar(keys%{ ${$RegBed_r}{$NM} }) == 0) {
			delete(${$Gene_r}{${$NMgene_r}{$NM}});
			delete(${$Region_r}{$chr}{$NM});
			delete(${$NM_Ex_r}{$NM});
			delete(${$RegBed_r}{$NM});
			}
		}
	}

}


#########################
## %{ $RegMut{"raw"} } = linkMut(\%{ $Regions{"raw"} },\%Mut);
sub linkMut {

my($NM_interval_r,$Mut_r) = @_;
# ${$NM_interval_r}{chr}{NM}{start of region} = end of region
# ${$Mut_r}{$chr}{$startMut} = $infoMut;

my%RegMut;				#$RegMut{$NM}{$startRegion}{$startMut} = $MutInfo;
my$c=0;				#idx of $Starts = sort{$a<=>$b}(keys%{ ${$NM_interval_r}{$chr} });
foreach my$chr (keys%{$NM_interval_r}) {
	foreach my$NM (keys%{ ${$NM_interval_r}{$chr} }) {
		my%interval = %{ ${$NM_interval_r}{$chr}{$NM} };
		my@Starts = sort{$a<=>$b}(keys%interval);
		my$c=0;
		foreach my$mut (sort{$a<=>$b}keys%{ ${$Mut_r}{$chr} }) {
			if ( $mut < $Starts[$c] )
				{ next; }
			while ( ($c < (scalar(@Starts)-1)) && ($mut > $interval{$Starts[$c]}) )
				{ $c++; }
			if ( ($mut <= $interval{$Starts[$c]}) && ($mut >= $Starts[$c]) )
				{ $RegMut{$NM}{$Starts[$c]}{$mut} = ${$Mut_r}{$chr}{$mut}; }
			}
		}
	}

return(\%RegMut);

}


#########################
#transposition and splicing, for exons of transcripts: 
#1st position of 1st exon = 0, 
#if end of exon = x, start of next exon = x

#@hashSub = transposeReg($spacer,\%{ $Regions{"raw"} });
#$spacer = $hashSub[0];
#%NMlength = %{$hashSub[1]};			#$NMlength{$NM} = $end of region (for start of region = 0)
#%{ $Regions{"coord0"} } = %{$hashSub[2]};	#$Regions{"coord0"}{$NM}{start of region} = $end of region (for start of region = 0)

sub transposeReg {

my($spacer,$Region_r) = @_;
my(%Reg_00,%NMlength);

#$Reg_00{$NM}{$startReg} = endReg (for start of region = 0)
foreach my$chr (keys%{$Region_r}) {
	#$gene or $NM
	foreach my$gene ( keys%{ ${$Region_r}{$chr} } ) {
#		print "transpose and splice $gene\n";
		my@Starts = sort{$a<=>$b}(keys%{ ${$Region_r}{$chr}{$gene} });
		#find $spacer as $NMlength/100
		#unless((exists $opts{S}) && ($opts{S} eq "N")) {
			$NMlength{$gene}=0;
			for (my$reg=0;$reg<scalar@Starts;$reg++)
				{ $NMlength{$gene} += (${$Region_r}{$chr}{$gene}{$Starts[$reg]} - $Starts[$reg] +1); }
			$spacer = int($NMlength{$gene}/100);
		#	}
		#transpose $Regions{$chr}{$gene}
		my$start00 = 0;
		my$end00 =  ${$Region_r}{$chr}{$gene}{$Starts[0]} - $Starts[0] + $start00 +1;
		$Reg_00{$gene}{$start00} = $end00;
		for (my$reg=1;$reg<scalar@Starts;$reg++) {
			$start00 = $end00 + $spacer;
			$end00 = ${$Region_r}{$chr}{$gene}{$Starts[$reg]} - $Starts[$reg] + $start00 +1;
			$Reg_00{$gene}{$start00} = $end00;
			}
		#re-calculates $NMlength
		my@Starts00 = sort{$a<=>$b}(keys%{ $Reg_00{$gene} });
		$NMlength{$gene} = $Reg_00{$gene}{$Starts00[-1]};
		}
	}

return($spacer,\%NMlength,\%Reg_00);

}

#########################
#transposition and splicing, for exons of transcripts: 
#1st position of 1st exon = 0, 
#if end of exon = x, start of next exon = x

#@hashSub = transposeNM($wUTR,$NMstartCod{"raw"}{$NM},$NMendCod{"raw"}{$NM},\%{ $Regions{"raw"}{$NMchr{$NM}}{$gene} },\%{ $Regions{"coord0"}{$gene} },\%{ $NM_Ex{"raw"}{$NM} }); 
#$NMstartCod{"coord0"}{$NM} = $hashSub[0];
#$NMendCod{"coord0"}{$NM} = $hashSub[1];
#%{ $NM_Ex{"coord0"}{$NM} } = %{$hashSub[2]};

sub transposeNM {

my($wUTR,$NMstartCod,$NMendCod,$Reg_raw,$Reg_00,$NM_Ex_r) = @_;

my(%NM_Ex00,$NMstartCod00,$NMendCod00);	
#NM_Ex00{$startReg}{$startExon} = endExon (for start of region = 0)
my@Starts = sort{$a<=>$b}(keys%{$Reg_raw});
my@Starts00 = sort{$a<=>$b}(keys%{$Reg_00});

my$r=0;
if ($wUTR) {
	#transpose $NMstartCod and $NMendCod
	if (($NMstartCod-1)!=$NMendCod) { 
		while ($NMstartCod > ${$Reg_raw}{$Starts[$r]}) { $r++; }
		$NMstartCod00 = $NMstartCod - $Starts[$r] + $Starts00[$r];
		while ( ($NMendCod > ${$Reg_raw}{$Starts[$r]}) ) { $r++; }
		$NMendCod00 = $NMendCod - $Starts[$r] + $Starts00[$r] +1;
		}
	else { $NMstartCod00=0; $NMendCod00=0; }
	#transpose %NM_Ex
	for ($r=0;$r<scalar@Starts;$r++) {
		foreach my$startEx (keys%{ ${$NM_Ex_r}{$Starts[$r]} }) {
			$NM_Ex00{$Starts00[$r]}{($startEx-$Starts[$r]+$Starts00[$r])} = ${$NM_Ex_r}{$Starts[$r]}{$startEx}-$Starts[$r]+$Starts00[$r]+1;
			}
		}
	}
else {
	my($rIni,$rFin);
	my@startsNM = sort{$a<=>$b}keys%{$NM_Ex_r};
	while ($r <= $#Starts && $Starts[$r] != $startsNM[0]) { $r++; }
	$rIni = $r;
	#transpose $NMstartCod and $NMendCod
	my@StartsEx1 = sort{$a<=>$b}keys%{ ${$NM_Ex_r}{$startsNM[0]} };
	$NMstartCod00 = $StartsEx1[0]-$startsNM[0] + $Starts00[$r];
	$r = $#Starts;
	while ($r > 0 && $Starts[$r] != $startsNM[-1]) { $r--; }
	$rFin = $r;
	my@StartsExf = sort{$a<=>$b}(keys%{ ${$NM_Ex_r}{$startsNM[-1]} });	
	$NMendCod00 = ${$Reg_raw}{$startsNM[-1]}-${$NM_Ex_r}{$startsNM[-1]}{$StartsExf[-1]} + $Starts00[$r] +1;
	#transpose %NM_Ex
	for ($r=$rIni;$r<=$rFin;$r++) {
		foreach my$startEx (keys%{ ${$NM_Ex_r}{$Starts[$r]} }) {
			$NM_Ex00{$Starts00[$r]}{($startEx-$Starts[$r]+$Starts00[$r])} = ${$NM_Ex_r}{$Starts[$r]}{$startEx}-$Starts[$r]+$Starts00[$r]+1;
			}
		}
	}

return($NMstartCod00,$NMendCod00,\%NM_Ex00);

}


#########################
#if ($All):
#make arrays for starts of non coding exons, before and after coding exons
#and arrays for starts of extragenic regions, before and after exons
#@hashSub = designExons1($wUTR,\%Genes,\%geneNM,\%{ $Regions{"coord0"} },\%{ $NM_Ex{"coord0"} },\%{ $NMstartCod{"coord0"} },\%{ $NMendCod{"coord0"} });

sub designExons1 {

my($wUTR,$Gene_r,$geneNM_r,$Reg_00,$NM_Ex_00,$NMstartCod00,$NMendCod00) = @_;
my(%intron,%UTR,%Cod);
#$intron{$NM}{$startInt} = $endInt;
#$UTR{$NM}{$startUTR} = $endUTR;
#$Cod{$NM}{$startCod} = $endCod;

foreach my$gene (keys%{$Gene_r}) {
#	print "design exons of $gene\n";
	my@startReg = sort{$a<=>$b}keys%{ ${$Reg_00}{$gene} } ; 
	foreach my$NM (@{ ${$geneNM_r}{$gene} }) {
#		print "\tfor $NM\n";
		my(%Exons,@startEx);
		#introns:
		foreach my$startR (@startReg) {
			my@Starts = sort{$a<=>$b}(keys%{ ${$NM_Ex_00}{$NM}{$startR} });
			if (@Starts) {
				$intron{$NM}{$startR} = $Starts[0];
				for (my$ex=0;$ex<(scalar@Starts-1);$ex++)
					{ $intron{$NM}{ ${$NM_Ex_00}{$NM}{$startR}{$Starts[$ex]} } = $Starts[$ex+1]; }
				$intron{$NM}{ ${$NM_Ex_00}{$NM}{$startR}{$Starts[-1]} } = ${$Reg_00}{$gene}{$startR};
				push (@startEx , @Starts); 
				for (my$ex=0;$ex<scalar@Starts;$ex++)
					{ $Exons{$Starts[$ex]} = ${$NM_Ex_00}{$NM}{$startR}{$Starts[$ex]}; }
				}
			else
				{ $intron{$NM}{$startR} = ${$Reg_00}{$gene}{$startR}; }
			}
		#exons:
		my$i=0;
		if ($wUTR) {
			if (${$NMstartCod00}{$NM} == ${$NMendCod00}{$NM}) {
				while ($i<scalar@startEx) {
					$UTR{$NM}{$startEx[$i]} = $Exons{$startEx[$i]};
					$i++;
					}
				}
			else {
					#UTRpreCod
				while ( $Exons{$startEx[$i]} < ${$NMstartCod00}{$NM} && $i<(scalar@startEx-1) ) {
					$UTR{$NM}{$startEx[$i]} = $Exons{$startEx[$i]};
					$i++;						
					}
				if ( $startEx[$i] < ${$NMstartCod00}{$NM} ) 
					{ $UTR{$NM}{$startEx[$i]} = ${$NMstartCod00}{$NM}; }
					#Cod
				if( $Exons{$startEx[$i]} >= ${$NMendCod00}{$NM} )
					{ $Cod{$NM}{${$NMstartCod00}{$NM}} = ${$NMendCod00}{$NM}; }
				else {
					$Cod{$NM}{${$NMstartCod00}{$NM}} = $Exons{$startEx[$i]}; 
					$i++;
					while ( $Exons{$startEx[$i]} < ${$NMendCod00}{$NM} && $i<(scalar@startEx-1) ) {
						$Cod{$NM}{$startEx[$i]} = $Exons{$startEx[$i]};
						$i++;	
						}
					if ( $Exons{$startEx[$i]} >= ${$NMendCod00}{$NM} ) 
						{ $Cod{$NM}{$startEx[$i]} = ${$NMendCod00}{$NM}; }
					}
					#UTRpostCod
				if ( $Exons{$startEx[$i]} > ${$NMendCod00}{$NM} ) 
					{ $UTR{$NM}{${$NMendCod00}{$NM}} = $Exons{$startEx[$i]}; }
				while ($i<(scalar@startEx-1)) {
					$i++;
					$UTR{$NM}{$startEx[$i]} = $Exons{$startEx[$i]};
					}
				}
			}
		else {
			foreach (@startEx)
				{ $Cod{$NM}{$_} = $Exons{$_}; }
			}
		}
	}

return(\%intron,\%UTR,\%Cod);

}

#########################
#foreach $NM
#make arrays for starts of non coding exons, before and after coding exons
#and arrays for starts of extragenic regions, before and after exons
#@hashSub = designExons2($wUTR,\%{ $Regions{"coord0"} },\%{ $NM_Ex{"coord0"} },\%{ $NMstartCod{"coord0"} },\%{ $NMendCod{"coord0"} });

sub designExons2 {

my($wUTR,$Reg_00,$NM_Ex00,$NMstartCod00,$NMendCod00) = @_;
my(%intron,%UTR,%Cod);
#$intron{$NM}{$startInt} = $endInt;
#$UTR{$NM}{$startUTR} = $endUTR;
#$Cod{$NM}{$startCod} = $endCod;

foreach my$NM (keys%{$Reg_00}) {
#	print "design exons of $NM\n";
	my@startReg = sort{$a<=>$b}keys%{ ${$Reg_00}{$NM} } ;
	my(%Exons,@startEx);
	#introns:
	foreach my$startR (@startReg) {
		my@Starts = sort{$a<=>$b}(keys%{ ${$NM_Ex00}{$NM}{$startR} });
		if (@Starts)
			{
			$intron{$NM}{$startR} = $Starts[0];
			for (my$ex=0;$ex<(scalar@Starts-1);$ex++)
				{ $intron{$NM}{ ${$NM_Ex00}{$NM}{$startR}{$Starts[$ex]} } = $Starts[$ex+1]; }
			$intron{$NM}{ ${$NM_Ex00}{$NM}{$startR}{$Starts[-1]} } = ${$Reg_00}{$NM}{$startR};
			push (@startEx , @Starts);
			for (my$ex=0;$ex<scalar@Starts;$ex++)
				{ $Exons{$Starts[$ex]} = ${$NM_Ex00}{$NM}{$startR}{$Starts[$ex]}; }
			}
		else
			{ $intron{$NM}{$startR} = ${$Reg_00}{$NM}{$startR}; }
		}
	#exons:
	my$i=0;
	if ($wUTR) {
		if (${$NMstartCod00}{$NM}==${$NMendCod00}{$NM}) {
			while ($i<scalar@startEx) {
				$UTR{$NM}{$startEx[$i]} = $Exons{$startEx[$i]};
				$i++;
				}
			}
		else {
				#UTRpreCod
			while ( $Exons{$startEx[$i]}<${$NMstartCod00}{$NM} && $i<(scalar@startEx-1) ) {
				$UTR{$NM}{$startEx[$i]} = $Exons{$startEx[$i]};
				$i++;						
				}
			if ( $startEx[$i] < ${$NMstartCod00}{$NM} ) 
				{ $UTR{$NM}{$startEx[$i]} = ${$NMstartCod00}{$NM}; }
				#Cod
			if( $Exons{$startEx[$i]} >= ${$NMendCod00}{$NM} )
				{ $Cod{$NM}{${$NMstartCod00}{$NM}} = ${$NMendCod00}{$NM}; }
			else {
				$Cod{$NM}{${$NMstartCod00}{$NM}} = $Exons{$startEx[$i]}; 
				$i++;
				while ( $Exons{$startEx[$i]}<${$NMendCod00}{$NM} && $i<(scalar@startEx-1) ) {
					$Cod{$NM}{$startEx[$i]} = $Exons{$startEx[$i]};
					$i++;	
					}
				if ( $Exons{$startEx[$i]} >= ${$NMendCod00}{$NM} ) 
					{ $Cod{$NM}{$startEx[$i]} = ${$NMendCod00}{$NM}; }
				}
				#UTRpostCod
			if ( $Exons{$startEx[$i]} > ${$NMendCod00}{$NM} ) 
				{ $UTR{$NM}{${$NMendCod00}{$NM}} = $Exons{$startEx[$i]}; }
			while ($i<(scalar@startEx-1)) {
				$i++;
				$UTR{$NM}{$startEx[$i]} = $Exons{$startEx[$i]};
				}
			}
		}
	else {
		foreach (@startEx)
			{ $Cod{$NM}{$_} = $Exons{$_}; }
		}
	}

return(\%intron,\%UTR,\%Cod);

}


#########################
#transposition and splicing of RegBed and RegMut
#@hashSub = transposeBed(\%{ $Regions{"raw"} },\%{ $Regions{"coord0"} },\%{ $RegBed{"raw"} },\%{ $RegMut{"raw"} });
#%{ $RegBed{"coord0"} } = %{$hashSub[0]};
#%{ $RegMut{"coord0"} } = %{$hashSub[1]};

sub transposeBed {

my($Reg_raw,$Reg_00,$RegBed_raw,$RegMut_raw) = @_;

my(%NMbed00,%NMmut00);		#$NMbed00{$NM}{$startBed00} = $endBed00 (start of Region=0 )
					#$NMmut00{$NM}{$mut00} = $infoMut (start of Region=0 )
foreach my$chr (keys%{$Reg_raw}) {
	foreach my$NM ( keys%{ ${$Reg_raw}{$chr} } )	#$gene or $NM
		{
		my@startReg = sort{$a<=>$b}(keys%{ ${$Reg_raw}{$chr}{$NM} });
		my@start00 = sort{$a<=>$b}(keys%{ ${$Reg_00}{$NM} });
		for (my$i=0;$i<scalar@startReg;$i++)
			{
			my$sub = $startReg[$i] - $start00[$i];
			foreach my$startBed ( keys%{ ${$RegBed_raw}{$NM}{$startReg[$i]} } )
				{ $NMbed00{$NM}{($startBed-$sub)} = ${$RegBed_raw}{$NM}{$startReg[$i]}{$startBed}-$sub+1; }
			foreach my$mut ( keys%{ ${$RegMut_raw}{$NM}{$startReg[$i]} } )
				{ $NMmut00{$NM}{($mut-$sub)} = ${$RegMut_raw}{$NM}{$startReg[$i]}{$mut}; }
			}
		}
	}

return(\%NMbed00,\%NMmut00);

}


##########################
#@hashSub = ReverseGene1($gene,$NMlength{$gene},\%{ $Regions{"coord0"}{$gene} },\%{ $RegBed{"coord0"}{$gene} },\%{ $RegMut{"coord0"}{$gene} });
#%{ $Regions{"rev"}{$gene} } = %{$hashSub[0]};
#%{ $RegBed{"rev"}{$gene} } = %{$hashSub[1]};
#%{ $RegMut{"rev"}{$gene} } = %{$hashSub[2]};

sub ReverseGene1 {

my($gene,$length,$Reg_00,$RegBed_00,$RegMut_00) = @_;

#print "reverse $gene\n";
my(%Reg_rev,%NMbed_rev,%NMmut_rev);

foreach my$startReg (keys%{$Reg_00}) {
	$Reg_rev{($length - ${$Reg_00}{$startReg})} = $length - $startReg;
	}

%NMbed_rev = reverseEx($length,$RegBed_00);

#$NMmut00{$NM}{$mut00} = $infoMut (start of Region=0 )
foreach (keys%{$RegMut_00}) {
	$NMmut_rev{($length - $_)} = ${$RegMut_00}{$_};
	}

return(\%Reg_rev,\%NMbed_rev,\%NMmut_rev);

}

##########################
#@hashSub = ReverseGene2($NM,$NMlength{$gene},\%{ $Regions{"coord0"}{$gene} },\%{ $NM_Ex{"coord0"}{$NM} },\%{ $intron{"coord0"}{$NM} },\%{ $UTR{"coord0"}{$NM} },\%{ $Cod{"coord0"}{$NM} });
#%{ $NM_Ex{"rev"}{$NM} } = %{$hashSub[0]};
#%{ $intron{"rev"}{$NM} } = %{$hashSub[1]};
#%{ $UTR{"rev"}{$NM} } = %{$hashSub[2]};
#%{ $Cod{"rev"}{$NM} } = %{$hashSub[3]};

sub ReverseGene2 {

my($NM,$length,$Reg_00,$NM_Ex_00,$intron_00,$UTR_00,$CDS_00) = @_;

#print "reverse $NM\n";
my(%NM_Ex_rev,%intron_rev,%UTR_rev,%CDS_rev);

foreach my$startReg (keys%{$Reg_00}) {
	%{ $NM_Ex_rev{($length-${$Reg_00}{$startReg})} } = reverseEx( $length,\%{ ${$NM_Ex_00}{$startReg} } );
	}
	
# $intron{$NM}{$startIntron} = $endIntron (for start of region = 0)
%intron_rev = reverseEx($length,$intron_00);
# $UTR{$NM}{$startUTR} = $endUTR (for start of region = 0)
%UTR_rev = reverseEx($length,$UTR_00);
# $Cod{$NM}{$startCod} = $endCod (for start of region = 0)}
%CDS_rev = reverseEx($length,$CDS_00);

return(\%NM_Ex_rev,\%intron_rev,\%UTR_rev,\%CDS_rev);

}

##########################
#@hashSub = ReverseNMs($NM,$NMlength{$NM},\%{ $Regions{"coord0"}{$NM} },\%{ $NM_Ex{"coord0"}{$NM} },\%{ $intron{"coord0"}{$NM} },\%{ $UTR{"coord0"}{$NM} },\%{ $Cod{"coord0"}{$NM} },\%{ $RegBed{"coord0"}{$NM} },\%{ $RegMut{"coord0"}{$NM} });
#%{ $Regions{"rev"}{$NM} } = %{$hashSub[0]};
#%{ $NM_Ex{"rev"}{$NM} } = %{$hashSub[1]};
#%{ $intron{"rev"}{$NM} } = %{$hashSub[2]};
#%{ $UTR{"rev"}{$NM} } = %{$hashSub[3]};
#%{ $Cod{"rev"}{$NM} }= %{$hashSub[4]};
#%{ $RegBed{"rev"}{$NM} } = %{$hashSub[5]};
#%{ $RegMut{"rev"}{$NM} } = %{$hashSub[6]};

sub ReverseNMs {

my($NM,$length,$Reg_00,$NM_Ex_00,$intron_00,$UTR_00,$CDS_00,$RegBed_00,$RegMut_00) = @_;

#print "reverse $NM\n";
my(%Reg_rev,%NM_Ex_rev,%intron_rev,%UTR_rev,%CDS_rev,%RegBed_rev,%RegMut_rev);

foreach my$startReg (keys%{$Reg_00}) {
	$Reg_rev{($length - ${$Reg_00}{$startReg})} = $length-$startReg;
	%{ $NM_Ex_rev{($length - ${$Reg_00}{$startReg})} } = reverseEx( $length,\%{ ${$NM_Ex_00}{$startReg} } );
	}

# $intron{$NM}{$startIntron} = $endIntron (for start of region = 0)
%intron_rev = reverseEx($length,$intron_00);
# $UTR{$NM}{$startUTR} = $endUTR (for start of region = 0)
%UTR_rev = reverseEx($length,$UTR_00);
# $Cod{$NM}{$startCod} = $endCod (for start of region = 0)}
%CDS_rev = reverseEx($length,$CDS_00);

#$NMbed00{$NM}{$startBed00} = $endBed00 (start of Region=0 )
%RegBed_rev = reverseEx($length,$RegBed_00);

#$NMmut00{$NM}{$mut00} = $infoMut (start of Region=0 )
foreach (keys%{$RegMut_00}) {
	$RegMut_rev{($length-$_)} = ${$RegMut_00}{$_};
	}

return(\%Reg_rev,\%NM_Ex_rev,\%intron_rev,\%UTR_rev,\%CDS_rev,\%RegBed_rev,\%RegMut_rev);

}

##########################

sub reverseEx {
my($NMlength,$coord_r) = @_;
my(%coordRev);
foreach (keys%{$coord_r}) {
	$coordRev{($NMlength-${$coord_r}{$_})} = $NMlength-$_;
	}
return(%coordRev);
}


##########################


1;




