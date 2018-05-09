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

my($NM_r,$interval1,$NM_Ex_r,$Bed_r) = @_;
my(%interval2,%NM_Ex2,%RegBed1);

my@Starts = sort{$a<=>$b}(keys%{$interval1});
my$c=0;	#$c: count of @StartInterval			
my@startBed = sort{$a<=>$b}keys%{$Bed_r};
my$endBed = ${$Bed_r}{$startBed[0]};
for (my$i=0;$i<scalar@startBed;$i++) {
	if ($startBed[$i] > ${$interval1}{$Starts[-1]})	#if bed keep going on after end of region
		{ last; }
	$endBed = ${$Bed_r}{$startBed[$i]};	
	if ( $endBed < $Starts[0] )
		{ next; }
	if ( $endBed < $Starts[$c] ) {
		if ($startBed[$i] > ${$interval1}{$Starts[$c-1]})
			{ 
			$interval2{$startBed[$i]} = $endBed; 
			$RegBed1{$startBed[$i]}{$startBed[$i]} = $endBed;
			next;
			}
		else { $c--; }
		}
	while ( ($c < (scalar@Starts -1)) && ($startBed[$i] > ${$interval1}{$Starts[$c]}) ) { 
		$interval2{$Starts[$c]} = ${$interval1}{$Starts[$c]};
		foreach my$NM (@{$NM_r}) {
			foreach my$startEx (keys%{ ${$NM_Ex_r}{$NM}{$Starts[$c]} })
				{ $NM_Ex2{$NM}{$Starts[$c]}{$startEx} = ${$NM_Ex_r}{$NM}{$Starts[$c]}{$startEx}; }
			}
		$c++; 
		}
	if ( $startBed[$i] > ${$interval1}{$Starts[$c]} ) { #for current $c
		$interval2{$Starts[$c]} = ${$interval1}{$Starts[$c]};
		foreach my$NM (@{$NM_r}) {
			foreach my$startEx (keys%{ ${$NM_Ex_r}{$NM}{$Starts[$c]} })
				{ $NM_Ex2{$NM}{$Starts[$c]}{$startEx} = ${$NM_Ex_r}{$NM}{$Starts[$c]}{$startEx}; }
			} 
		}
	while ( ($c < (scalar@Starts -1)) && ($startBed[$i] <= ${$interval1}{$Starts[$c]}) && ($endBed >= $Starts[$c]) ) {
		if ($startBed[$i] < $Starts[$c]) {
			if ($endBed > ${$interval1}{$Starts[$c]})
				{ $interval2{$startBed[$i]} = $endBed; }
			else
				{ $interval2{$startBed[$i]} = ${$interval1}{$Starts[$c]}; }
			$RegBed1{$startBed[$i]}{$startBed[$i]} = $endBed;
			foreach my$NM (@{$NM_r}) {
				foreach my$startEx (keys%{ ${$NM_Ex_r}{$NM}{$Starts[$c]} })
					{ $NM_Ex2{$NM}{$startBed[$i]}{$startEx} = ${$NM_Ex_r}{$NM}{$Starts[$c]}{$startEx}; }
				}
			}
		else {
			if ($endBed > ${$interval1}{$Starts[$c]})
				{ $interval2{$Starts[$c]} = $endBed; }
			else
				{ $interval2{$Starts[$c]} = ${$interval1}{$Starts[$c]}; }
			$RegBed1{$Starts[$c]}{$startBed[$i]} = $endBed;
			foreach my$NM (@{$NM_r}) {
				foreach my$startEx (keys%{ ${$NM_Ex_r}{$NM}{$Starts[$c]} })
					{ $NM_Ex2{$NM}{$Starts[$c]}{$startEx} = ${$NM_Ex_r}{$NM}{$Starts[$c]}{$startEx}; }
				}
			}
		$c++;
		}
	if ( ($startBed[$i] <= ${$interval1}{$Starts[$c]}) && ($endBed >= $Starts[$c]) ) {	#for current $c
		if ($startBed[$i] < $Starts[$c]) {
			if ($endBed > ${$interval1}{$Starts[$c]})
				{ $interval2{$startBed[$i]} = $endBed; }
			else
				{ $interval2{$startBed[$i]} = ${$interval1}{$Starts[$c]}; }
			$RegBed1{$startBed[$i]}{$startBed[$i]} = $endBed;
			foreach my$NM (@{$NM_r}) {
				foreach my$startEx (keys%{ ${$NM_Ex_r}{$NM}{$Starts[$c]} })
					{ $NM_Ex2{$NM}{$startBed[$i]}{$startEx} = ${$NM_Ex_r}{$NM}{$Starts[$c]}{$startEx}; }
				}
			}
		else {
			if ($endBed > ${$interval1}{$Starts[$c]})
				{ $interval2{$Starts[$c]} = $endBed; }
			else
				{ $interval2{$Starts[$c]} = ${$interval1}{$Starts[$c]}; }
			$RegBed1{$Starts[$c]}{$startBed[$i]} = $endBed;
			foreach my$NM (@{$NM_r}) {
				foreach my$startEx (keys%{ ${$NM_Ex_r}{$NM}{$Starts[$c]} })
					{ $NM_Ex2{$NM}{$Starts[$c]}{$startEx} = ${$NM_Ex_r}{$NM}{$Starts[$c]}{$startEx}; }
				}
			}
		}
	}
#if region keep going on after end of bed
while ( ($c < scalar@Starts) && ($Starts[$c] > $endBed) ) {
	$interval2{$Starts[$c]} = ${$interval1}{$Starts[$c]};
	foreach my$NM (@{$NM_r}) {
		foreach my$startEx (keys%{ ${$NM_Ex_r}{$NM}{$Starts[$c]} })
			{ $NM_Ex2{$NM}{$Starts[$c]}{$startEx} = ${$NM_Ex_r}{$NM}{$Starts[$c]}{$startEx}; }
		}
	$c++; 
	}

#merge intervals
my(%interval3,%NM_Ex3,%RegBed2);
my@Starts2 = sort{$a<=>$b}(keys%interval2);
my$startReg = $Starts2[0];
my$endReg = $interval2{$startReg};
foreach my$NM (keys%NM_Ex2) {
	foreach my$startEx (keys%{ $NM_Ex2{$NM}{$startReg} })
		{ $NM_Ex3{$NM}{$startReg}{$startEx} = $NM_Ex2{$NM}{$startReg}{$startEx}; }
	}
foreach my$startBed (keys%{ $RegBed1{$startReg} })
	{ $RegBed2{$startReg}{$startBed} = $RegBed1{$startReg}{$startBed}; }
for (my$i=1;$i<scalar(@Starts2);$i++) {
	if ($Starts2[$i] <= $endReg) { 
		if ($interval2{$Starts2[$i]} > $endReg)
			{ $endReg = $interval2{$Starts2[$i]}; }
		foreach my$NM (keys%NM_Ex2) {
			foreach my$startEx (keys%{ $NM_Ex2{$NM}{$Starts2[$i]} })
				{ $NM_Ex3{$NM}{$startReg}{$startEx} = $NM_Ex2{$NM}{$Starts2[$i]}{$startEx}; }
			}
		foreach my$startBed (keys%{ $RegBed1{$Starts2[$i]} })
			{ $RegBed2{$startReg}{$startBed} = $RegBed1{$Starts2[$i]}{$startBed}; }
		}
	else { 
		$interval3{$startReg} = $endReg;
		$startReg = $Starts2[$i];
		$endReg = $interval2{$startReg};
		foreach my$NM (keys%NM_Ex2)
			{
			foreach my$startEx (keys%{ $NM_Ex2{$NM}{$startReg} })
				{ $NM_Ex3{$NM}{$startReg}{$startEx} = $NM_Ex2{$NM}{$startReg}{$startEx}; }
			}
		foreach my$startBed (keys%{ $RegBed1{$startReg} })
			{ $RegBed2{$startReg}{$startBed} = $RegBed1{$startReg}{$startBed}; }
		}
	}
$interval3{$startReg} = $endReg;

%{$interval1} = %interval3;

foreach my$NM (keys%NM_Ex3)
	{ %{ ${$NM_Ex_r}{$NM} } = %{ $NM_Ex3{$NM} }; }

return(\%RegBed2);
}


#########################
#extract bed from file, and intersect with %Regions (regions of selected genes) 
#create 1 hash :
#@hashSub = changeRegion2N(\%{ $Regions{"raw"}{$NMchr{$NM}}{$NM} },\%{ $NM_Ex{"raw"}{$NM} },\%{ $Bed{$NMchr{$NM}} });
#%{ $RegBed{"raw"}{$NM} } = %{$hashSub[0]};	#$RegBed{$NM}{$startReg}{$startBed} = $endBed;

sub changeRegionN {

my($interval1,$NM_Ex_r,$Bed_r) = @_;
my(%interval2,%NM_Ex2,%RegBed1);

my@Starts = sort{$a<=>$b}(keys%{$interval1});
my$c=0;	#$c: count of @StartInterval			
my@startBed = sort{$a<=>$b}keys%{$Bed_r};
my$endBed = ${$Bed_r}{$startBed[0]};
for (my$i=0;$i<scalar@startBed;$i++)
	{
	if ($startBed[$i] > ${$interval1}{$Starts[-1]})	#if bed keep going on after end of region
		{ last; }
	$endBed = ${$Bed_r}{$startBed[$i]};	
	if ( $endBed < $Starts[0] )
		{ next; }
	if ( $endBed < $Starts[$c] )
		{
		if ($startBed[$i] > ${$interval1}{$Starts[$c-1]})
			{ 
			$interval2{$startBed[$i]} = $endBed; 
			$RegBed1{$startBed[$i]}{$startBed[$i]} = $endBed;
			next;
			}
		else { $c--; }
		}
	while ( ($c < (scalar@Starts -1)) && ($startBed[$i] > ${$interval1}{$Starts[$c]}) )
		{ 
		$interval2{$Starts[$c]} = ${$interval1}{$Starts[$c]};
		foreach my$startEx (keys%{ ${$NM_Ex_r}{$Starts[$c]} })
			{ $NM_Ex2{$Starts[$c]}{$startEx} = ${$NM_Ex_r}{$Starts[$c]}{$startEx}; }
		$c++; 
		}
	if ( $startBed[$i] > ${$interval1}{$Starts[$c]} )	#for current $c
		{ 
		$interval2{$Starts[$c]} = ${$interval1}{$Starts[$c]};
		foreach my$startEx (keys%{ ${$NM_Ex_r}{$Starts[$c]} })
			{ $NM_Ex2{$Starts[$c]}{$startEx} = ${$NM_Ex_r}{$Starts[$c]}{$startEx}; }
		}
	while ( ($c < (scalar@Starts -1)) && ($startBed[$i] <= ${$interval1}{$Starts[$c]}) && ($endBed >= $Starts[$c]) )
		{
		if ($startBed[$i] < $Starts[$c])
			{
			if ($endBed > ${$interval1}{$Starts[$c]})
				{ $interval2{$startBed[$i]} = $endBed; }
			else
				{ $interval2{$startBed[$i]} = ${$interval1}{$Starts[$c]}; }
			$RegBed1{$startBed[$i]}{$startBed[$i]} = $endBed;
			foreach my$startEx (keys%{ ${$NM_Ex_r}{$Starts[$c]} })
				{ $NM_Ex2{$startBed[$i]}{$startEx} = ${$NM_Ex_r}{$Starts[$c]}{$startEx}; }
			}
		else
			{
			if ($endBed > ${$interval1}{$Starts[$c]})
				{ $interval2{$Starts[$c]} = $endBed; }
			else
				{ $interval2{$Starts[$c]} = ${$interval1}{$Starts[$c]}; }
			$RegBed1{$Starts[$c]}{$startBed[$i]} = $endBed;
			foreach my$startEx (keys%{ ${$NM_Ex_r}{$Starts[$c]} })
				{ $NM_Ex2{$Starts[$c]}{$startEx} = ${$NM_Ex_r}{$Starts[$c]}{$startEx}; }
			}
		$c++;
		}
	if ( ($startBed[$i] <= ${$interval1}{$Starts[$c]}) && ($endBed >= $Starts[$c]) )	#for current $c
		{
		if ($startBed[$i] < $Starts[$c])
			{
			if ($endBed > ${$interval1}{$Starts[$c]})
				{ $interval2{$startBed[$i]} = $endBed; }
			else
				{ $interval2{$startBed[$i]} = ${$interval1}{$Starts[$c]}; }
			$RegBed1{$startBed[$i]}{$startBed[$i]} = $endBed;
			foreach my$startEx (keys%{ ${$NM_Ex_r}{$Starts[$c]} })
				{ $NM_Ex2{$startBed[$i]}{$startEx} = ${$NM_Ex_r}{$Starts[$c]}{$startEx}; }
			}
		else
			{
			if ($endBed > ${$interval1}{$Starts[$c]})
				{ $interval2{$Starts[$c]} = $endBed; }
			else
				{ $interval2{$Starts[$c]} = ${$interval1}{$Starts[$c]}; }
			$RegBed1{$Starts[$c]}{$startBed[$i]} = $endBed;
			foreach my$startEx (keys%{ ${$NM_Ex_r}{$Starts[$c]} })
				{ $NM_Ex2{$Starts[$c]}{$startEx} = ${$NM_Ex_r}{$Starts[$c]}{$startEx}; }
			}
		}
	}
#if region keep going on after end of bed
while ( ($c < scalar@Starts) && ($Starts[$c] > $endBed) )
	{
	$interval2{$Starts[$c]} = ${$interval1}{$Starts[$c]};
	foreach my$startEx (keys%{ ${$NM_Ex_r}{$Starts[$c]} })
		{ $NM_Ex2{$Starts[$c]}{$startEx} = ${$NM_Ex_r}{$Starts[$c]}{$startEx}; }
	$c++; 
	}

#merge intervals
my(%interval3,%NM_Ex3,%RegBed2);
my@Starts2 = sort{$a<=>$b}(keys%interval2);
my$startReg = $Starts2[0];
my$endReg = $interval2{$startReg};
foreach my$startEx (keys%{ $NM_Ex2{$startReg} })
	{ $NM_Ex3{$startReg}{$startEx} = $NM_Ex2{$startReg}{$startEx}; }
foreach my$startBed (keys%{ $RegBed1{$startReg} })
	{ $RegBed2{$startReg}{$startBed} = $RegBed1{$startReg}{$startBed}; }
for (my$i=1;$i<scalar(@Starts2);$i++)
	{
	if ($Starts2[$i] <= $endReg)
		{ 
		if ($interval2{$Starts2[$i]} > $endReg)
			{ $endReg = $interval2{$Starts2[$i]}; }
		foreach my$startEx (keys%{ $NM_Ex2{$Starts2[$i]} })
			{ $NM_Ex3{$startReg}{$startEx} = $NM_Ex2{$Starts2[$i]}{$startEx}; }
		foreach my$startBed (keys%{ $RegBed1{$Starts2[$i]} })
			{ $RegBed2{$startReg}{$startBed} = $RegBed1{$Starts2[$i]}{$startBed}; }
		}
	else
		{ 
		$interval3{$startReg} = $endReg;
		$startReg = $Starts2[$i];
		$endReg = $interval2{$startReg};
		foreach my$startEx (keys%{ $NM_Ex2{$startReg} })
			{ $NM_Ex3{$startReg}{$startEx} = $NM_Ex2{$startReg}{$startEx}; }
		foreach my$startBed (keys%{ $RegBed1{$startReg} })
			{ $RegBed2{$startReg}{$startBed} = $RegBed1{$startReg}{$startBed}; }
		}
	}
$interval3{$startReg} = $endReg;

%{$interval1} = %interval3;

%{$NM_Ex_r} = %NM_Ex3;

return(\%RegBed2);
}


#########################
#extract bed from file, and intersect with %Regions (regions of selected genes) 
#no -L option: shorten bed eventually
#create 1 hash :
#%{ $RegBed{"raw"} } = linkBed(\%{ $Regions{"raw"} },\%Bed);

sub linkBed {

my($interval,$Bed_r) = @_;

#my%covBed2;			
my%RegBed;			#$RegBed{$NM}{$startRegion}{$startBed} = $endBed;
my$c=0;			#idx of $startByChr = sort{$a<=>$b}(keys%{ $allInterval{$chr} });
foreach my$chr (keys%{$interval}) {
	foreach my$NM (keys%{ ${$interval}{$chr} }) {
		my@Starts = sort{$a<=>$b}(keys%{ ${$interval}{$chr}{$NM} });
		my$c=0;
		foreach my$start (sort{$a<=>$b}keys%{ ${$Bed_r}{$chr} }) {
			my$end = ${$Bed_r}{$chr}{$start};
			if ( $end < $Starts[$c] )
				{ next; }
			while ( ($c < (scalar(@Starts)-1)) && ($start > ${$interval}{$chr}{$NM}{$Starts[$c]}) )
				{ $c++; }
			my$c2=$c;
			while ( ($c2 < scalar@Starts) && ($start <= ${$interval}{$chr}{$NM}{$Starts[$c2]}) && ($end >= $Starts[$c2]) ) { 
				#eventually restrict %Bed
				if ($start < $Starts[$c2]) {
					if ($end > ${$interval}{$chr}{$NM}{$Starts[$c2]}) { 
						#$covBed2{$chr}{$Starts[$c2]} = ${$interval}{$chr}{$NM}{$Starts[$c2]};
						$RegBed{$NM}{$Starts[$c2]}{$Starts[$c2]} = ${$interval}{$chr}{$NM}{$Starts[$c2]};
						}
					else { 
						#$covBed2{$chr}{$Starts[$c2]} = $end;
						$RegBed{$NM}{$Starts[$c2]}{$Starts[$c2]} = $end;
						}
					}
				else
					{
					if ($end > ${$interval}{$chr}{$NM}{$Starts[$c2]}) { 
						#$covBed2{$chr}{$start} = ${$interval}{$chr}{$NM}{$Starts[$c2]}; 
						$RegBed{$NM}{$Starts[$c2]}{$start} = ${$interval}{$chr}{$NM}{$Starts[$c2]};
						}
					else {
						#$covBed2{$chr}{$start} = $end; 
						$RegBed{$NM}{$Starts[$c2]}{$start} = $end;
						}
					}
				$c2++;
				}
			}
		}
	}
#return(\%covBed2,\%RegBed);
return(\%RegBed);
}


#########################
#eliminates $NM with no corresponding cov bed
#if ($All) { notAnalysedG($geneNM_r,\%Genes,\%{ $Regions{"raw"} },\%{ $NM_Ex{"raw"} },\%{ $RegBed{"raw"} }); }

sub notAnalysedG {

my($geneNM_r,$Gene_r,$Region_r,$NM_Ex_r,$RegBed_r) = @_;
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

}

#############
#eliminates $NM with no corresponding cov bed
#else  { notAnalysedN(\%Genes,\%NMgene,\%{ $Regions{"raw"} },\%{ $NM_Ex{"raw"} },\%{ $RegBed{"raw"} }); }

sub notAnalysedN {

my($Gene_r,$NMgene_r,$Region_r,$NM_Ex_r,$RegBed_r) = @_;
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




