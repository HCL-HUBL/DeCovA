package DeCovA::covByRegion;


use strict;
use warnings;



#########################

#if bedtools covfile and interval not identical (in previous versions)
#scan bedtools file, and intersect with %allInterval (coding regions of selected genes) 
#create 2 hash : 
# - %allDepth: depth foreach position in intervals from %allInterval : $allCov{$chr}{$loc} = $cov
# - %notCov: start and end of each not covered domain within %allInterval

sub intersectCovFile {

my($Files,$smplIdx,$depthFile,$maxD,$Thresholds_r,$Intervals_r,$doNotCov) = @_;

my@Starts = sort{$a<=>$b}(keys%{$Intervals_r});
my%allDepth;		#for depth line : $allDepth{$loc} = $cov
my$c = 0;			#idx of $startByChr

my$ok = 0;
unless ( -s "$depthFile" ) { die "$depthFile not found or empty\n"; }
open(my$fh, "<", "$depthFile") || die "can't read file $depthFile\n";
while (my$line = <$fh>) {
	chomp $line;
	my@tab = split(/\t/,$line);
	my$loc = $tab[0];
	while ( ($loc > ${$Intervals_r}{$Starts[$c]}) && ($c < (scalar@Starts-1)) ) { $c++; }
	if ( ($loc >= $Starts[$c]) && ($loc <= ${$Intervals_r}{$Starts[$c]}) ) {
		foreach my$f (@{$Files}) {
			my$cov = $tab[${$smplIdx}{$f}];
			if ($maxD) { if ($cov > $maxD) { $cov = $maxD; } }
			$allDepth{$f}{$loc} = $cov;
			}
		}
	if ( $loc > ${$Intervals_r}{$Starts[-1]} )
		{ last; }
	}
close $fh;

#cov-domains:
my%notCov;		#for cov-domains: foreach threshold, for regions < threshold: $notCov{$chr}{$start} = $end
my%notCovStarts;	#for cov-domains: foreach threshold, ordered starts of regions < threshold
if ($doNotCov) {
	my$noCovStart=0; my$noCovEnd=0;
	foreach my$threshold (@{$Thresholds_r}) {
		foreach my$start (@Starts) {
			foreach my$f (@{$Files}) {
				my$pos = $start;
				while ($pos<=${$Intervals_r}{$start}) {
					if ( $allDepth{$f}{$pos} < $threshold ) {
						if ( $pos == ($noCovEnd + 1) ) { 
							$notCov{$f}{$threshold}{$noCovStart} = $pos ;
							$noCovEnd = $pos;
							}
						else { 
							$notCov{$f}{$threshold}{$pos} = $pos ;
							$noCovStart = $pos;
							$noCovEnd = $pos;
							}
						}
					$pos++;
					}
				}
			}
		foreach my$f (@{$Files}) {
			@{ $notCovStarts{$f}{$threshold} } = sort{$a<=>$b}keys%{ $notCov{$f}{$threshold} };
			}
		}
	return(\%allDepth,\%notCov,\%notCovStarts);
	}
else { return(\%allDepth); }

}


#########################

# covByThreshold($threshold,\%{ $Regions{"raw"}{$chr}{$gene} },\%allDepth,\%{ $NMCov{$gene}{$threshold} });
sub covByThreshold {

my($threshold,$Regions_r,$Depth_r,$Cov_r) = @_;
# ${$Regions_r}{start of region} = end of region
# ${$Depth_r}{$loc} = $cov
# @{ ${$Cov_r}{$reg} } = [nbr of covered samples foreach pos of reg], for defined threshold
#print"for threshold $threshold, add covered positions\n";
my@Pos = sort{$a<=>$b}keys%{$Depth_r};
my@startReg = sort{$a<=>$b}keys%{$Regions_r};
my$d = 0;		#index of @posByChr=keys%{ ${$Depth_r}{$chr} }
for (my$r=0;$r<scalar@startReg;$r++) {
	my$c=0;		#index of @{ $NMCov{$threshold}{$NM}{$r} }
	while ( ( $startReg[$r] > $Pos[$d] ) && ( $d < $#Pos ) ) { $d++; }
	while ( ( ${$Regions_r}{$startReg[$r]} >= $Pos[$d] ) && ( $d < $#Pos ) ) {
		if ((${$Depth_r}{$Pos[$d]} ne "") && (${$Depth_r}{$Pos[$d]} >= $threshold)) { ${$Cov_r}{$r}[$c]++; }
		$d++; $c++;
		}
	}

}


#########################

#for cov-line :
#create a hash : foreach file, foreach NM, foreach exon, return an array of values for ordered positions
# %NMcov: key = $file, value = %(key = $NM, value = %( key = $ex, value = @(cov foreach ordered bp of exon) ) )
sub depthLine {

my($Regions_r,$allDepth_r)=@_; 
# ${$Regions_r}{start of region} = end of region
# ${$allDepth_r}{$loc} = $cov

my%depthLine;		#@{ $depthLine{$startReg} } = [$allDepth{$posByChr[$c]}, ]
my@Pos = sort{$a<=>$b}keys%{$allDepth_r};
my@startReg = sort{$a<=>$b}keys%{$Regions_r};
my$c=0;
for (my$r=0;$r<scalar@startReg;$r++) { 
	while ( ( $startReg[$r] > $Pos[$c] ) && ( $c < (scalar@Pos -1) ) )
		{ $c++; }
	while ( ( ${$Regions_r}{$startReg[$r]} >= $Pos[$c] ) && ( $c < (scalar@Pos -1) ) ) {
		push( @{ $depthLine{$r} } , ${$allDepth_r}{$Pos[$c]} );
		$c++;
		}
	}
return(\%depthLine);

}	


########################

#for cov-domains:
#intersections not covered regions and exons 
	#starts and ends of not covered regions within exons
sub notCovDomains1 {

#@hashSub = notCovDomains1($All,$fName{$file},$outdir,$outfile,$gene,$chr,$NM,$NMsens{$NM},\@Thresholds,$pThreshold,\%{ $Regions{$chr}{$gene} },\%{ $NM_Ex{$NM} },\%notCov,\%notCovStarts,\%NMnotCov);
my($All,$file2,$sName2,$outdir,$outfile,$gene,$chr,$NM,$sens,$Thresholds_r,$printThreshold,$Regions_r,$NM_Ex_r,$NMstartCod,$NMendCod,$notCov_r,$notCovStarts_r,$NMnotCov_r,$printReports) = @_;
# ${$Regions_r}{start of region} = end of region
# ${$NM_Ex_r}{start of region}{start of exon} = end of exon
# ${$notCov_r}{$threshold}{$start} = $end
# ${$notCovStarts_r}
# ${$NMnotCov_r}

my(%NMcovStart,%NMcovEnd);
#$NMcovEnd{$NM}{$threshold}{$ex}{$notCovStart} = $notCovEnd;
#@{ $NMcovStart{$NM}{$threshold}{$ex} } = [ordered $notCovStarts]

my$nExon=0;	#total nber of $Exons
foreach my$startR (keys%{$NM_Ex_r})
	{ $nExon+=scalar(keys%{ ${$NM_Ex_r}{$startR} }); }
my@startReg = sort{$a<=>$b}keys%{$Regions_r};

#my$txt1 = "\nsample: $sName2\ngene: $gene\ntranscript: $NM\ncoding sequence: chr $chr : $NMstartCod - $NMendCod\n";
my$txt1 = "\ngene: $gene\ntranscript: $NM\ncoding sequence: chr $chr : $NMstartCod - $NMendCod\n";
my%txt;
my($geneL,%geneUC,%NM_UC);
foreach my$threshold (@{$Thresholds_r}) {

	#$txt{$threshold} .= "\ndepth threshold = $threshold\n\n";
	my$nEx=0;	#idx of $Exons (for printing)
	my$posix=0;	#idx of $notCovStarts
	$geneUC{$threshold}=0; $geneL=0;
	for (my$r=0;$r<scalar@startReg;$r++)
		{
		my$startEx = $startReg[$r]; 
		my$endEx = ${$Regions_r}{$startReg[$r]};
		$geneL += ${$Regions_r}{$startReg[$r]}-$startReg[$r];

		unless ( scalar@{ ${$notCovStarts_r}{$threshold} } == 0 )
			{
			while ( ( $startEx > ${$notCov_r}{$threshold}{${$notCovStarts_r}{$threshold}[$posix]} ) && ( $posix < (scalar@{ ${$notCovStarts_r}{$threshold} }-1) ) )
				{ $posix++; }
	
			while ( ( $endEx >= ${$notCovStarts_r}{$threshold}[$posix] ) && ( $posix < (scalar@{ ${$notCovStarts_r}{$threshold} }-1) ) )
				{
				if ($startEx >= ${$notCovStarts_r}{$threshold}[$posix]) 
					{
					if ($endEx > ${$notCov_r}{$threshold}{${$notCovStarts_r}{$threshold}[$posix]})
						{ $NMcovEnd{$threshold}{$r}{$startEx} = ${$notCov_r}{$threshold}{${$notCovStarts_r}{$threshold}[$posix]}; }
					else	#$endEx <= ${$notCov_r}{$chr}{${$notCovStarts_r}{$chr}[$posix]}
						{ $NMcovEnd{$threshold}{$r}{$startEx} = $endEx; }
					}
				else	#$startEx < ${$notCovStarts_r}{$chr}[$posix]
					{
					if ($endEx <= ${$notCov_r}{$threshold}{${$notCovStarts_r}{$threshold}[$posix]})
						{ $NMcovEnd{$threshold}{$r}{${$notCovStarts_r}{$threshold}[$posix]} = $endEx; }
					else	#$end > $notCov{$chr}{${$notCovStarts_r}{$chr}[$posix]}
						{ $NMcovEnd{$threshold}{$r}{${$notCovStarts_r}{$threshold}[$posix]} = ${$notCov_r}{$threshold}{${$notCovStarts_r}{$threshold}[$posix]}; }
					}
				$posix++;
				}
			#for last [posix]:
			if ( ($startEx <= ${$notCov_r}{$threshold}{${$notCovStarts_r}{$threshold}[$posix]}) && ($endEx >= ${$notCovStarts_r}{$threshold}[$posix]) )
				{
				if ($startEx >= ${$notCovStarts_r}{$threshold}[$posix]) 
					{
					if ($endEx > ${$notCov_r}{$threshold}{${$notCovStarts_r}{$threshold}[$posix]})
						{ $NMcovEnd{$threshold}{$r}{$startEx} = ${$notCov_r}{$threshold}{${$notCovStarts_r}{$threshold}[$posix]}; }
					else	#$endEx <= ${$notCov_r}{$chr}{${$notCovStarts_r}{$chr}[$posix]}
						{ $NMcovEnd{$threshold}{$r}{$startEx} = $endEx; }
					}
				else	#$startEx < ${$notCovStarts_r}{$chr}[$posix]
					{
					if ($endEx <= ${$notCov_r}{$threshold}{${$notCovStarts_r}{$threshold}[$posix]})
						{ $NMcovEnd{$threshold}{$r}{${$notCovStarts_r}{$threshold}[$posix]} = $endEx; }
					else	#$end > ${$notCov_r}{$chr}{${$notCovStarts_r}{$chr}[$posix]}
						{ $NMcovEnd{$threshold}{$r}{${$notCovStarts_r}{$threshold}[$posix]} = ${$notCov_r}{$threshold}{${$notCovStarts_r}{$threshold}[$posix]}; }
					}
				}
			}
		if ($printThreshold) {
			if ($sens eq "+") {
				foreach my$ex (sort{$a<=>$b}keys%{ ${$NM_Ex_r}{$startReg[$r]} }) {
					$nEx++;
					if (exists $NMcovEnd{$threshold}{$r})
						{ $txt{$threshold} .= "exon $nEx : ".$ex." - ".${$NM_Ex_r}{$startReg[$r]}{$ex}."\n"; } 
					}
				}
			else {
				foreach my$ex (sort{$a<=>$b}keys%{ ${$NM_Ex_r}{$startReg[$r]} }) {
					if (exists $NMcovEnd{$threshold}{$r}) 
						{ $txt{$threshold} .= "exon ".($nExon-$nEx)." : ".$ex." - ".${$NM_Ex_r}{$startReg[$r]}{$ex}."\n"; }
					$nEx++; 
					}
				}
			if ((exists $NMcovEnd{$threshold}{$r}) && (scalar(keys%{ ${$NM_Ex_r}{$startReg[$r]} })!=0)) {
				${$NMnotCov_r}{$NM}{$file2}{$threshold} = 1;
				$txt{$threshold} .= "\tnot covered from :\n";
				foreach my$pos (sort{$a<=>$b}(keys%{ $NMcovEnd{$threshold}{$r} })) { 
					push(@{ $NMcovStart{$threshold}{$r} }, $pos);
					$txt{$threshold} .= "\t$pos\tto\t".$NMcovEnd{$threshold}{$r}{$pos}."\n";
					$geneUC{$threshold} += $NMcovEnd{$threshold}{$r}{$pos} - $pos;
					}
				}
			}
		else {
			if ($sens eq "+") {
				foreach my$ex ( sort{$a<=>$b}keys%{ ${$NM_Ex_r}{$startReg[$r]} } ) {
					$nEx++;
					$txt{$threshold} .= "exon $nEx : ".$ex." - ".${$NM_Ex_r}{$startReg[$r]}{$ex}."\n"; 
					}
				}
			else {
				foreach my$ex ( sort{$a<=>$b}keys%{ ${$NM_Ex_r}{$startReg[$r]} } ) { 
					$txt{$threshold} .= "exon ".($nExon-$nEx)." : ".$ex." - ".${$NM_Ex_r}{$startReg[$r]}{$ex}."\n"; 
					$nEx++; 
					}
				}
			if ((exists $NMcovEnd{$threshold}{$r}) && (scalar(keys%{ ${$NM_Ex_r}{$startReg[$r]} })!=0)) {
				$txt{$threshold} .= "\tnot covered from :\n";
				foreach my$pos (sort{$a<=>$b}(keys%{ $NMcovEnd{$threshold}{$r} })) { 
					push(@{ $NMcovStart{$threshold}{$r} }, $pos);
					$txt{$threshold} .= "\t$pos\tto\t".$NMcovEnd{$threshold}{$r}{$pos}."\n";
					$geneUC{$threshold} += $NMcovEnd{$threshold}{$r}{$pos} - $pos;
					}
				}
			}
		}
	$NM_UC{$threshold} = ($geneUC{$threshold}/$geneL)*100;
#	print "$gene\_$NM\t".sprintf("%.1f",$NM_UC{$threshold})."\n";
	}
if($printReports) {
	if ($printThreshold) {
		if (exists ${$NMnotCov_r}{$NM}{$file2}{$printThreshold}) {
			#open(OUT, ">$outdir/cov\_$sName2/$outfile$sName2\_$gene\_$NM.txt") || die "can't create file $outdir/$outfile$file2\_$NM.txt\n";
			open(OUT, ">>", "$outdir/cov\_$sName2/$outfile$sName2\_geneReport.txt") || die "can't create file $outdir/cov\_$sName2/$outfile$sName2\_geneReport.txt\n";
			print OUT $txt1;
			#foreach my$threshold (@{$Thresholds_r}) {
			#	print OUT "\ndepth threshold = $threshold\n\n";
			#	print OUT "not covered over ".sprintf("%.1f",$NM_UC{$threshold})."% of exons\n\n";
			#	if (exists $txt{$threshold}) { print OUT $txt{$threshold}; }
			#	}
			print OUT "not covered over ".sprintf("%.1f",$NM_UC{$printThreshold})."% of exons\n\n";
			if (exists $txt{$printThreshold}) { print OUT $txt{$printThreshold}; }
			close OUT;
			}
		}
	else {
		#open(OUT, ">$outdir/cov\_$sName2/$outfile$sName2\_$gene\_$NM.txt") || die "can't create file $outdir/$outfile$file2\_$NM.txt\n";
		open(OUT, ">>", "$outdir/cov\_$sName2/$outfile$sName2\_geneReport.txt") || die "can't create file $outdir/cov\_$sName2/$outfile$sName2\_geneReport.txt\n";
		print OUT $txt1;
		foreach my$threshold (@{$Thresholds_r}) {
			print OUT "\ndepth threshold = $threshold\n\n";
			print OUT "not covered over ".sprintf("%.1f",$NM_UC{$threshold})."% of exons:\n\n";
			print OUT $txt{$threshold};
			}
		 close OUT;
		}
	}

return(\%NMcovStart,\%NMcovEnd,\%NM_UC);

}


########################

#for cov-domains:
#intersections not covered regions and exons 
	#starts and ends of not covered regions within exons
# @hashSub = notCovDomains2(\@Thresholds,\%{ $Regions{"raw"}{$chr}{$gene} },\%notCov,\%notCovStarts);

sub notCovDomains2 {

my($Thresholds_r,$Regions_r,$notCov_r,$notCovStarts_r) = @_;
# ${$Thresholds_r}
# ${$Regions_r}{start of region} = end of region
# ${$notCov_r}{$threshold}{$start} = $end
# ${$notCovStarts_r}

my(%NMcovStart, %NMcovEnd);
#$NMcovEnd{$NM}{$threshold}{$ex}{$notCovStart} = $notCovEnd;
#@{ $NMcovStart{$NM}{$threshold}{$ex} } = [ordered $notCovStarts]

my@startReg = sort{$a<=>$b}keys%{$Regions_r};
foreach my$threshold (@{$Thresholds_r})
	{
	my$posix=0;	#idx of $notCovStarts
	for (my$r=0;$r<scalar@startReg;$r++)
		{
		my$startEx = $startReg[$r]; 
		my$endEx = ${$Regions_r}{$startReg[$r]}; 

		unless ( scalar@{ ${$notCovStarts_r}{$threshold} } == 0 )
			{
			while ( ( $startEx > ${$notCov_r}{$threshold}{${$notCovStarts_r}{$threshold}[$posix]} ) && ( $posix < (scalar@{ ${$notCovStarts_r}{$threshold} }-1) ) )
				{ $posix++; }
	
			while ( ( $endEx >= ${$notCovStarts_r}{$threshold}[$posix] ) && ( $posix < (scalar@{ ${$notCovStarts_r}{$threshold} }-1) ) )
				{
				if ($startEx >= ${$notCovStarts_r}{$threshold}[$posix]) 
					{
					if ($endEx > ${$notCov_r}{$threshold}{${$notCovStarts_r}{$threshold}[$posix]})
						{ $NMcovEnd{$threshold}{$r}{$startEx} = ${$notCov_r}{$threshold}{${$notCovStarts_r}{$threshold}[$posix]}; }
					else	#$endEx <= ${$notCov_r}{$chr}{${$notCovStarts_r}{$chr}[$posix]}
						{ $NMcovEnd{$threshold}{$r}{$startEx} = $endEx; }
					}
				else	#$startEx < ${$notCovStarts_r}{$chr}[$posix]
					{
					if ($endEx <= ${$notCov_r}{$threshold}{${$notCovStarts_r}{$threshold}[$posix]})
						{ $NMcovEnd{$threshold}{$r}{${$notCovStarts_r}{$threshold}[$posix]} = $endEx; }
					else	#$end > ${$notCov_r}{$chr}{${$notCovStarts_r}{$chr}[$posix]}
						{ $NMcovEnd{$threshold}{$r}{${$notCovStarts_r}{$threshold}[$posix]} = ${$notCov_r}{$threshold}{${$notCovStarts_r}{$threshold}[$posix]}; }
					}
				$posix++;
				}
			#for last [posix]:
			if ( ($startEx <= ${$notCov_r}{$threshold}{${$notCovStarts_r}{$threshold}[$posix]}) && ($endEx >= ${$notCovStarts_r}{$threshold}[$posix]) )
				{
				if ($startEx >= ${$notCovStarts_r}{$threshold}[$posix]) 
					{
					if ($endEx > ${$notCov_r}{$threshold}{${$notCovStarts_r}{$threshold}[$posix]})
						{ $NMcovEnd{$threshold}{$r}{$startEx} = ${$notCov_r}{$threshold}{${$notCovStarts_r}{$threshold}[$posix]}; }
					else	#$endEx <= ${$notCov_r}{$chr}{${$notCovStarts_r}{$chr}[$posix]}
						{ $NMcovEnd{$threshold}{$r}{$startEx} = $endEx; }
					}
				else	#$startEx < ${$notCovStarts_r}{$chr}[$posix]
					{
					if ($endEx <= ${$notCov_r}{$threshold}{${$notCovStarts_r}{$threshold}[$posix]})
						{ $NMcovEnd{$threshold}{$r}{${$notCovStarts_r}{$threshold}[$posix]} = $endEx; }
					else	#$end > ${$notCov_r}{$chr}{${$notCovStarts_r}{$chr}[$posix]}
						{ $NMcovEnd{$threshold}{$r}{${$notCovStarts_r}{$threshold}[$posix]} = ${$notCov_r}{$threshold}{${$notCovStarts_r}{$threshold}[$posix]}; }
					}
				}
			}
		if (exists $NMcovEnd{$threshold}{$r})
			{
			foreach my$pos (sort{$a<=>$b}(keys%{ $NMcovEnd{$threshold}{$r} }))
				{ push(@{ $NMcovStart{$threshold}{$r} }, $pos); }
			}
		}
	}
close OUT;
return(\%NMcovStart,\%NMcovEnd);

}


#########################

#transposition and splicing, for %NMcovStart and %NMcovEnd
#for $gene and $NM
#@hashSub = transposeCov(\@Thresholds,\%{ $Regions{"raw"}{$chr}{$gene} },\%{ $Regions{"coord0"}{$gene} },\%{ $NMcovStart{$gene} },\%{ $NMcovEnd{$gene} });

sub transposeCov {

my($Thresholds_r,$Regions_r,$Reg_00_r,$NMcovStart_r,$NMcovEnd_r) = @_;
# ${$Thresholds_r}
# ${$Regions_r}{start of region} = end of region
# ${$Reg_00_r}{start of region} = end of region (start of region =0)
# ${$NMcovStart_r}
# ${$NMcovEnd_r}

my(%covStart01,%covEnd01);

my@startReg = sort{$a<=>$b}keys%{$Regions_r};
my@start00 = sort{$a<=>$b}keys%{$Reg_00_r};
foreach my$threshold (@{$Thresholds_r})
	{
	my$sub = $startReg[0];
	for (my$i=0;$i<scalar(keys%{ ${$NMcovEnd_r}{$threshold}{0} });$i++)
		{ $covEnd01{$threshold}{0}{(${$NMcovStart_r}{$threshold}{0}[$i]-$sub)} = ${$NMcovEnd_r}{$threshold}{0}{${$NMcovStart_r}{$threshold}{0}[$i]}-$sub+1; }

	for (my$ex=1;$ex<scalar@startReg;$ex++)
		{
		$sub = $startReg[$ex] - $start00[$ex];
		for (my$i=0;$i<scalar(keys%{ ${$NMcovEnd_r}{$threshold}{$ex} });$i++)
			{ $covEnd01{$threshold}{$ex}{(${$NMcovStart_r}{$threshold}{$ex}[$i]-$sub)} = ${$NMcovEnd_r}{$threshold}{$ex}{${$NMcovStart_r}{$threshold}{$ex}[$i]}-$sub+1; }
		}

	for (my$ex=0;$ex<scalar@startReg;$ex++)
		{
		foreach my$pos ( sort{$a<=>$b}(keys%{ $covEnd01{$threshold}{$ex} }) )
			{  push(@{ $covStart01{$threshold}{$ex} }, $pos); }
		}
	}

return(\%covStart01,\%covEnd01);

}


############################

#foreach NM, find max cov, among all samples
# @hashSub = maxCov(\@Files,\%{ $NMdepth{$gene} },\%{ $Regions{"coord0"}{$gene} });

sub maxCov {

my($Files_r,$NMdepth_r,$Reg_00_r) = @_;
# @{ ${$NMdepth_r}{$file}{$startReg} } = [ cov foreach ordered bp of exon ]
# ${$Reg_00_r}{start of region} = $end of region (for start of region = 0)

my$maxCovA = 0;		#for all samples
my(%maxCovS);		#for each sample
foreach my$file (@{$Files_r}) {
	my$maxD = 0;
	for (my$ex=0;$ex<scalar(keys%{$Reg_00_r});$ex++) {
		foreach my$cov (@{ ${$NMdepth_r}{$file}{$ex} }) {
			if ($cov > $maxD) { $maxD = $cov; }
			if ($cov > $maxCovA) { $maxCovA = $cov; }
			}
		}
	$maxCovS{$file} = $maxD;
	}
return($maxCovA,\%maxCovS);

}


##########################

# ReverseCov_Gene($Sum,$allS,$byS,$gene,$NMlength{$gene},\@Thresholds,\%{ $NMdepth{$gene} },\%{ $covEnd01{$gene} },\%{ $covStart01{$gene} },\%{ $NMCov{$gene} });

sub ReverseCov_Gene {

my($Sum,$allS,$byS,$gene,$length,$Thresholds_r,$NMdepth_r,$covEnd01_r,$covStart01_r,$NMcov_r) = @_;

#print "reverse cov for $gene\n";

use DeCovA::geneRegion;

if ($allS || $byS) {		
	# @{ ${$NMdepth_r}{$file}{$startReg} } = [ cov foreach ordered bp of exon ]
	foreach my$file (keys%{$NMdepth_r}) {
		my$nEx=scalar(keys%{ ${$NMdepth_r}{$file} })-1;
		my%tmp;
		foreach my$ex (sort{$a<=>$b}(keys%{ ${$NMdepth_r}{$file} })) {
			@{ $tmp{($nEx-$ex)} } = reverse@{ ${$NMdepth_r}{$file}{$ex} };
			}
		%{ ${$NMdepth_r}{$file} } = %tmp;
		}
	# ${$covEnd01_r}{$file}{$threshold}{$ex}{$notCovStart00} = $notCovEnd00;
	# @{ ${$covStart01_r}{$file}{$threshold}{$ex} } = [ordered $notCovStarts00]
	foreach my$file (keys%{$covEnd01_r}) {
		foreach my$threshold (keys%{ ${$covEnd01_r}{$file} }) {
			my%tmp;
			my$nEx=scalar(keys%{ ${$covEnd01_r}{$file}{$threshold} })-1;
			foreach my$ex (keys%{ ${$covEnd01_r}{$file}{$threshold} }) {
				%{ $tmp{($nEx-$ex)} } = DeCovA::geneRegion::reverseEx($length,\%{ ${$covEnd01_r}{$file}{$threshold}{$ex} });
				}
			%{ ${$covEnd01_r}{$file}{$threshold} } = %tmp;
			foreach my$ex (keys%{ ${$covEnd01_r}{$file}{$threshold} }) {
				@{ ${$covStart01_r}{$file}{$threshold}{$ex} } = sort{$a<=>$b}keys%{ ${$covEnd01_r}{$file}{$threshold}{$ex} };
				}
			}
		}
	}

if ($Sum) {
	# @{ ${$NMcov_r}{$threshold}{$NM}{$nReg} } = [nber of covered samples foreach pos of Regions{$NM}{$startReg}]
	foreach my$threshold (@{$Thresholds_r}) {
		my$nEx = scalar(keys%{ ${$NMcov_r}{$threshold} })-1;
		my%tmp;
		foreach my$ex (sort{$a<=>$b}(keys%{ ${$NMcov_r}{$threshold} })) {
			@{ $tmp{($nEx-$ex)} } = reverse@{ ${$NMcov_r}{$threshold}{$ex} };
			}
		%{ ${$NMcov_r}{$threshold} } = %tmp;
		}
	}

}

##########################

# ReverseCov_NMs($Sum,$allS,$byS,$NM,$NMlength{$NM},\@Thresholds,\%{ $NMdepth{$NM} },\%{ $covEnd01{$NM} },\%{ $covStart01{$NM} },\%{ $NMCov{$NM} });

sub ReverseCov_NMs {

my($Sum,$allS,$byS,$NM,$length,$Thresholds_r,$NMdepth_r,$covEnd01_r,$covStart01_r,$NMCov_r) = @_;

#print "reverse cov for $NM\n";

use DeCovA::geneRegion;

if ($allS || $byS) {		
	# @{ ${$NMdepth_r}{$file}{$startReg} } = [ cov foreach ordered bp of exon ]
	foreach my$file (keys%{$NMdepth_r}) {
		my$nEx=scalar(keys%{ ${$NMdepth_r}{$file} })-1;
		my%tmp;
		foreach my$ex (sort{$a<=>$b}(keys%{ ${$NMdepth_r}{$file} })) {
			@{ $tmp{($nEx-$ex)} } = reverse@{ ${$NMdepth_r}{$file}{$ex} };
			}
		%{ ${$NMdepth_r}{$file} } = %tmp;
		}
	# ${$covEnd01_r}{$file}{$NM}{$threshold}{$ex}{$notCovStart00} = $notCovEnd00;
	# @{ ${$covStart01_r}{$file}{$NM}{$threshold}{$ex} } = [ordered $notCovStarts00]
	foreach my$file (keys%{$covEnd01_r}) {
		foreach my$threshold (keys%{ ${$covEnd01_r}{$file} }) {
			my%tmp;
			my$nEx=scalar(keys%{ ${$covEnd01_r}{$file}{$threshold} })-1;
			foreach my$ex (keys%{ ${$covEnd01_r}{$file}{$threshold} }) {
				%{ $tmp{($nEx-$ex)} } = DeCovA::geneRegion::reverseEx($length,\%{ ${$covEnd01_r}{$file}{$threshold}{$ex} });
				}
			%{ ${$covEnd01_r}{$file}{$threshold} } = %tmp;
			foreach my$ex (keys%{ ${$covEnd01_r}{$file}{$threshold} }) {
				@{ ${$covStart01_r}{$file}{$threshold}{$ex} } = sort{$a<=>$b}keys%{ ${$covEnd01_r}{$file}{$threshold}{$ex} };
				}
			}
		}
	}

if ($Sum) {
	# @{ ${$NMCov_r}{$threshold}{$NM}{$nReg} } = [nber of covered samples foreach pos of Regions{$NM}{$startReg}]
	foreach my$threshold (@{$Thresholds_r}) {
		my$nEx = scalar(keys%{ ${$NMCov_r}{$threshold} })-1;
		my%tmp;
		foreach my$ex (sort{$a<=>$b}(keys%{ ${$NMCov_r}{$threshold} })) {
			@{ $tmp{($nEx-$ex)} } = reverse@{ ${$NMCov_r}{$threshold}{$ex} };
			}
		%{ ${$NMCov_r}{$threshold} } = %tmp;
		}
	}

}


##########################

#@hashSub= covDomains(\@Thresholds,\%{ $Reg_00{$NM} },\%{ $NMCov{$NM} });
#$NMCovEnd{$NM}{$threshold}{$ex}{$start}=$end
#$NMCovVal{$NM}{$threshold}{$ex}{$start}=$value
sub covDomains {

my($Thresholds_r,$Reg_00_r,$NMCov_r) = @_;
# ${$Reg_00_r}{start of region} = $end of region (for start of region = 0)
# @{ ${$NMCov_r}{$threshold}{$nReg} } = [nber of covered samples foreach pos of Regions{$NM}{$startReg}]

my(%NMCovEnd,%NMCovVal);
# $NMCovEnd{$threshold}{$ex}{$start}=$end
# $NMCovVal{$threshold}{$ex}{$start}=$value

my@startR = sort{$a<=>$b}keys%{$Reg_00_r};
foreach my$threshold (@{$Thresholds_r}) {
	foreach my$r (sort{$a<=>$b}(keys%{ ${$NMCov_r}{$threshold} })) {
		my$start = $startR[$r];
		my$end = $start+1; 
		my$cov = ${$NMCov_r}{$threshold}{$r}[0];
		for (my$i=1;$i<scalar@{ ${$NMCov_r}{$threshold}{$r} };$i++) {
			if (${$NMCov_r}{$threshold}{$r}[$i] == $cov)
				{ $end++; }
			else {
				$NMCovEnd{$threshold}{$r}{$start} = $end;
				$NMCovVal{$threshold}{$r}{$start} = $cov;
				$start = $startR[$r]+$i;
				$end = $start+1;
				$cov = ${$NMCov_r}{$threshold}{$r}[$i];
				}
			$NMCovEnd{$threshold}{$r}{$start} = $end;
			$NMCovVal{$threshold}{$r}{$start} = $cov;
			}
		}
	}
return(\%NMCovEnd,\%NMCovVal);

}

##########################



1;


