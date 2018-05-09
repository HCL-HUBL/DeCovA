package DeCovA::geneGraph;


use strict;
use warnings;



########################
#for $All $NM on same graph
# graphAllSampleG($nGraf,$suff,$outdir,$gene,$sens,$maxGr,$Rev,$chr,$maxCovA{$gene},$NMlength{$gene},\@NMs,\@Files,\@colors,\@Thresholds,$Regions{"raw"},$Regions{"rev"},$NM_Ex{"rev"},$RegBed{"rev"},$intron{"rev"},$UTR{"rev"},$Cod{"rev"},$NMdepth{$gene},$covStart01{$gene},$covEnd01{$gene},$RegMut{"rev"},\%sName2);

sub graphAllSampleG_old {	# = all in 1 sheet, in png
my($nGraf,$suff,$outdir,$gene,$sens,$maxD,$Rev,$chr,$maxY,$NMlength,$NMs_r,$Files_r,$colors_r,$Thresholds_r,$Regions_r,$Reg_00_r,$NM_Ex00_r,$NMbed00_r,$intron_r,$UTR_r,$Cod_r,$NMdepth_r,$covStart01_r,$covEnd01_r,$NMmut00_r,$sampleName_r)= @_;

print "print CMDR graphAllSampleG: $gene\n";

unless ($maxD) {
	if ($maxY<10) { $maxY = 10; }
	}

my$nNMs = scalar@{$NMs_r};
my@Starts = sort{$a<=>$b}(keys%{ $Regions_r->{$chr}{$gene} });
my$startReg = $Starts[0];
my$endReg = $Regions_r->{$chr}{$gene}{$Starts[-1]};

if ($nGraf eq "max") { $nGraf = scalar@{$Files_r}; }
my$cmdR = "";
my$i=1;my$n=1;my$N=1;
foreach my$file (@{$Files_r}) {

	my($Y,$Y1,$Y2);

	#$line1: cmdR for plot cov line
	#@{ $NMdepth{$file}{$ex} } = [ cov foreach ordered bp of exon ]
	$Y1=(0.5+0.25*$nNMs);$Y2=1.25;
	my($line1a,$line1b) = line1($maxY,$NMlength,$Y1,$Y2,\%{ $Reg_00_r->{$gene} },\%{ $NMdepth_r->{$file} });
	
	#$line2: cmdR for cmdR for lines(intron)	#$intron{$startInt} = $endInt;
	my$line2="";
	for(my$i=0;$i<$nNMs;$i++) { 
		$Y=(0.425+0.25*$i);
		$line2 .= line2($maxY,$Y,\%{ $intron_r->{${$NMs_r}[$i]} });
		}

	#$line3a: cmdR for rect(bed, empty)	#$NMbed00{$startBed00} = $endBed00 (start of Region=0 )
	$Y1=0.1;$Y2=0.25;
	my$line3a = line3a($maxY,$Y1,$Y2,\%{ $NMbed00_r->{$gene} });
	my$line3b = line3a($maxY,$Y1,$Y2,\%{ $Reg_00_r->{$gene} });

	#$line3: cmdR for rect(not covered domains)
	#$covEnd01{$file}{$threshold}{$ex}{$notCovStart00} = $notCovEnd00;
	#@{ $covStart01{$file}{$threshold}{$ex} } = [ordered $notCovStarts00]
	$Y1=0.1;$Y2=0.25;
	my%line3 = line3($maxY,$Y1,$Y2,$Thresholds_r,$colors_r,scalar(keys%{ $Reg_00_r->{$gene} }),\%{ $covStart01_r->{$file} },\%{ $covEnd01_r->{$file} });

	#$line4a: cmdR for cmdR for rect(UTR exons)	#$UTR{$NM}{$startUTR} = $endUTR;
	my$line4a="";
	for(my$i=0;$i<$nNMs;$i++) {
		if (scalar(keys%{ $UTR_r->{${$NMs_r}[$i]} }) != 0) {
			$Y1=(0.375+0.25*$i);$Y2=(0.475+0.25*$i);
			$line4a .= line4a($maxY,$Y1,$Y2,\%{ $UTR_r->{${$NMs_r}[$i]} });
			}
		}

	#$line4b: cmdR for rect(coding exons)	#$Cod{$NM}{$startCod} = $endCod;
	my$line4b="";
	for(my$i=0;$i<$nNMs;$i++) {
		if (scalar(keys%{ $Cod_r->{${$NMs_r}[$i]} }) != 0) {
			$Y1=(0.35+0.25*$i);$Y2=(0.5+0.25*$i);
			$line4b .= line4b($maxY,$Y1,$Y2,\%{ $Cod_r->{${$NMs_r}[$i]} });
			}
		}
	
	#$line4c: Nb exon
	#%NM_Ex00{NM}{start of region}{start of exon} = end of exon (for start of region = 0)
	my$line4c="";
	for(my$i=0;$i<$nNMs;$i++) {
		$Y=(0.3+0.25*$i);
		$line4c .= line4c($sens,$Rev,$maxY,$Y,\%{ ${$NM_Ex00_r}{${$NMs_r}[$i]} });
		}
	
	#$line5: cmdR for legends
	my$comp = "<";
	my$line5 = line5a($comp,$colors_r,$Thresholds_r);
	
	#$line6: cmdR for mutations (symbol)
	#$NMmut00{$NM}{$mut00} = $infoMut (start of Region=0 )
	my$line6a = ""; my$line6b = "";
	if (scalar(keys%{ $NMmut00_r->{$gene} }) != 0) {
		$Y1=0.05;$Y2=(0.3+0.25*$nNMs);
		$line6a = line6a($maxY,$Y1,$Y2,\%{ $NMmut00_r->{$gene} });
		#$line6b: cmdR for mutations (text)
		#$NMmut00{$NM}{$mut00} = $infoMut (start of Region=0 )
		$Y=(0.3+0.25*$nNMs);
		$line6b = line6b($maxY,$Y,\%{ $NMmut00_r->{$gene} });
		}

	#$line7: NM names
	my$line7 = line7($NMlength,$maxY,$NMs_r);

	#print CMDR:
	$cmdR .= 
"plot (c(0,0), xlim=c(-".(0.8*$NMlength/10).",".$NMlength."), ylim=c(-".((0.5+0.25*$nNMs)*$maxY).",".(1.25*$maxY)."), type =\"n\", xlab = \"sample : ".${$sampleName_r}{$file}."\", ylab = \"depth\", cex.lab=1.5, cex.axis=1.2, yaxt='n')
axis(2, at=seq(0, $maxY, length = 5), labels=seq(0, $maxY, length = 5), las=2)\n";

	$cmdR .= $line1a.$line1b;

	$cmdR .= 
"abline(h=0, col=\"darkgrey\", lty = \"dotted\", lwd=2)\n";
	my$col=0;
	foreach my$threshold (@{$Thresholds_r}) {
		$cmdR .= 
"abline(h=$threshold, col=\"".${$colors_r}[$col]."\", lty = \"dotted\", lwd=2)\n"; 
		$col++;
		}

	$cmdR .= $line2.$line4a.$line4b.$line4c.$line5.$line6a.$line6b.$line7;

	$cmdR .= $line3b." lwd=1, col=\"green\", border =\"green\")\n";
	foreach my$threshold (@{$Thresholds_r}) {
		if ($line3{$threshold}) 
			{ $cmdR .= $line3{$threshold}; }
		}
	$cmdR .= $line3a." lwd=1, col=NA, border =\"black\")\n";
	
	if ($i==(scalar@{$Files_r}) || $n==$nGraf) {
		open (CMDR, ">$outdir/$gene\_temp.R") || die;
		print CMDR "#!/usr/bin/env Rscript\n\n" ;
		if ($nGraf==scalar@{$Files_r}) { print CMDR
"png(\"".$outdir."/cov_All/".$gene."_bySample.png\", 1500, ".($nGraf*(400*(0.875+(0.125*$nNMs)))).")\n
par(mfrow=c($nGraf,1))\n"; }
		else {
			if ($N>1) { print CMDR
"png(\"".$outdir."/cov_All/".$gene."_bySample_$N.png\", 1500, ".($nGraf*(400*(0.875+(0.125*$nNMs)))).")\n
par(mfrow=c($nGraf,1))\n"; }
			else { print CMDR 
"png(\"".$outdir."/cov_All/".$gene."_bySample_$N.png\", 1500, ".($n*(400*(0.875+(0.125*$nNMs)))).")\n
par(mfrow=c($n,1))\n"; }
			}
		print CMDR "$cmdR";
		if ($sens eq "+") { print CMDR
"mtext(\"$gene\tchr $chr\t".$startReg." >>> ".$endReg."\", side = 3, outer=TRUE, line=-3, cex=1.5)\n"; }
		else { print CMDR 
"mtext(\"$gene\tchr $chr\t".$startReg." <<< ".$endReg."\", side = 3, outer=TRUE, line=-3, cex=1.5)\n"; }
		print CMDR "dev.off();\n";
		close CMDR;
		system "Rscript $outdir/$gene\_temp.R";
		unlink "$outdir/$gene\_temp.R";
		$cmdR="";
		$n=0;
		$N++;
		}

	$i++;$n++;

	}

}

########################
#for $All $NM on same graph
# graphAllSampleG($nGraf,$suff,$outdir,$gene,$sens,$maxGr,$Rev,$chr,$maxCovA{$gene},$NMlength{$gene},\@NMs,\@Files,\@colors,\@Thresholds,$Regions{"raw"},$Regions{"rev"},$NM_Ex{"rev"},$RegBed{"rev"},$intron{"rev"},$UTR{"rev"},$Cod{"rev"},$NMdepth{$gene},$covStart01{$gene},$covEnd01{$gene},$RegMut{"rev"},\%sName2);

sub graphAllSampleG {		#1sample by sheet, 1 doc in pdf

my($nGraf,$suff,$outdir,$gene,$sens,$maxD,$Rev,$chr,$maxY,$NMlength,$NMs_r,$Files_r,$colors_r,$Thresholds_r,$Regions_r,$Reg_00_r,$NM_Ex00_r,$NMbed00_r,$intron_r,$UTR_r,$Cod_r,$NMdepth_r,$covStart01_r,$covEnd01_r,$NMmut00_r,$sampleName_r) = @_;

print "print CMDR graphAllSampleG: $gene\n";

unless ($maxD) {
	if ($maxY<10) { $maxY = 10; }
	}

my$nNMs = scalar@{$NMs_r};
my@Starts = sort{$a<=>$b}(keys%{ $Regions_r->{$chr}{$gene} });
my$startReg = $Starts[0];
my$endReg = $Regions_r->{$chr}{$gene}{$Starts[-1]};

open (CMDR, ">$outdir/$gene\_temp.R") || die;
print CMDR "#!/usr/bin/env Rscript\n\n" ;
print CMDR "pdf(\"".$outdir."/cov_All/".$gene."_bySample.pdf\", width=11.69, height=".(4.135*(0.75+0.875+(0.125*$nNMs))/1.75).")\n";

foreach my$file (@{$Files_r}) {

	my$cmdR = "";

	my($Y,$Y1,$Y2);

	#$line1: cmdR for plot cov line
	#@{ $NMdepth{$file}{$ex} } = [ cov foreach ordered bp of exon ]
	$Y1=(0.5+0.25*$nNMs);$Y2=1.25;
	my($line1a,$line1b) = line1($maxY,$NMlength,$Y1,$Y2,\%{ $Reg_00_r->{$gene} },\%{ $NMdepth_r->{$file} });
	
	#$line2: cmdR for cmdR for lines(intron)	#$intron{$startInt} = $endInt;
	my$line2="";
	for(my$i=0;$i<$nNMs;$i++) { 
		$Y=(0.425+0.25*$i);
		$line2 .= line2($maxY,$Y,\%{ $intron_r->{${$NMs_r}[$i]} });
		}

	#$line3a: cmdR for rect(bed, empty)	#$NMbed00{$startBed00} = $endBed00 (start of Region=0 )
	$Y1=0.1;$Y2=0.25;
	my$line3a = line3a($maxY,$Y1,$Y2,\%{ $NMbed00_r->{$gene} });
	my$line3b = line3a($maxY,$Y1,$Y2,\%{ $Reg_00_r->{$gene} });

	#$line3: cmdR for rect(not covered domains)
	#$covEnd01{$file}{$threshold}{$ex}{$notCovStart00} = $notCovEnd00;
	#@{ $covStart01{$file}{$threshold}{$ex} } = [ordered $notCovStarts00]
	$Y1=0.1;$Y2=0.25;
	my%line3 = line3($maxY,$Y1,$Y2,$Thresholds_r,$colors_r,scalar(keys%{ $Reg_00_r->{$gene} }),\%{ $covStart01_r->{$file} },\%{ $covEnd01_r->{$file} });

	#$line4a: cmdR for cmdR for rect(UTR exons)	#$UTR{$NM}{$startUTR} = $endUTR;
	my$line4a="";
	for(my$i=0;$i<$nNMs;$i++) {
		if (scalar(keys%{ $UTR_r->{${$NMs_r}[$i]} }) != 0) {
			$Y1=(0.375+0.25*$i);$Y2=(0.475+0.25*$i);
			$line4a .= line4a($maxY,$Y1,$Y2,\%{ $UTR_r->{${$NMs_r}[$i]} });
			}
		}

	#$line4b: cmdR for rect(coding exons)	#$Cod{$NM}{$startCod} = $endCod;
	my$line4b="";
	for(my$i=0;$i<$nNMs;$i++) {
		if (scalar(keys%{ $Cod_r->{${$NMs_r}[$i]} }) != 0) {
			$Y1=(0.35+0.25*$i);$Y2=(0.5+0.25*$i);
			$line4b .= line4b($maxY,$Y1,$Y2,\%{ $Cod_r->{${$NMs_r}[$i]} });
			}
		}
	
	#$line4c: Nb exon
	#%NM_Ex00{NM}{start of region}{start of exon} = end of exon (for start of region = 0)
	my$line4c="";
	for(my$i=0;$i<$nNMs;$i++) {
		$Y=(0.3+0.25*$i);
		$line4c .= line4c($sens,$Rev,$maxY,$Y,\%{ ${$NM_Ex00_r}{${$NMs_r}[$i]} });
		}
	
	#$line5: cmdR for legends
	my$comp = "<";
	my$line5 = line5a($comp,$colors_r,$Thresholds_r);
	
	#$line6: cmdR for mutations (symbol)
	#$NMmut00{$NM}{$mut00} = $infoMut (start of Region=0 )
	my$line6a = ""; my$line6b = "";
	if (scalar(keys%{ $NMmut00_r->{$gene} }) != 0) {
		$Y1=0.05;$Y2=(0.3+0.25*$nNMs);
		$line6a = line6a($maxY,$Y1,$Y2,\%{ $NMmut00_r->{$gene} });
		#$line6b: cmdR for mutations (text)
		#$NMmut00{$NM}{$mut00} = $infoMut (start of Region=0 )
		$Y=(0.3+0.25*$nNMs);
		$line6b = line6b($maxY,$Y,\%{ $NMmut00_r->{$gene} });
		}

	#$line7: NM names
	my$line7 = line7($NMlength,$maxY,$NMs_r);

	#print CMDR:
	$cmdR .= "plot (c(0,0), xlim=c(-".(0.8*$NMlength/10).",".$NMlength."), ylim=c(-".((0.5+0.25*$nNMs)*$maxY).",".(1.25*$maxY)."), type =\"n\", xlab = \"sample : ".${$sampleName_r}{$file}."\", ylab = \"Depth\", cex.lab=1, cex.axis=0.8, yaxt='n')
axis(2, at=seq(0, $maxY, length = 5), labels=seq(0, $maxY, length = 5), las=2, cex.axis=0.8)\n";

	$cmdR .= $line1a.$line1b;

	$cmdR .= "abline(h=0, col=\"darkgrey\", lty = \"dotted\", lwd=1.5)\n";
	my$col=0;
	foreach my$threshold (@{$Thresholds_r}) {
		$cmdR .= "abline(h=$threshold, col=\"".${$colors_r}[$col]."\", lty = \"dotted\", lwd=1.5)\n"; 
		$col++;
		}

	$cmdR .= $line2.$line4a.$line4b.$line4c.$line5.$line6a.$line6b.$line7;

	$cmdR .= $line3b." lwd=1, col=\"green\", border =\"green\")\n";
	foreach my$threshold (@{$Thresholds_r}) {
		if ($line3{$threshold}) 
			{ $cmdR .= $line3{$threshold}; }
		}
	$cmdR .= $line3a." lwd=1, col=NA, border =\"black\")\n";
	
	if ($sens eq "+") {
		$cmdR .= "mtext(\"$gene : chr $chr:$startReg >>> $endReg\", side = 3, outer=TRUE, line=-3, cex=1)\n";
		}
	else {
		$cmdR .= "mtext(\"$gene : chr $chr:$startReg <<< $endReg\", side = 3, outer=TRUE, line=-3, cex=1)\n";
		}
	
	print CMDR "$cmdR\n";

	}

print CMDR "dev.off();\n";
close CMDR;
system "Rscript $outdir/$gene\_temp.R";
unlink "$outdir/$gene\_temp.R";

}


########################
#graph foreach $NM 
# graphAllSampleN($nGraf,$NM,$suff,$outdir,$gene,$sens,$maxGr,$Rev,$chr,$maxCovA{$NM},$NMlength{$NM},\@Files,\@colors,\@Thresholds,$Regions{"raw"},$Regions{"rev"},$NM_Ex{"rev"},$RegBed{"rev"},$intron{"rev"},$UTR{"rev"},$Cod{"rev"},$NMdepth{$NM},$covStart01{$NM},$covEnd01{$NM},$RegMut{"rev"},\%sName2);

sub graphAllSampleN {	#still in png, but not used anymore

my($nGraf,$NM,$suff,$outdir,$gene,$sens,$maxD,$Rev,$chr,$maxY,$NMlength,$Files_r,$colors_r,$Thresholds_r,$Regions_r,$Reg_00_r,$NM_Ex00_r,$NMbed00_r,$intron_r,$UTR_r,$Cod_r,$NMdepth_r,$covStart01_r,$covEnd01_r,$NMmut00_r,$sampleName_r)= @_;

print "print CMDR graphAllSampleN: $gene\t$NM\n";

unless ($maxD) {
	if ($maxY<10) { $maxY = 10; }
	}

my@Starts = sort{$a<=>$b}(keys%{ $Regions_r->{$chr}{$NM} });
my$startReg = $Starts[0];
my$endReg = $Regions_r->{$chr}{$NM}{$Starts[-1]};

if ($nGraf eq "max") { $nGraf = scalar@{$Files_r}; }
my$cmdR = "";
my$i=1;my$n=1;my$N=1;
foreach my$file (@{$Files_r}) {

	#$line1: cmdR for plot cov line
	#@{ $NMdepth{$file}{$NM}{$ex} } = [ cov foreach ordered bp of exon ]
	my$Y1=0.75;my$Y2=1.25;
	my($line1a,$line1b) = line1($maxY,$NMlength,$Y1,$Y2,\%{ $Reg_00_r->{$NM} },\%{ $NMdepth_r->{$file} });
	
	#$line2: cmdR for cmdR for lines(intron)	#$intron{$startInt} = $endInt;
	my$Y=0.425;
	my$line2 = line2($maxY,$Y,\%{ $intron_r->{$NM} });

	#$line3a: cmdR for rect(bed, empty)	#$NMbed00{$startBed00} = $endBed00 (start of Region=0 )
	$Y1=0.1;$Y2=0.25;
	my$line3a = line3a($maxY,$Y1,$Y2,\%{ $NMbed00_r->{$NM} });	#empty black squares
	my$line3b = line3a($maxY,$Y1,$Y2,\%{ $Reg_00_r->{$NM} });	#full green squares

	#$line3: cmdR for rect(not covered domains)
	#$covEnd01{$file}{$threshold}{$ex}{$notCovStart00} = $notCovEnd00;
	#@{ $covStart01{$file}{$threshold}{$ex} } = [ordered $notCovStarts00]
	$Y1=0.1;$Y2=0.25;
	my%line3 = line3($maxY,$Y1,$Y2,$Thresholds_r,$colors_r,scalar(keys%{ $Reg_00_r->{$NM} }),\%{ $covStart01_r->{$file} },\%{ $covEnd01_r->{$file} });

	#$line4a: cmdR for cmdR for rect(UTR exons)	#$UTR{$startUTR} = $endUTR;
	$Y1=0.375;$Y2=0.475;
	my$line4a = "";
	if (scalar(keys%{ $UTR_r->{$NM} }) != 0) { $line4a .= line4a($maxY,$Y1,$Y2,\%{ $UTR_r->{$NM} }); }

	#$line4b: cmdR for rect(coding exons)	#$Cod{$startCod} = $endCod;
	$Y1=0.35;$Y2=0.5;
	my$line4b = "";
	if (scalar(keys%{ $Cod_r->{$NM} }) != 0) { $line4b .= line4b($maxY,$Y1,$Y2,\%{ $Cod_r->{$NM} }); }

	#$line4c: Nb exon
	#%NM_Ex00{NM}{start of region}{start of exon} = end of exon (for start of region = 0)
	$Y=0.3;
	my$line4c = line4c($sens,$Rev,$maxY,$Y,\%{ $NM_Ex00_r->{$NM} });
	
	#$line5: cmdR for legends
	my$comp = "<";
	my$line5 = line5a($comp,$colors_r,$Thresholds_r);
	
	#$line6: cmdR for mutations (symbol)
	#$NMmut00{$NM}{$mut00} = $infoMut (start of Region=0 )
	my$line6a = ""; my$line6b = "";
	if (scalar(keys%{ $NMmut00_r->{$NM} }) != 0) {
		$Y1=0.05;$Y2=0.55;
		$line6a = line6a($maxY,$Y1,$Y2,\%{ $NMmut00_r->{$NM} });
		#$line6b: cmdR for mutations (text)
		#$NMmut00{$NM}{$mut00} = $infoMut (start of Region=0 )
		$Y=0.55;
		$line6b = line6b($maxY,$Y,\%{ $NMmut00_r->{$NM} });
		}

	#print CMDR:
	$cmdR .= 
"plot (c(0,0), xlim=c(0,".$NMlength."), ylim=c(-".(0.75*$maxY).",".(1.25*$maxY)."), type =\"n\", xlab = \"sample : ".${$sampleName_r}{$file}."\", ylab = \"depth\", cex.lab=1.5, cex.axis=1.2, yaxt='n')
axis(2, at=seq(0, $maxY, length = 5),labels=seq(0, $maxY, length = 5), las=2)\n";

	$cmdR .= $line1a.$line1b;

	$cmdR .= 
"abline(h=0, col=\"darkgrey\", lty = \"dotted\", lwd=2)\n";
	my$col=0;
	foreach my$threshold (@{$Thresholds_r}) {
		$cmdR .= 
"abline(h=$threshold, col=\"".${$colors_r}[$col]."\", lty = \"dotted\", lwd=2)\n"; 
		$col++;
		}

	$cmdR .= $line2.$line4a.$line4b.$line4c.$line5.$line6a.$line6b;

	$cmdR .= $line3b." lwd=1, col=\"green\", border =\"green\")\n";
	foreach my$threshold (@{$Thresholds_r}) {
		if ($line3{$threshold}) 
			{ $cmdR .= $line3{$threshold}; }
		}
	$cmdR .= $line3a." lwd=1, col=NA, border =\"black\")\n";
	

	if ($i==(scalar@{$Files_r}) || $n==$nGraf) {
		open (CMDR, ">$outdir/$NM\_temp.R") || die;
		print CMDR "#!/usr/bin/env Rscript\n\n" ;
		if ($nGraf==scalar@{$Files_r}) { print CMDR
"png(\"".$outdir."/cov_All/".$gene."_".$NM."_bySample.png\", 1500, ".($nGraf*400).")\n
par(mfrow=c($nGraf,1))\n"; }
		else {
			if ($N>1) { print CMDR
"png(\"".$outdir."/cov_All/".$gene."_".$NM."_bySample_$N.png\", 1500, ".($nGraf*400).")\n
par(mfrow=c($nGraf,1))\n"; }
			else { print CMDR 
"png(\"".$outdir."/cov_All/".$gene."_".$NM."_bySample_$N.png\", 1500, ".($n*400).")\n
par(mfrow=c($n,1))\n"; }
			}
		print CMDR "$cmdR";
		if ($sens eq "+") { print CMDR 
"mtext(\"$gene\t$NM\tchr $chr\t".$startReg." >>> ".$endReg."\", side = 3, outer=TRUE, line=-3, cex=1.5)\n"; }
		else { print CMDR 
"mtext(\"$gene\t$NM\tchr $chr\t".$startReg." <<< ".$endReg."\", side = 3, outer=TRUE, line=-3, cex=1.5)\n"; }
		print CMDR "dev.off();\n";
		close CMDR;
		system "Rscript $outdir/$NM\_temp.R";
		unlink "$outdir/$NM\_temp.R";
		$cmdR="";
		$n=0;
		$N++;
		}

	$i++;$n++;

	}

}


########################
#for $All $NM on same graph
# graphBySampleG($nGraf,$suff,$outdir,$gene,$sens,$maxGr,$Rev,$chr,$NMlength{$gene},$maxCovS{$gene},\@NMs,\@Files,\@colors,\@Thresholds,$Regions{"raw"},$Regions{"rev"},$NM_Ex{"rev"},$RegBed{"rev"},$intron{"rev"},$UTR{"rev"},$Cod{"rev"},$NMdepth{$gene},$covStart01{$gene},$covEnd01{$gene},$RegMut{"rev"},\%sName2);

sub graphBySampleG {

my($nGraf,$suff,$outdir,$gene,$sens,$maxD,$Rev,$chr,$NMlength,$maxCov_r,$NMs_r,$Files_r,$colors_r,$Thresholds_r,$Regions_r,$Reg_00_r,$NM_Ex00_r,$NMbed00_r,$intron_r,$UTR_r,$Cod_r,$NMdepth_r,$covStart01_r,$covEnd01_r,$NMmut00_r,$sampleName_r)= @_;

my$nNMs = scalar@{$NMs_r};
my@Starts = sort{$a<=>$b}(keys%{ $Regions_r->{$chr}{$gene} });
my$startReg = $Starts[0];
my$endReg = $Regions_r->{$chr}{$gene}{$Starts[-1]};

foreach my$file (@{$Files_r}) {

	print "print CMDR graphBySampleG: $gene\n";
	open (CMDR, ">$outdir/$gene\_temp.R") || die;
	print CMDR "#!/usr/bin/env Rscript\n\n" ;
	#print CMDR "png(\"$outdir/cov\_".${$sampleName_r}{$file}."/$gene\_".${$sampleName_r}{$file}.".png\", 1500, ".(400*(0.875+(0.125*$nNMs))).")\n";
	print CMDR "pdf(\"$outdir/cov\_".${$sampleName_r}{$file}."/$gene\_".${$sampleName_r}{$file}.".pdf\", width=11.69, height=".(4.135*(0.75+0.875+(0.125*$nNMs))/1.75).")\n";

	my$maxY = $maxCov_r->{$file};
	unless ($maxD) {
		if ($maxY<10) { $maxY = 10; }
		}

	my($Y,$Y1,$Y2);
	#$line1: cmdR for plot depth line
	#@{ $NMdepth{$file}{$ex} } = [ cov foreach ordered bp of exon ]
	$Y1=(0.5+0.25*$nNMs);$Y2=1.25;
	my($line1a,$line1b) = line1($maxY,$NMlength,$Y1,$Y2,\%{ $Reg_00_r->{$gene} },\%{ $NMdepth_r->{$file} });
	
	#$line2: cmdR for cmdR for lines(intron)	#$intron{$startInt} = $endInt;
	my$line2="";
	for(my$i=0;$i<$nNMs;$i++) { 
		$Y=(0.425+0.25*$i);
		$line2 .= line2($maxY,$Y,\%{ $intron_r->{${$NMs_r}[$i]} });
		}

	#$line3a: cmdR for rect(bed, empty)	#$NMbed00{$startBed00} = $endBed00 (start of Region=0 )
	$Y1=0.1;$Y2=0.25;
	my$line3a = line3a($maxY,$Y1,$Y2,\%{ $NMbed00_r->{$gene} });
	my$line3b = line3a($maxY,$Y1,$Y2,\%{ $Reg_00_r->{$gene} });

	#$line3: cmdR for rect(not covered domains)
	#$covEnd01{$file}{$threshold}{$ex}{$notCovStart00} = $notCovEnd00;
	#@{ $covStart01{$file}{$threshold}{$ex} } = [ordered $notCovStarts00]
	$Y1=0.1;$Y2=0.25;
	my%line3 = line3($maxY,$Y1,$Y2,$Thresholds_r,$colors_r,scalar(keys%{ $Reg_00_r->{$gene} }),\%{ $covStart01_r->{$file} },\%{ $covEnd01_r->{$file} });

	#$line4a: cmdR for cmdR for rect(UTR exons)	#$UTR{$NM}{$startUTR} = $endUTR;
	my$line4a="";
	for(my$i=0;$i<$nNMs;$i++) {
		if (scalar(keys%{ $UTR_r->{${$NMs_r}[$i]} }) != 0) {
			$Y1=(0.375+0.25*$i);$Y2=(0.475+0.25*$i);
			$line4a .= line4a($maxY,$Y1,$Y2,\%{ $UTR_r->{${$NMs_r}[$i]} });
			}
		}

	#$line4b: cmdR for rect(coding exons)	#$Cod{$NM}{$startCod} = $endCod;
	my$line4b="";
	for(my$i=0;$i<$nNMs;$i++) {
		if (scalar(keys%{ $Cod_r->{${$NMs_r}[$i]} }) != 0) {
			$Y1=(0.35+0.25*$i);$Y2=(0.5+0.25*$i);
			$line4b .= line4b($maxY,$Y1,$Y2,\%{ $Cod_r->{${$NMs_r}[$i]} });
			}
		}
	
	#$line4c: Nb exon
	#%NM_Ex00{NM}{start of region}{start of exon} = end of exon (for start of region = 0)
	my$line4c="";
	for(my$i=0;$i<$nNMs;$i++) {
		$Y=(0.3+0.25*$i);
		$line4c .= line4c($sens,$Rev,$maxY,$Y,\%{ ${$NM_Ex00_r}{${$NMs_r}[$i]} });
		}
	
	#$line5: cmdR for legends
	my$comp = "<";
	my$line5 = line5a($comp,$colors_r,$Thresholds_r);
	
	#$line6: cmdR for mutations (symbol)
	#$NMmut00{$NM}{$mut00} = $infoMut (start of Region=0 )
	my$line6a = ""; my$line6b = "";
	if (scalar(keys%{ $NMmut00_r->{$gene} }) != 0) {
		$Y1=0.05;$Y2=(0.3+0.25*$nNMs);
		$line6a = line6a($maxY,$Y1,$Y2,\%{ $NMmut00_r->{$gene} });
		#$line6b: cmdR for mutations (text)
		#$NMmut00{$NM}{$mut00} = $infoMut (start of Region=0 )
		$Y=(0.3+0.25*$nNMs);
		$line6b = line6b($maxY,$Y,\%{ $NMmut00_r->{$gene} });
		}

	#$line7: NM names
	my$line7 = line7($NMlength,$maxY,$NMs_r);

	#print CMDR:
	print CMDR "plot (c(0,0), xlim=c(-".(0.8*$NMlength/10).",".$NMlength."), ylim=c(-".((0.5+0.25*$nNMs)*$maxY).",".(1.25*$maxY)."), type =\"n\", xlab = \"sample : ".${$sampleName_r}{$file}."\", ylab = \"depth\", cex.lab=1, cex.axis=0.8, yaxt='n')
axis(2, at=seq(0, $maxY, length = 5),labels=seq(0, $maxY, length = 5), las=2, cex.axis=0.8)\n";
	#png:"plot (c(0,0), xlim=c(-".(0.8*$NMlength/10).",".$NMlength."), ylim=c(-".((0.5+0.25*$nNMs)*$maxY).",".(1.25*$maxY)."), type =\"n\", xlab = \"sample : ".${$sampleName_r}{$file}."\", ylab = \"depth\", cex.lab=1.5, cex.axis=1.2, yaxt='n')
#axis(2, at=seq(0, $maxY, length = 5),labels=seq(0, $maxY, length = 5), las=2, cex.axis=1.2)\n";

	print CMDR $line1a.$line1b;

	print CMDR "abline(h=0, col=\"darkgrey\", lty = \"dotted\", lwd=1.5)\n";
	my$col=0;
	foreach my$threshold (@{$Thresholds_r}) {
		print CMDR "abline(h=$threshold, col=\"".${$colors_r}[$col]."\", lty = \"dotted\", lwd=1.5)\n"; 
		$col++;
		}

	print CMDR $line2.$line4a.$line4b.$line4c.$line5.$line6a.$line6b.$line7;

	print CMDR $line3b." lwd=1, col=\"green\", border =\"green\")\n";
	foreach my$threshold (@{$Thresholds_r}) {
		if ($line3{$threshold}) 
			{ print CMDR $line3{$threshold}; }
		}
	print CMDR $line3a." lwd=1, col=NA, border =\"black\")\n";
	
	if ($sens eq "+") {
		#png: print CMDR "mtext(\"$gene\tchr $chr\t".$startReg." >>> ".$endReg."\", side = 3, outer=TRUE, line=-3, cex=1.5)\n";
		print CMDR "mtext(\"$gene : chr $chr:$startReg >>> $endReg\", side = 3, outer=TRUE, line=-3, cex=1)\n";
		}
	else {
		#png: print CMDR "mtext(\"$gene\tchr $chr\t".$startReg." <<< ".$endReg."\", side = 3, outer=TRUE, line=-3, cex=1.5)\n";
		print CMDR "mtext(\"$gene : chr $chr:$startReg <<< $endReg\", side = 3, outer=TRUE, line=-3, cex=1)\n";
		}

	print CMDR "dev.off();\n";
	close CMDR;
	system "Rscript $outdir/$gene\_temp.R";
	unlink "$outdir/$gene\_temp.R";
	}

}


########################
#graph foreach $NM
# graphBySampleN($nGraf,$NM,$suff,$outdir,$gene,$sens,$maxGr,$Rev,$chr,$NMlength{$NM},$maxCovS{$NM},\@Files,\@colors,\@Thresholds,$Regions{"raw"},$Regions{"rev"},$NM_Ex{"rev"},$RegBed{"rev"},$intron{"rev"},$UTR{"rev"},$Cod{"rev"},$NMdepth{$NM},$covStart01{$NM},$covEnd01{$NM},$RegMut{"rev"},\%sName2);

sub graphBySampleN {	#still in png, but not used anymore

my($nGraf,$NM,$suff,$outdir,$gene,$sens,$maxD,$Rev,$chr,$NMlength,$maxCov_r,$Files_r,$colors_r,$Thresholds_r,$Regions_r,$Reg_00_r,$NM_Ex00_r,$NMbed00_r,$intron_r,$UTR_r,$Cod_r,$NMdepth_r,$covStart01_r,$covEnd01_r,$NMmut00_r,$sampleName_r)= @_;

my@Starts = sort{$a<=>$b}(keys%{ $Regions_r->{$chr}{$NM} });
my$startReg = $Starts[0];
my$endReg = $Regions_r->{$chr}{$NM}{$Starts[-1]};

foreach my$file (@{$Files_r}) {

	print "print CMDR graphBySampleN: $gene\t$NM\n";
	open (CMDR, ">$outdir/$NM\_temp.R") || die;
	print CMDR "#!/usr/bin/env Rscript\n\n" ;
	print CMDR 
	"png(\"$outdir/cov\_".${$sampleName_r}{$file}."/$gene\_$NM\_".${$sampleName_r}{$file}.".png\", 1500, 400)\n"; 

	my$maxY = $maxCov_r->{$file};
	unless ($maxD) {
		if ($maxY<10) { $maxY = 10; }
		}

	#$line1: cmdR for plot cov line
	#@{ $NMdepth{$file}{$NM}{$ex} } = [ cov foreach ordered bp of exon ]
	my$Y1=0.75;my$Y2=1.25;
	my($line1a,$line1b) = line1($maxY,$NMlength,$Y1,$Y2,\%{ $Reg_00_r->{$NM} },\%{ $NMdepth_r->{$file} });
	
	#$line2: cmdR for cmdR for lines(intron)	#$intron{$NM}{$startInt} = $endInt;
	my$Y=0.425;
	my$line2 = line2($maxY,$Y,\%{ $intron_r->{$NM} });

	#$line3a: cmdR for rect(bed, empty)	#$NMbed00{$startBed00} = $endBed00 (start of Region=0 )
	$Y1=0.1;$Y2=0.25;
	my$line3a = line3a($maxY,$Y1,$Y2,\%{ $NMbed00_r->{$NM} });	#empty black squares
	my$line3b = line3a($maxY,$Y1,$Y2,\%{ $Reg_00_r->{$NM} });	#full green squares

	#$line3: cmdR for rect(not covered domains)
	#$covEnd01{$file}{$threshold}{$ex}{$notCovStart00} = $notCovEnd00;
	#@{ $covStart01{$file}{$threshold}{$ex} } = [ordered $notCovStarts00]
	$Y1=0.1;$Y2=0.25;
	my%line3 = line3($maxY,$Y1,$Y2,$Thresholds_r,$colors_r,scalar(keys%{ $Reg_00_r->{$NM} }),\%{ $covStart01_r->{$file} },\%{ $covEnd01_r->{$file} });

	#$line4a: cmdR for cmdR for rect(UTR exons)	#$UTR{$NM}{$startUTR} = $endUTR;
	$Y1=0.375;$Y2=0.475;
	my$line4a = "";
	if (scalar(keys%{ $UTR_r->{$NM} }) !=0 ) { $line4a .= line4a($maxY,$Y1,$Y2,\%{ $UTR_r->{$NM} }); }

	#$line4b: cmdR for rect(coding exons)	#$Cod{$NM}{$startCod} = $endCod;
	$Y1=0.35;$Y2=0.5;
	my$line4b = "";
	if (scalar(keys%{ $Cod_r->{$NM} }) !=0 ) { $line4b .= line4b($maxY,$Y1,$Y2,\%{ $Cod_r->{$NM} }); }

	#$line4c: Nb exon
	#%NM_Ex00{NM}{start of region}{start of exon} = end of exon (for start of region = 0)
	$Y=0.3;
	my$line4c = line4c($sens,$Rev,$maxY,$Y,\%{ $NM_Ex00_r->{$NM} });
	
	#$line5: cmdR for legends
	my$comp = "<";
	my$line5 = line5a($comp,$colors_r,$Thresholds_r);
	
	#$line6: cmdR for mutations (symbol)
	#$NMmut00{$NM}{$mut00} = $infoMut (start of Region=0 )
	my$line6a = ""; my$line6b = "";
	if (scalar(keys%{ $NMmut00_r->{$NM} }) != 0) {
		$Y1=0.05;$Y2=0.55;
		$line6a = line6a($maxY,$Y1,$Y2,\%{ $NMmut00_r->{$NM} });
		#$line6b: cmdR for mutations (text)
		#$NMmut00{$NM}{$mut00} = $infoMut (start of Region=0 )
		$Y=0.55;
		$line6b = line6b($maxY,$Y,\%{ $NMmut00_r->{$NM} });
		}

	#print CMDR:
	print CMDR 
"plot (c(0,0), xlim=c(0,".$NMlength."), ylim=c(-".(0.75*$maxY).",".(1.25*$maxY)."), type =\"n\", xlab = \"sample : ".${$sampleName_r}{$file}."\", ylab = \"depth\", cex.lab=1.5, cex.axis=1.2, yaxt='n')
axis(2, at=seq(0, $maxY, length = 5), labels=seq(0, $maxY, length = 5), las=2)\n";

	print CMDR $line1a.$line1b;


	print CMDR 
"abline(h=0, col=\"darkgrey\", lty = \"dotted\", lwd=2)\n"; 
	my$col=0;
	foreach my$threshold (@{$Thresholds_r}) {
		print CMDR 
"abline(h=$threshold, col=\"".${$colors_r}[$col]."\", lty = \"dotted\", lwd=2)\n"; 
		$col++;
		}

	print CMDR $line2.$line4a.$line4b.$line4c.$line5.$line6a.$line6b;

	print CMDR $line3b." lwd=1, col=\"green\", border =\"green\")\n";
	foreach my$threshold (@{$Thresholds_r}) {
		if ($line3{$threshold}) 
			{ print CMDR $line3{$threshold}; }
		}
	print CMDR $line3a." lwd=1, col=NA, border =\"black\")\n";
	
	if ($sens eq "+") {
		print CMDR 
	"mtext(\"$gene\t$NM\tchr $chr\t".$startReg." >>> ".$endReg."\", side = 3, outer=TRUE, line=-3, cex=1.5)\n";
		}
	else {
		print CMDR 
	"mtext(\"$gene\t$NM\tchr $chr\t".$startReg." <<< ".$endReg."\", side = 3, outer=TRUE, line=-3, cex=1.5)\n";
		}
	print CMDR "dev.off();\n";
	close CMDR;
	system "Rscript $outdir/$NM\_temp.R";
	unlink "$outdir/$NM\_temp.R";
	}

}


############################
#$line1: cmdR for plot cov line
#@{ $NMdepth{$NM}{$file}{$ex} } = [ cov foreach ordered bp of exon ]
#my($line1a,$line1b) = line1($maxY,$NMlength,$Y1,$Y2,\%{ $Reg_00_r->{$gene} },\%{ $NMdepth_r->{$file} });

sub line1 {
my($maxY,$NMlength,$Y1,$Y2,$Reg_00_r,$NMdepth_r) = @_;
my@startR = sort{$a<=>$b}keys%{$Reg_00_r};
my$line1a="";my$line1b="";
my$p = 0; my$i = 0;
for (my$n=0;$n<scalar@startR;$n++) {
	$line1a .= "par(new=TRUE)\nplot( c(";
	$p=$startR[$n]; $i=0;
	while ($i < scalar@{ ${$NMdepth_r}{$n} })
		{ $line1a .= $p.","; $p++; $i++; }
	chop $line1a;
	$line1a .= "), c(";
	$p=$startR[$n]; $i=0;
	while ($i < scalar@{ ${$NMdepth_r}{$n} }) {
		if (${$NMdepth_r}{$n}[$i] eq "") { 
			$line1a .= "0,";
# png: 			$line1b .= "lines (c($p,".($p+1)."), c(0,0), lwd=2, col=\"grey\")\n";
# pdf:
 				$line1b .= "lines (c($p,".($p+1)."), c(0,0), lwd=1.5, col=\"grey\")\n";
			}
		else { $line1a .= ${$NMdepth_r}{$n}[$i].","; }
		$p++; $i++;
		}
	chop $line1a;
# png:	$line1a .= "), xlim=c(-".(0.8*$NMlength/10).",".$NMlength."), ylim=c(-".($Y1*$maxY).",".($Y2*$maxY)."), type =\"l\", lwd=2, col=\"black\", axes=FALSE, ann=FALSE)\n";
# pdf:
		$line1a .= "), xlim=c(-".(0.8*$NMlength/10).",".$NMlength."), ylim=c(-".($Y1*$maxY).",".($Y2*$maxY)."), type =\"l\", lwd=1.5, col=\"black\", axes=FALSE, ann=FALSE)\n";
	}
return($line1a,$line1b);
}

###########################
#$line2: cmdR for cmdR for lines(intron)
#$line2 .= line2($maxY,$Y,\%{ ${$intron_r}{${$NMs_r}[$i]} });
#${$intron_r}{$startInt} = $endInt;
sub line2 {
my($maxY,$Y,$intron_r) = @_;
my$line="";
my@startIntron = sort{$a<=>$b}keys%{$intron_r}; 
for (my$ex=0;$ex<scalar@startIntron;$ex++) {
	$line .= "lines (c(".$startIntron[$ex].",".${$intron_r}{$startIntron[$ex]}."),c(-".($Y*$maxY).",-".($Y*$maxY)."), lwd=3)\n";
	}
return($line);
}

###########################
#cmdR for rect(bed) empty
#$line3a = line3a($maxY,$Y1,$Y2,\%NMbed00);	#empty black squares
#$line3b = line3a($maxY,$Y1,$Y2,\%Reg_00);	#full green squares
#${$NMbed00_r}{$startBed00} = $endBed00 (start of Region=0 )
sub line3a {
my($maxY,$Y1,$Y2,$NMbed00_r) = @_;
my@startBed = sort{$a<=>$b}keys%{$NMbed00_r};
my$line = "rect (c(";
foreach my$start (@startBed)
	{ $line .= $start.","; }
chop $line;
$line .= "), c(";
for (my$i=0;$i<scalar(@startBed);$i++)
	{ $line .= "-".($Y1*$maxY).","; }			
chop $line;
$line .= "), c(";
foreach my$start (@startBed)
	{ $line .= ${$NMbed00_r}{$start}.","; }
chop $line;
$line .= "), c(";
for (my$i=0;$i<scalar(@startBed);$i++)
	{ $line .= "-".($Y2*$maxY).","; }
chop $line;
$line .= "),";
return($line);
}

##########################
#$line3: cmdR for rect(not covered domains)
#$covEnd01{$file}{$NM}{$threshold}{$ex}{$notCovStart00} = $notCovEnd00;
#@{ $covStart01{$file}{$NM}{$threshold}{$ex} } = [ordered $notCovStarts00]
#my%line3 = line3($maxY,$Y1,$Y2,$Thresholds_r,$colors_r,scalar(keys%{ $Reg_00_r->{$gene} }),\%{ $covStart01_r->{$file} },\%{ $covEnd01_r->{$file} });
sub line3 {
my($maxY,$Y1,$Y2,$Thresholds_r,$colors_r,$N_Reg,$covStart01_r,$covEnd01_r) = @_;
my$col=0;
my%line;
foreach my$threshold (@{$Thresholds_r}) {
	my$ok=0;
	for (my$ex=0;$ex<$N_Reg;$ex++) {
		if ( scalar(keys(%{ ${$covEnd01_r}{$threshold}{$ex} })) != 0) 
			{ $ok = 1; last; } 
		}
	if ($ok) {
		$line{$threshold} = "rect (c(";
		for (my$ex=0;$ex<$N_Reg;$ex++) {			
			foreach my$pos (@{ ${$covStart01_r}{$threshold}{$ex} })
				{ $line{$threshold} .= $pos.","; }
			}
		chop $line{$threshold};
		$line{$threshold} .= "), c(";
		for (my$ex=0;$ex<$N_Reg;$ex++) {			
			for (my$i=0;$i<scalar@{ ${$covStart01_r}{$threshold}{$ex} };$i++)
				{ $line{$threshold} .= "-".($Y1*$maxY).","; }
			}			
		chop $line{$threshold};
		$line{$threshold} .= "), c(";
		for (my$ex=0;$ex<$N_Reg;$ex++) {
			foreach my$pos (@{ ${$covStart01_r}{$threshold}{$ex} })
				{ $line{$threshold} .= ${$covEnd01_r}{$threshold}{$ex}{$pos}.","; }
			}
		chop $line{$threshold};
		$line{$threshold} .= "), c(";
		for (my$ex=0;$ex<$N_Reg;$ex++) {
			for (my$i=0;$i<scalar@{ ${$covStart01_r}{$threshold}{$ex} };$i++)
				{ $line{$threshold} .= "-".($Y2*$maxY).","; }
			}
		chop $line{$threshold};
		$line{$threshold} .= "), lwd=0.01, col=\"".${$colors_r}[$col]."\", border =\"".${$colors_r}[$col]."\")\n";
		}
	$col++;
	}
return(%line);
}

#############################
#$line4a: cmdR for cmdR for rect(UTR exons)
#$line4a .= line4a($maxY,$Y1,$Y2,\%{ ${$UTR_r}{${$NMs_r}[$i]} });	
#${$UTR_r}{$startUTR} = $endUTR;
sub line4a {
my($maxY,$Y1,$Y2,$UTR_r) = @_;
my@startUTR = sort{$a<=>$b}keys%{$UTR_r};
my$line = "rect (c(";
for (my$ex=0;$ex<scalar@startUTR;$ex++)
	{ $line .= $startUTR[$ex].","; }
chop $line;
$line .= "), c(";
for (my$ex=0;$ex<scalar@startUTR;$ex++)
	{ $line .= "-".($Y1*$maxY).","; }
chop $line;
$line .= "), c(";
for (my$ex=0;$ex<scalar@startUTR;$ex++)
	{ $line .= ${$UTR_r}{$startUTR[$ex]}.","; }
chop $line;
$line .= "), c(";
for (my$ex=0;$ex<scalar@startUTR;$ex++)
	{ $line .= "-".($Y2*$maxY).","; }
chop $line;
#png: $line .= "), lwd=2, col=\"blue\")\n";
#pdf:
$line .= "), lwd=1.5, col=\"blue\")\n";
return($line);
}

##########################
#$line4b: cmdR for rect(coding exons)
#$line4b .= line4b($maxY,$Y1,$Y2,\%{ ${$Cod_r}{${$NMs_r}[$i]} });	
#${$Cod_r}{$startCod} = $endCod;
sub line4b {
my($maxY,$Y1,$Y2,$Cod_r) = @_;
my@startCod = sort{$a<=>$b}keys%{$Cod_r};
my$line = "rect (c(";
foreach (my$ex=0;$ex<scalar@startCod;$ex++) 
	{ $line .= $startCod[$ex].","; }
chop $line;
$line .= "), c(";
foreach (my$ex=0;$ex<scalar@startCod;$ex++) 
	{ $line .= "-".($Y1*$maxY).","; }
chop $line;
$line .= "), c(";
foreach (my$ex=0;$ex<scalar@startCod;$ex++) 
	{ $line .= ${$Cod_r}{$startCod[$ex]}.","; }
chop $line;
$line .= "), c(";
foreach (my$ex=0;$ex<scalar@startCod;$ex++) 
	{ $line .= "-".($Y2*$maxY).","; }
chop $line;
#png: $line .= "), lwd=2, col=\"blue\")\n";
#pdf:
$line .= "), lwd=1.5, col=\"blue\")\n";
return($line);
}

###########################
#$line4c: text Nb exon
#$line4c .= line4c($sens,$Rev,$maxY,$Y,\%{ ${$NM_Ex00_r}{${$NMs_r}[$i]} });
#${$NM_Ex00_r}{start of region}{start of exon} = end of exon (for start of region = 0)
sub line4c {
my($sens,$Rev,$maxY,$Y,$NM_Ex00_r) = @_;
my$line = "text( c(";
foreach my$startR (sort{$a<=>$b}keys%{$NM_Ex00_r}) {
	foreach my$startBed (sort{$a<=>$b}keys%{ ${$NM_Ex00_r}{$startR} })
		{ $line .= ((${$NM_Ex00_r}{$startR}{$startBed}+$startBed)/2).","; }
	}
chop $line;
$line .= "), c(";
my$nEx=0;
foreach my$startR (keys%{$NM_Ex00_r}) {
	for (my$i=0;$i<scalar(keys%{ ${$NM_Ex00_r}{$startR} });$i++)
		{ $line .= "-".($Y*$maxY).","; $nEx++; }
	}		
chop $line;
$line .= "), c(";
my$i=1;
if ($Rev) {
	foreach my$startR (keys%{$NM_Ex00_r}) {
		while ($i<=$nEx)
			{ $line .= "\"$i\" ,"; $i++; }
		}			
	}
else {
	if ($sens eq "-") {
		$i=$nEx;
		foreach my$startR (keys%{$NM_Ex00_r}) {
			while ($i>0)
				{ $line .= "\"$i\" ,"; $i--; }
			}			
		}
	}
chop $line;
## png: $line .= "), cex = 1)\n";
##pdf:
$line .= "), cex = 0.66)\n";
return($line);
}

##########################
#$line5: cmdR for legends
#my$line5 = line5a($comp,$colors_r,$Thresholds_r);
sub line5a {
my($comp,$colors_r,$Thresholds_r) = @_;
my$line = "legend(\"topright\", legend = c(\">=".${$Thresholds_r}[0]."\",";
foreach my$threshold (@{$Thresholds_r})
	{ $line .="\"$comp$threshold\","; }
chop $line;
$line .= "), col = c(\"green\",";
for (my$i=0;$i<scalar(@{$Thresholds_r});$i++)
	{ $line .="\"".${$colors_r}[$i]."\","; }
chop $line;
# png: $line .= "), pch = 15, bty = \"n\", pt.cex = 2.5, cex = 1, horiz = TRUE, inset = c(0, 0))\n";
# pdf:
$line .= "), pch = 15, bty = \"n\", pt.cex = 2, cex = 0.75, horiz = TRUE, inset = c(0, 0))\n";
return($line);
}

##########################
#$line5: cmdR for legends
sub line5b {
my($comp,$colors_r,$Thresholds_r) = @_;
my$line = "legend(\"topright\", legend = c(";
foreach my$threshold (@{$Thresholds_r})
	{ $line .="\"$comp$threshold\","; }
chop $line;
$line .= "), col = c(";
for (my$i=0;$i<scalar(@{$Thresholds_r});$i++)
	{ $line .="\"".${$colors_r}[$i]."\","; }
chop $line;
# png: $line .= "), pch = 15, bty = \"n\", pt.cex = 2.5, cex = 1, horiz = TRUE, inset = c(0, 0))\n";
# pdf:
$line .= "), pch = 15, bty = \"n\", pt.cex = 2, cex = 0.75, horiz = TRUE, inset = c(0, 0))\n";
return($line);
}

#############################
#$line6: cmdR for mutations (symbol)
#$line6a = line6a($maxY,$Y1,$Y2,$NMmut00_r);
#$NMmut00{$NM}{$mut00} = $infoMut (start of Region=0 )
sub line6a {
my($maxY,$Y1,$Y2,$NMmut00_r) = @_;
my$line="";
#above
$line .= "points(c(";
foreach my$mut (keys%{$NMmut00_r})
	{ $line .= $mut.","; }
chop $line;
$line .= "), c(";
for (my$i=0;$i<scalar(keys%{$NMmut00_r});$i++)
	{ $line .= "-".($Y1*$maxY).","; }			
chop $line;
#png: $line .= "), pch = 6, col = \"black\", lwd = 2, cex = 1)\n";
#pdf:
$line .= "), pch = 6, col = \"black\", lwd = 1, cex = 0.75)\n";
#below
$line .= "points(c(";
foreach my$mut (keys%{$NMmut00_r})
	{ $line .= $mut.","; }
chop $line;
$line .= "), c(";
for (my$i=0;$i<scalar(keys%{$NMmut00_r});$i++)
	{ $line .= "-".($Y2*$maxY).","; }			
chop $line;
#png: $line .= "), pch = 2, col = \"black\", lwd = 2, cex = 1)\n";
#pdf:
$line .= "), pch = 2, col = \"black\", lwd = 1, cex = 0.75)\n";
return($line);
}

##############################
#$line6b: cmdR for mutations (text)
#$line6b = line6b($maxY,$Y,$NMmut00_r);
#$NMmut00{$NM}{$mut00} = $infoMut (start of Region=0 )
sub line6b {
my($maxY,$Y,$NMmut00_r) = @_;
my$line="";
$line = "text( c(";
foreach my$mut (keys%{$NMmut00_r})
	{ $line .= $mut.","; }
chop $line;
$line .= "), c(";
for (my$i=0;$i<scalar(keys%{$NMmut00_r});$i++)
	{ $line .= "-".($Y*$maxY).","; }			
chop $line;
$line .= "), c(";
foreach my$mut (keys%{$NMmut00_r})
	{ $line .= "\"".${$NMmut00_r}{$mut}."\" ,"; }			
chop $line;
#png: $line .= "), adj = c(0.5,1.75), cex = 1.25, srt = -15)\n";
#pdf:
$line .= "), adj = c(0.5,1.75), cex = 0.66, srt = -15)\n";
#$line .= "), pos = 1, cex = 1)\n";
#$line .= "), cex = 1)\n";
return($line);
}

##############################
#$line7: cmdR for text NM
#my$line7 = line7($length,$maxY,$NMs_r);
sub line7 {
my($NMlength,$maxY,$NMs_r) = @_;
my$line="";
$line = "text( c(";
for (my$i=0;$i<scalar@{$NMs_r};$i++)
	{ $line .= "-".(0.6*$NMlength/10).","; }
chop $line;
$line .= "), c(";
for (my$i=0;$i<scalar@{$NMs_r};$i++)
	{ $line .= "-".((0.425+0.25*$i)*$maxY).","; }		
chop $line;
$line .= "), c(";
for (my$i=0;$i<scalar@{$NMs_r};$i++)
	{ $line .= "\"".${$NMs_r}[$i]."\" ,"; }			
chop $line;
#$line .= "), adj = c(0,2), cex = 1)\n";
## png: $line .= "), cex = 1.2)\n";
##pdf:
$line .= "), cex = 0.75)\n";
return($line);
}


##############################
#graph with all transcripts of a gene
# graphSumG($nGraf,$gene,$suff,$outdir,$sens,$chr,$NMlength{$gene},\@NMs,\@Files,\@colors,\@Thresholds,$Regions{"raw"},$NM_Ex{"rev"},$RegBed{"rev"},$intron{"rev"},$UTR{"rev"},$Cod{"rev"},$RegMut{"rev"},$NMCovEnd{$gene},$NMCovVal{$gene});

sub graphSumG {

my($nGraf,$gene,$suff,$outdir,$sens,$Rev,$chr,$length,$NMs_r,$Files_r,$colors_r,$Thresholds_r,$Regions_r,$NM_Ex00_r,$NMbed00_r,$intron_r,$UTR_r,$Cod_r,$NMmut00_r,$NMCovEnd_r,$NMCovVal_r) = @_;

print "print CMDR graphSumG: $gene\n";
open (CMDR, ">$outdir/$gene\_temp.R") || die;
print CMDR "#!/usr/bin/env Rscript\n\n" ;
my$nNMs=scalar(@{$NMs_r});
my$maxY=scalar(@{$Files_r});
my@Starts = sort{$a<=>$b}(keys%{ $Regions_r->{$chr}{$gene} });
my$startReg = $Starts[0];
my$endReg = $Regions_r->{$chr}{$gene}{$Starts[-1]};

print CMDR "pdf(\"".$outdir."/cov_All/".$gene."_covSum.pdf\", width=11.69, height=".(4.135*(0.75+0.875+(0.125*$nNMs))/1.75).")\n";

if ($sens eq "+")
	{ print CMDR 
#png: "plot (c(0,0), xlim=c(-".(0.8*$length/10).",".$length."), ylim=c(-".((0.5+0.25*$nNMs)*$maxY).",".(1.25*$maxY)."), type =\"l\", xlab = \"chr $chr\t$startReg >>> $endReg\", ylab = \"nber of covered samples\", cex.lab=1.5, cex.axis=1.2, yaxt='n')
#axis(2, at=seq(0, $maxY, length = 2),labels=seq(0, $maxY, length = 2), las=2)\n";
#pdf:
"plot (c(0,0), xlim=c(-".(0.8*$length/10).",".$length."), ylim=c(-".((0.5+0.25*$nNMs)*$maxY).",".(1.25*$maxY)."), type =\"l\", xlab = \"chr $chr: $startReg >>> $endReg\", ylab = \"nber of covered samples\", cex.lab=1, cex.axis=0.8, yaxt='n')
axis(2, at=seq(0, $maxY, length = 2),labels=seq(0, $maxY, length = 2), las=2, cex.axis=0.8)\n";
	}
else
	{ print CMDR 
#png: "plot (c(0,0), xlim=c(-".(0.8*$length/10).",".$length."), ylim=c(-".((0.5+0.25*$nNMs)*$maxY).",".(1.25*$maxY)."), type =\"l\", xlab = \"chr $chr\t$startReg <<< $endReg\", ylab = \"nber of covered samples\", cex.lab=1.5, cex.axis=1.2, yaxt='n')
#axis(2, at=seq(0, $maxY, length = 2),labels=seq(0, $maxY, length = 2), las=2)\n";
#pdf:
"plot (c(0,0), xlim=c(-".(0.8*$length/10).",".$length."), ylim=c(-".((0.5+0.25*$nNMs)*$maxY).",".(1.25*$maxY)."), type =\"l\", xlab = \"chr $chr: $startReg <<< $endReg\", ylab = \"nber of covered samples\", cex.lab=1, cex.axis=0.8, yaxt='n')
axis(2, at=seq(0, $maxY, length = 2),labels=seq(0, $maxY, length = 2), las=2, cex.axis=0.8)\n";
	}

my($Y,$Y1,$Y2);
#$line1
#$NMCovEnd{$threshold}{$r}{$start}=$end
#$NMCovVal{$threshold}{$r}{$start}=$value
my$col=0;
my%line1;
#$line1: cmdR for plot cov line
foreach my$threshold (@{$Thresholds_r}) {
	$line1{$threshold} = line1g2(${$colors_r}[$col],\%{ ${$NMCovEnd_r}{$threshold} },\%{ ${$NMCovVal_r}{$threshold} });
	$col++;
	}

#$line2: cmdR for cmdR for lines(intron)	#$intron{$NM}{$startInt} = $endInt;
my$line2="";
for(my$i=0;$i<$nNMs;$i++) { 
	$Y=(0.425+0.25*$i);
	$line2 .= line2($maxY,$Y,\%{ ${$intron_r}{${$NMs_r}[$i]} });
	}

#$line3a: cmdR for rect(bed, empty)	#$NMbed00{$NM}{$startBed00} = $endBed00 (start of Region=0 )
$Y1=0.1;$Y2=0.25;
my$line3a = line3a($maxY,$Y1,$Y2,\%{ $NMbed00_r->{$gene} });			#\%{ ${$NMbed00_r}{$gene} }	${$NMbed00_r}{$gene}

#$line4a: cmdR for cmdR for rect(UTR exons)	#$UTR{$NM}{$startUTR} = $endUTR;
my$line4a="";
for(my$i=0;$i<$nNMs;$i++) {
	if (scalar(keys%{ ${$UTR_r}{${$NMs_r}[$i]} }) != 0) {
		$Y1=(0.375+0.25*$i);$Y2=(0.475+0.25*$i);
		$line4a .= line4a($maxY,$Y1,$Y2,\%{ ${$UTR_r}{${$NMs_r}[$i]} });
		}
	}

#$line4b: cmdR for rect(coding exons)	#$Cod{$NM}{$startCod} = $endCod;
my$line4b="";
for(my$i=0;$i<$nNMs;$i++) {
	if (scalar(keys%{ ${$Cod_r}{${$NMs_r}[$i]} }) != 0) {
		$Y1=(0.35+0.25*$i);$Y2=(0.5+0.25*$i);
		$line4b .= line4b($maxY,$Y1,$Y2,\%{ ${$Cod_r}{${$NMs_r}[$i]} });
		}
	}

#$line4c: Nb exon
#%NM_Ex00{NM}{start of region}{start of exon} = end of exon (for start of region = 0)
my$line4c="";
for(my$i=0;$i<$nNMs;$i++) {
	$Y=(0.3+0.25*$i);
	$line4c .= line4c($sens,$Rev,$maxY,$Y,\%{ ${$NM_Ex00_r}{${$NMs_r}[$i]} });
	}

#$line5: cmdR for legends
my$comp = ">=";
my$line5 = line5b($comp,$colors_r,$Thresholds_r);

#$line6: cmdR for mutations (symbol)
#$NMmut00{$NM}{$mut00} = $infoMut (start of Region=0 )
my$line6a = ""; my$line6b = "";
if (scalar(keys%{ $NMmut00_r->{$gene} }) != 0) {
	$Y1=0.05;$Y2=(0.3+0.25*$nNMs);
	$line6a = line6a($maxY,$Y1,$Y2,\%{ $NMmut00_r->{$gene} });
	#$line6b: cmdR for mutations (text)
	#$NMmut00{$NM}{$mut00} = $infoMut (start of Region=0 )
	$Y=(0.3+0.25*$nNMs);
	$line6b = line6b($maxY,$Y,\%{ $NMmut00_r->{$gene} });
	}

#$line7: NM names
my$line7 = line7($length,$maxY,$NMs_r);

print CMDR "abline(h=0, col=\"black\", lty = \"dotted\", lwd=1)\n"; 
print CMDR "abline(h=$maxY, col=\"darkgrey\", lty = \"dotted\", lwd=1)\n"; 

foreach my$threshold (reverse@{$Thresholds_r})
	{ print CMDR $line1{$threshold}; }

print CMDR $line2.$line4a.$line4b.$line4c.$line5.$line6a.$line6b.$line7;

print CMDR $line3a." lwd=1, col=NA, border =\"black\")\n";


#png: print CMDR "mtext(\"$gene\", side = 3, outer=TRUE, line=-3, cex=1.5)\n";
#pdf:
print CMDR "mtext(\"$gene\", side = 3, outer=TRUE, line=-3, cex=1)\n";
print CMDR "dev.off();\n";
close CMDR;
system "Rscript $outdir/$gene\_temp.R";
unlink "$outdir/$gene\_temp.R";

}


##############################
#graph foreach transcript of a gene
#graphSumN($nGraf,$gene,$suff,$outdir,$sens,$chr,\%NMlength,\@NMs,\@Files,\@colors,\@Thresholds,$Regions{"raw"},$NM_Ex{"rev"},$RegBed{"rev"},$intron{"rev"},$UTR{"rev"},$Cod{"rev"},$RegMut{"rev"},\%NMCovEnd,\%NMCovVal);

sub graphSumN {	#still in png, but not used anymore

my($nGraf,$gene,$suff,$outdir,$sens,$Rev,$chr,$lengths_r,$NMs_r,$Files_r,$colors_r,$Thresholds_r,$Regions_r,$NM_Ex00_r,$NMbed00_r,$intron_r,$UTR_r,$Cod_r,$NMmut00_r,$NMCovEnd_r,$NMCovVal_r)= @_;

print "print CMDR graphSumN: $gene\n";

my$nNM=scalar@{$NMs_r};
my$maxY=scalar(@{$Files_r});
if ($nGraf eq "max") { $nGraf = $nNM; }
my$cmdR = "";
my$i=1;my$n=1;my$N=1;
foreach my$NM (sort(@{$NMs_r})) {
	my@Starts = sort{$a<=>$b}(keys%{ ${$Regions_r}{$chr}{$NM} });
	my$startReg = $Starts[0];
	my$endReg = ${$Regions_r}{$chr}{$NM}{$Starts[-1]};

	if ($sens eq "+") {
		$cmdR .= 
"plot (c(0,0), xlim=c(0,".${$lengths_r}{$NM}."), ylim=c(-".(0.6*$maxY).",".(1.25*$maxY)."), type =\"l\", xlab = \"transcript : $NM\t\tchr $chr\t$startReg >>> $endReg\", ylab = \"nber of covered samples\", cex.lab=1.5, cex.axis=1.2, yaxt='n')
axis(2, at=seq(0, $maxY, length = 2), labels=seq(0, $maxY, length = 2), las=2)\n";
		}
	else
		{ $cmdR .= 
"plot (c(0,0), xlim=c(0,".${$lengths_r}{$NM}."), ylim=c(-".(0.6*$maxY).",".(1.25*$maxY)."), type =\"l\", xlab = \"transcript : $NM\t\tchr $chr\t$startReg <<< $endReg\", ylab = \"nber of covered samples\", cex.lab=1.5, cex.axis=1.2, yaxt='n')
axis(2, at=seq(0, $maxY, length = 2), labels=seq(0, $maxY, length = 2), las=2)\n";
		}

	#$line1: cmdR for plot cov line
	#$NMCovEnd{$threshold}{$r}{$start}=$end
	#$NMCovVal{$threshold}{$r}{$start}=$value
	my$col=0;my($Y,$Y1,$Y2);
	my%line1;
	foreach my$threshold (@{$Thresholds_r}) {
		$line1{$threshold} = line1g2(${$colors_r}[$col],\%{ ${$NMCovEnd_r}{$NM}{$threshold} },\%{ ${$NMCovVal_r}{$NM}{$threshold} });
		$col++;
		}

	#$line2: cmdR for cmdR for lines(intron)	#$intron{$NM}{$startInt} = $endInt;
	$Y=0.425;
	my$line2 = line2($maxY,$Y,\%{ ${$intron_r}{$NM} });

	#$line3a: cmdR for rect(bed, empty)	#$NMbed00{$NM}{$startBed00} = $endBed00 (start of Region=0 )
	$Y1=0.1;$Y2=0.25;
	my$line3a = line3a($maxY,$Y1,$Y2,\%{ $NMbed00_r->{$NM} });

	#$line4a: cmdR for cmdR for rect(UTR exons)	#$UTR{$NM}{$startUTR} = $endUTR;
	$Y1=0.375;$Y2=0.475;
	my$line4a = "";
	if (scalar(keys%{ ${$UTR_r}{$NM} }) != 0 ) { $line4a .= line4a($maxY,$Y1,$Y2,\%{ ${$UTR_r}{$NM} }); }

	#$line4b: cmdR for rect(coding exons)	#$Cod{$NM}{$startCod} = $endCod;
	$Y1=0.35;$Y2=0.5;
	my$line4b = "";
	if (scalar(keys%{ ${$Cod_r}{$NM} }) != 0) { $line4b .= line4b($maxY,$Y1,$Y2,\%{ ${$Cod_r}{$NM} }); }

	#$line4c: Nb exon
	#%NM_Ex00{NM}{start of region}{start of exon} = end of exon (for start of region = 0)
	$Y=0.3;
	my$line4c = line4c($sens,$Rev,$maxY,$Y,\%{ ${$NM_Ex00_r}{$NM} });

	#$line5: cmdR for legends
	my$comp = ">=";
	my$line5 = line5b($comp,$colors_r,$Thresholds_r);
	
	#$line6: cmdR for mutations (symbol)
	#$NMmut00{$NM}{$mut00} = $infoMut (start of Region=0 )
	my$line6a = ""; my$line6b = "";
	if (scalar(keys%{ ${$NMmut00_r}{$NM} }) != 0) {
		$Y1=0.05;$Y2=0.55;
		$line6a = line6a($maxY,$Y1,$Y2,\%{ ${$NMmut00_r}{$NM} });
		#$line6b: cmdR for mutations (text)
		#$NMmut00{$NM}{$mut00} = $infoMut (start of Region=0 )
		$Y=0.5;
		$line6b = line6b($maxY,$Y,\%{ ${$NMmut00_r}{$NM} });
		}
	
	$cmdR .= "abline(h=0, col=\"black\", lty = \"dotted\", lwd=1)\n";
	$cmdR .= "abline(h=$maxY, col=\"darkgrey\", lty = \"dotted\", lwd=1)\n"; 

	foreach my$threshold (reverse@{$Thresholds_r})
		{ $cmdR .= $line1{$threshold}; }

	$cmdR .= $line2.$line4a.$line4b.$line4c.$line5.$line6a.$line6b;

	$cmdR .= $line3a." lwd=1, col=NA, border =\"black\")\n";
	

	if ($i==($nNM) || $n==$nGraf) {

		open (CMDR, ">$outdir/$gene\_temp.R") || die;
		print CMDR "#!/usr/bin/env Rscript\n\n" ;
		if ($nGraf==$nNM) { print CMDR
"png(\"".$outdir."/cov_All/".$gene."_covSum.png\", 1500, ".($nGraf*400).")\n
par(mfrow=c($nGraf,1))\n"; }
		else {
			if ($N>1) { print CMDR
"png(\"".$outdir."/cov_All/".$gene."_covSum_$N.png\", 1500, ".($nGraf*400).")\n
par(mfrow=c($nGraf,1))\n"; }
			else { print CMDR 
"png(\"".$outdir."/cov_All/".$gene."_covSum_$N.png\", 1500, ".($n*400).")\n
par(mfrow=c($n,1))\n"; }
			}
		print CMDR "$cmdR";
		print CMDR "mtext(\"$gene\", side = 3, outer=TRUE, line=-3, cex=1.5)\n";
		print CMDR "dev.off();\n";
		close CMDR;
		system "Rscript $outdir/$gene\_temp.R";
		unlink "$outdir/$gene\_temp.R";
		$cmdR="";
		$n=0;
		$N++;
		}

	$i++;$n++;

	}
}


###############################
#graph line (disabled)
sub line1g {

my($NM,$maxY,$Y1,$Y2,$NMlength,$color,$h1,$h2) = @_;
my%Reg_00=%$h1;
my%NMCov=%$h2;
my@startR = sort{$a<=>$b}keys%Reg_00;
my$line="";
for (my$n=0;$n<scalar@startR;$n++)
	{
	$line .= "par(new=TRUE)\nplot( c(";
	my$p=$startR[$n]; my$i=0;
	while ($i<scalar@{ $NMCov{$n} })
		{ $line .= $p.","; $p++; $i++; }
	chop $line;
	$line .= "), c(";
	$i=0;
	while ($i<scalar@{ $NMCov{$n} })
		{ $line .= $NMCov{$n}[$i].","; $i++; }
	chop $line;
	$line .= "), xlim=c(0,".$NMlength."), ylim=c(-".($Y1*$maxY).",".($Y2*$maxY)."), type =\"l\", lwd=2, col=\"$color\", axes=FALSE, ann=FALSE)\n";
	}
return($line);

}
####################################
#same but graph rect
sub line1g2 {
my($color,$NMCovEnd_r,$NMCovVal_r) = @_;

#$NMCovEnd{$NM}{$ex}{$start}=$end
#$NMCovVal{$NM}{$ex}{$start}=$value
my@startR = sort{$a<=>$b}keys%{$NMCovEnd_r};
my$line = "";
for (my$r=0;$r<scalar@startR;$r++) {
	my@Starts= sort{$a<=>$b}keys%{ ${$NMCovEnd_r}{$startR[$r]} };
	$line .= "rect (c(";
	for (my$i=0;$i<scalar@Starts;$i++)
		{ $line .= $Starts[$i].","; }
	chop $line;
	$line .= "), c(";
	for (my$i=0;$i<scalar@Starts;$i++)
		{ $line .= "0,"; }
	chop $line;
	$line .= "), c(";
	for (my$i=0;$i<scalar@Starts;$i++)
		{ $line .= ${$NMCovEnd_r}{$startR[$r]}{$Starts[$i]}.","; }
	chop $line;
	$line .= "), c(";
	for (my$i=0;$i<scalar@Starts;$i++)
		{ $line .= ${$NMCovVal_r}{$startR[$r]}{$Starts[$i]}.","; }
	chop $line;
	$line .= "), border=NA, col=\"$color\")\n";
	}
return($line);

}


########################
#for $All $NM on same graph
# CNVgraphByGene($file,$suff,"$outdir/CNV_analysis/$sName2{$file}",$gene,$sens,$maxGr,$Rev,$chr,$CNV_maxCov{$gene},$pThreshold,$NMlength{$gene},\@Files,\%sName2,$Regions{"raw"},$Regions{"rev"},\@NMs,$NM_Ex{"rev"},$RegBed{"rev"},$intron{"rev"},$UTR{"rev"},$Cod{"rev"},$CNV_NMdepth{$gene});

sub CNVgraphByGene {

my($file,$suff,$outdir,$gene,$sens,$maxD,$Rev,$chr,$maxY,$threshold,$NMlength,$Files_r,$sampleName_r,$Regions_r,$Reg_00_r,$NMs_r,$NM_Ex00_r,$NMbed00_r,$intron_r,$UTR_r,$Cod_r,$NMdepth_r) = @_;

print "print CMDR CNVgraphByGene: ".${$sampleName_r}{$file}." , in $gene\n";

unless ($maxD) {
	if ($maxY<10) { $maxY = 10; }
	}

my$nNMs = scalar@{$NMs_r};
my@Starts = sort{$a<=>$b}(keys%{ $Regions_r->{$chr}{$gene} });
my$startReg = $Starts[0];
my$endReg = $Regions_r->{$chr}{$gene}{$Starts[-1]};

my($Y,$Y1,$Y2);
#$line1: cmdR for plot cov line
my$line1="";
foreach my$f2 (@{$Files_r}) {
	unless ($f2 eq $file) {
		#@{ $NMdepth{$file}{$ex} } = [ cov foreach ordered bp of exon ]
		$line1 .= line1c("black",\%{ $Reg_00_r->{$gene} },\%{ $NMdepth_r->{$f2} });
		}
	}
$line1 .= line1c("green",\%{ $Reg_00_r->{$gene} },\%{ $NMdepth_r->{$file} });

#$line2: cmdR for cmdR for lines(intron)	#$intron{$startInt} = $endInt;
my$line2="";
for(my$i=0;$i<$nNMs;$i++) { 
	$Y=(0.425+0.25*$i);
	$line2 .= line2($maxY,$Y,\%{ $intron_r->{${$NMs_r}[$i]} });
	}

#$line3a: cmdR for rect(bed, empty)	#$NMbed00{$startBed00} = $endBed00 (start of Region=0 )
$Y1=0.1;$Y2=0.25;
my$line3a = line3a($maxY,$Y1,$Y2,\%{ $NMbed00_r->{$gene} });


#$line4a: cmdR for cmdR for rect(UTR exons)	#$UTR{$NM}{$startUTR} = $endUTR;
my$line4a="";
for(my$i=0;$i<$nNMs;$i++) {
	if (scalar(keys%{ $UTR_r->{${$NMs_r}[$i]} })!=0) {
		$Y1=(0.375+0.25*$i);$Y2=(0.475+0.25*$i);
		$line4a .= line4a($maxY,$Y1,$Y2,\%{ $UTR_r->{${$NMs_r}[$i]} });
		}
	}

#$line4b: cmdR for rect(coding exons)	#$Cod{$NM}{$startCod} = $endCod;
my$line4b="";
for(my$i=0;$i<$nNMs;$i++) {
	if (scalar(keys%{ $Cod_r->{${$NMs_r}[$i]} })) {
		$Y1=(0.35+0.25*$i);$Y2=(0.5+0.25*$i);
		$line4b .= line4b($maxY,$Y1,$Y2,\%{ $Cod_r->{${$NMs_r}[$i]} });
		}
	}

#$line4c: Nb exon
#%NM_Ex00{NM}{start of region}{start of exon} = end of exon (for start of region = 0)
my$line4c="";
for(my$i=0;$i<$nNMs;$i++) {
	$Y=(0.3+0.25*$i);
	$line4c .= line4c($sens,$Rev,$maxY,$Y,\%{ ${$NM_Ex00_r}{${$NMs_r}[$i]} });
	}

#$line7: NM names
my$line7 = line7($NMlength,$maxY,$NMs_r);

#print CMDR:
open (CMDR, ">$outdir/".${$sampleName_r}{$file}.".$gene.temp.R") || die "cannot create $outdir/".${$sampleName_r}{$file}.".$gene.temp.R";
print CMDR "#!/usr/bin/env Rscript\n\n" ;
print CMDR "pdf(\"$outdir/".${$sampleName_r}{$file}.".$gene.pdf\", width=11.69, height=".(4.135*(0.75+0.875+(0.125*$nNMs))/1.75).")\n";
 
#png: my$cmdR = "plot (c(0,0), xlim=c(-".(0.8*$NMlength/10).",".$NMlength."), ylim=c(-".((0.5+0.25*$nNMs)*$maxY).",".(1.25*$maxY)."), type =\"n\", xlab = \"sample : ".${$sampleName_r}{$file}."\", ylab = \"normalized_depth\", cex.lab=1.5, cex.axis=1.2)\n";
#pdf:
my$cmdR = "plot (c(0,0), xlim=c(-".(0.8*$NMlength/10).",".$NMlength."), ylim=c(-".((0.5+0.25*$nNMs)*$maxY).",".(1.25*$maxY)."), type =\"n\", xlab = \"sample : ".${$sampleName_r}{$file}."\", ylab = \"normalized_depth\", cex.lab=1, cex.axis=0.8)\n";

$cmdR .= $line1;

$cmdR .= "abline(h=0, col=\"darkgrey\", lty = \"dotted\", lwd=1.5)\n";
if ($threshold) { $cmdR .= "abline(h=$threshold, col=\"darkgrey\", lty = \"dotted\", lwd=1.5)\n"; }

$cmdR .= $line2.$line4a.$line4b.$line4c.$line7;

$cmdR .= $line3a." lwd=1, col=NA, border =\"black\")\n";

print CMDR "$cmdR";
if ($sens eq "+") {
	print CMDR "mtext(\"$gene : chr $chr:$startReg >>> $endReg\", side = 3, outer=TRUE, line=-3, cex=1)\n";
	}
else {
	print CMDR "mtext(\"$gene : chr $chr:$startReg <<< $endReg\", side = 3, outer=TRUE, line=-3, cex=1)\n";
	}
print CMDR "dev.off();\n";
close CMDR;
system "Rscript $outdir/".${$sampleName_r}{$file}.".$gene.temp.R";
unlink "$outdir/".${$sampleName_r}{$file}.".$gene.temp.R";

}


########################
#$line1 .= line1c("green",\%{ $Reg_00_r->{$gene} },\%{ $NMdepth_r->{$file} });
sub line1c {
my($col,$Reg_00_r,$NMdepth_r) = @_;
my@startR = sort{$a<=>$b}keys%{$Reg_00_r};
my$line1="";
my$p=0;my$i=0;
for (my$n=0;$n<scalar@startR;$n++) {
	$line1 .= "lines( c(";
	$p=$startR[$n]; $i=0;
	while ($i<scalar@{ ${$NMdepth_r}{$n} })
		{ $line1 .= $p.","; $p++; $i++; }
	chop $line1;
	$line1 .= "), c(";
	$p=$startR[$n]; $i=0;
	while ($i<scalar@{ ${$NMdepth_r}{$n} }) {
		if (${$NMdepth_r}{$n}[$i] eq "") { $line1 .= "0,"; }
		else { $line1 .= ${$NMdepth_r}{$n}[$i].","; }
		$p++; $i++;
		}
	chop $line1;
	$line1 .= "), type =\"l\", lwd=1.5, col=\"$col\")\n";
	}
return($line1);
}


############################



1;


