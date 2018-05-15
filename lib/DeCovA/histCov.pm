package DeCovA::histCov;


use strict;
use warnings;



############################
#	if (!$bedFile || $jobs{"genePlot"}) {
#		covPlot($bin,$maxPl,$outPlotDir",all_Genes_Regions",\@Files,\%sName2,$smplIdx,\%{ ${$depthFilePerChr}{"raw"} },$coordBed_r,$histDepth,$totLength);
#		}
#	if ($bedFile) {
#		covPlot($bin,$maxPl,$outPlotDir,$bedName,\@Files,\%fName,$smplIdx,\%{ ${$depthFilePerChr}{"raw"} },$Bed_r,$histDepth,$totLength);
#		}

sub covPlot {

my($bin,$maxPl,$outdir,$outName,$Files,$smplName,$smplIdx,$depthFilePerChr,$allInterval_r,$threshold,$lengthBed,$depthCount,$covBases,$totBases)=@_;

print "\n\n####\n\ndoing depthPlot for $outName (by samples and avg of all samples)\n\n####\n";

## ${$depthCount}{$file}{depth value} = nber pos with this depth
unless ($depthCount) {
	($lengthBed,$depthCount,$covBases,$totBases) = intersectForHist($Files,$smplIdx,$depthFilePerChr,$allInterval_r,$threshold);
	}

print "\n\t## printing per bin depth count and cumulative depth count tab files\n";

open(my$fhBin, ">", "$outdir/$outName\_binCovPlot.txt") || die "can't create $outdir/$outName\_binCovPlot.txt ($!)\n";
open(my$fhSum, ">", "$outdir/$outName\_sumCovPlot.txt") || die "can't create $outdir/$outName\_sumCovPlot.txt ($!)\n";
print $fhBin "depth\t";
print $fhSum "depth\t";
my$i=0;
while ( ($i*$bin) <= $maxPl ) {
	print $fhBin ($bin*$i)."-".(($bin*($i+1))-1)."\t";
	print $fhSum ">".($bin*$i)."\t";
	$i++;
	}
print $fhBin ">=$maxPl\n";
print $fhSum "\n";

foreach my$file (@{$Files}) {
	print "\t\tanalysing ".${$smplName}{$file}." for plot on $outName\n";

	my%histD;
	foreach my$d (keys%{ ${$depthCount}{$file} }) {
		if ($maxPl && $d > $maxPl) { $histD{$file}{($maxPl+1)} += ${$depthCount}{$file}{$d}; }
		else { $histD{$file}{$d} = ${$depthCount}{$file}{$d} / $lengthBed; }
		}
	if ($maxPl) {
		if (exists $histD{$file}{($maxPl+1)}) { $histD{$file}{($maxPl+1)} /= $lengthBed; }
		else { $histD{$file}{($maxPl+1)} = 0; }
		}

	print $fhBin "sample_".${$smplName}{$file}."\t";
	my$histCount;
	$i=0;
	while ( ($i*$bin) <= $maxPl ) {
		$histCount=0;
		for (my$j=($bin*$i);$j<($bin*($i+1));$j++) {
			if (exists $histD{$file}{$j}) { $histCount += $histD{$file}{$j}; }
			}
		print $fhBin $histCount."\t";
		$i++;
		}
	my@array=sort{$a<=>$b}keys%{ $histD{$file} };
	$histCount=0;
	for ($i=$maxPl;$i<=$#array;$i++) {
		$histCount += $histD{$file}{$array[$i]};
		}
	print $fhBin $histCount."\n";

	print $fhSum "sample_".${$smplName}{$file}."\t";
	my$totCount=0;
	if (exists $histD{$file}{0}) { $totCount += $histD{$file}{0}; }
	print $fhSum (1-$totCount)."\t";
	$i=0;
	while ( ($i*$bin) < $maxPl ) {
		$histCount=0;
		for (my$j=(($bin*$i)+1);$j<=($bin*($i+1));$j++) {
			if (exists $histD{$file}{$j}) { $histCount += $histD{$file}{$j}; }
			}
		$totCount += $histCount;
		print $fhSum (1-$totCount)."\t";
		$i++;
		}
	print $fhSum "\n";
	}

close($fhBin); 
close($fhSum);

#R barplot -sum, mean of all samples
if (scalar@{$Files} > 1) {
	meanCovPlot($bin,$maxPl,$outdir,$outName,$lengthBed);
	}
#R barplot -sum, 1graph / sample
allCovPlot($bin,$maxPl,$outdir,$outName,$lengthBed,$Files,$smplName);

depthBoxPlot($outdir,$outName,$Files,$smplName,$threshold,$lengthBed,$depthCount,$covBases,$totBases);


}


##########
## ($lengthBed,$hist_r) = intersectForHist($Files,$smplIdx,$depthFilePerChr,$maxPl,$allInterval_r);

sub intersectForHist {
print "\n\t##intersecting plot target regions with depth Files\n";
my($Files,$smplIdx,$depthFilePerChr,$Intervals_r,$threshold) = @_;

my$lengthBed = 0;
my(%allDepth,%covBases,%totBases);		#for hist depth : $allDepth{$file}{depth value} = nber pos with this depth

foreach my$chr (keys%{$Intervals_r}) {
	my@Starts = sort{$a<=>$b}(keys%{ ${$Intervals_r}{$chr} });
	my$c=0;
	open(my$fh, "<", ${$depthFilePerChr}{$chr}) || die "can't read file ".${$depthFilePerChr}{$chr}." ($!)\n";
	while (my$line = <$fh>) {
		chomp $line;
		my@tab = split(/\t/,$line);
		my$pos = $tab[0];
		while ( ($pos > ${$Intervals_r}{$chr}{$Starts[$c]}) && ($c < $#Starts) ) { $c++; }
		if ( ($pos >= $Starts[$c]) && ($pos <= ${$Intervals_r}{$chr}{$Starts[$c]}) ) {
			$lengthBed++;
			foreach my$f (@{$Files}) {
				my$depth = $tab[${$smplIdx}{$f}];
				$totBases{$f} += $depth;
				$allDepth{$f}{$depth}++;
				if ($threshold && $depth >= $threshold) { $covBases{$f}++; }
				}
			}
		if ( $pos > ${$Intervals_r}{$chr}{$Starts[-1]} ) { last; }
		}
	}

return($lengthBed,\%allDepth,\%covBases,\%totBases);

}

##########

sub meanCovPlot {
my($bin,$maxPl,$outdir,$outName,$lengthBed) = @_;
print "\n\t## doing meanCovPlot\n";
my$line="";
my$i=0;
while ( ($i*$bin) <= $maxPl ) {
	$line .= "\">".($bin*$i)."\",";
	$i++;
	}
chop $line;
open (CMDR, ">$outdir/plot_temp.R") || die;
print CMDR "#!/usr/bin/env Rscript\n\n" ;
print CMDR "pdf(\"$outdir/$outName\_meanCovPlot.pdf\", width=11.69, height=4.135)
mat <- read.table(\"$outdir/$outName\_sumCovPlot.txt\",header=T, row.names=1)
allMean = apply(mat,2,mean)
std = apply(mat,2,sd)
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop(\"vectors must be same length\")
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
barx <- barplot(allMean, names.arg=c(";
print CMDR $line;
print CMDR "), ylim=c(0,1), col=\"blue\", axis.lty=1, main=\"mean of all samples\", xlab=\"Depth\", ylab=\"fraction of $outName ($lengthBed bp)\")
error.bar(barx,allMean, 1.96*std/10)
dev.off();warnings();\n";
close CMDR;
system "Rscript $outdir/plot_temp.R";
unlink "$outdir/plot_temp.R";
}


##########
sub allCovPlot {
my($bin,$maxPl,$outdir,$outName,$lengthBed,$Files,$smplName) = @_;
print "\n\t## doing allCovPlot\n";
my$line="";
my$i=0;
while ( ($i*$bin) <= $maxPl ) {
	$line .= "\">".($bin*$i)."\",";
	$i++;
	}
chop $line;
open (CMDR, ">$outdir/hist_temp.R") || die;
print CMDR "#!/usr/bin/env Rscript\n\n" ;
##png:
#print CMDR "png(\"$outdir/covPlot/$outName\_allCovPlot.png\", 800, ".(scalar@Files*250).")
#par(mfrow=c(".scalar@Files.",1))
#mat <- read.table(\"$outdir/covPlot/$outName\_sumCovPlot.txt\",header=T, row.names=1)\n";
#foreach my$file (@Files) {
#	print CMDR "barplot(t(matrix(mat[\"sample_${$smplName}{$file}\",])), names.arg=c(";
#	print CMDR $line;
#	print CMDR "), ylim=c(0,1), col=\"blue\", axis.lty=1, main=\"$${$smplName}{$file}\",xlab=\"Depth\", ylab=\"fraction of $outName\", cex.lab=1.2, cex.main=1.5, #cex.axis=1.2, cex.names=1.2)\n";
#	}
##pdf:
print CMDR "pdf(\"$outdir/$outName\_allCovPlot.pdf\", width=11.69, height=4.135)
mat <- read.table(\"$outdir/$outName\_sumCovPlot.txt\",header=T, row.names=1)\n";
foreach my$file (@{$Files}) {
	print CMDR "barplot(t(matrix(mat[\"sample_".${$smplName}{$file}."\",])), names.arg=c(";
	print CMDR $line;
	print CMDR "), ylim=c(0,1), col=\"blue\", axis.lty=1, main=\"".${$smplName}{$file}."\",xlab=\"Depth\", ylab=\"fraction of $outName\", cex.lab=1.2, cex.main=1.5, cex.axis=1.2, cex.names=1.2)\n";
	}
print CMDR "dev.off();warnings();\n";
close CMDR;
system "Rscript $outdir/hist_temp.R";
unlink "$outdir/hist_temp.R";
}


##########
#R barplot -bin , mean of all samples
#binCovPlot($bin,$maxPl,$outdir,$outName,$lengthBed);

sub binCovPlot {
my($bin,$maxPl,$outdir,$outName,$lengthBed) = @_;
print "doing binCovPlot\n";
my$line="";
my$i=0;
while ( ($i*$bin) <= $maxPl ) {
	$line .= "\"".($bin*$i)."-\",";
	$i++;
	}
$line .= "\">100\"";
open (CMDR, ">$outdir/plot_temp.R") || die;
print CMDR "#!/usr/bin/env Rscript\n\n" ;
print CMDR "png(\"$outdir/covPlot/$outName\_binCovPlot.png\", 800, 500)
mat <- read.table(\"$outdir/covPlot/$outName\_binCovPlot.txt\",header=T, row.names=1)
mean1 = apply(mat,2,mean)
mean2 = mean1[1:length(mean1)-1] 
std1 = apply(mat,2,sd)
std2 = std1[1:length(std1)-1]
error.bar <- function(x, y, upper, lower=upper, length=0.1,...) {
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop(\"vectors must be same length\")
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...) }
barx <- barplot(mean1, names.arg=c(";
print CMDR $line;
print CMDR "), col=\"blue\", axis.lty=1, main=\"mean of all samples\", xlab=\"Depth\", ylab=\"fraction of $outName ($lengthBed bp)\")
error.bar(barx,mean1, 1.96*std1/10)
dev.off();warnings();\n";
close CMDR;
system "Rscript $outdir/plot_temp.R";
unlink "$outdir/plot_temp.R";
}


##########
#R barplot -bin , 1graph / sample
#allBinCovPlot($bin,$maxPl,$outdir,$outName,$lengthBed,\@Files,\%sName2);

sub allBinCovPlot {
my($bin,$maxPl,$outdir,$outName,$lengthBed,$h1,$h2) = @_;
my@Files = @$h1;
my%sampleName = %$h2;
print "doing allBinCovPlot\n";
my$line="";
my$i=0;
while ( ($i*$bin) <= $maxPl )
	{
	$line .= "\"".($bin*$i)."-\",";
	$i++;
	}
$line .= "\">100\"";
open (CMDR, ">$outdir/plot_temp.R") || die;
print CMDR "#!/usr/bin/env Rscript\n\n" ;
print CMDR "png(\"$outdir/covPlot/$outName\_allBinCovPlot.png\", 800, ".(scalar@Files*250).")
par(mfrow=c(".scalar@Files.",1))
mat <- read.table(\"$outdir/covPlot/$outName\_binCovPlot.txt\",header=T, row.names=1)\n";
foreach my$file (@Files)
	{
	print CMDR "barplot(t(matrix(mat[\"sample_$sampleName{$file}\",])), names.arg=c(";
	print CMDR $line;
	print CMDR "), ylim=c(0,1), col=\"blue\", axis.lty=1, main=\"$sampleName{$file}\",xlab=\"Depth\", ylab=\"fraction of $outName\", cex.lab=1.2, cex.main=1.5, cex.axis=1.2, cex.names=1.2)\n";
	}

print CMDR "dev.off();warnings();\n";
close CMDR;
system "Rscript $outdir/plot_temp.R";
unlink "$outdir/plot_temp.R";
}


##########
#depthBoxPlot($outdir,$outName,$Files,$smplName,$threshold,$lengthBed,$depthCount,$covBases,$totBases);

sub depthBoxPlot {

my($outdir,$outName,$Files,$smplName,$threshold,$lengthBed,$depthCount,$covBases,$totBases) = @_;

print "\n\t## doing depth boxPlot\n";
my%Vals;
foreach my$file (@{$Files}) {
	my@Depths = sort{$a<=>$b}(keys%{ ${$depthCount}{$file} });
	$Vals{"min"}{$file} = $Depths[0];
	$Vals{"max"}{$file} = $Depths[-1];
	$Vals{"median"}{$file} = medianFromHistog($lengthBed,\@Depths,\%{ ${$depthCount}{$file} });
	my($q1Length,@q1Set,$q3Length,@q3Set);
	foreach (@Depths) {
		if ($_ < $Vals{"median"}{$file}) {
			push(@q1Set, $_);
			$q1Length += ${$depthCount}{$file}{$_};
			}
		elsif ($_ > $Vals{"median"}{$file}) {
			push(@q3Set, $_);
			$q3Length += ${$depthCount}{$file}{$_};
			}
		}
	$Vals{"Q1"}{$file} = medianFromHistog($q1Length,\@q1Set,\%{ ${$depthCount}{$file} });
	$Vals{"Q3"}{$file} = medianFromHistog($q3Length,\@q3Set,\%{ ${$depthCount}{$file} });
	$Vals{"lower_wsk"}{$file} = $Vals{"Q1"}{$file} - 1.5*($Vals{"Q3"}{$file}-$Vals{"Q1"}{$file});
	if ($Vals{"lower_wsk"}{$file} < $Vals{"min"}{$file}) { $Vals{"lower_wsk"}{$file} = $Vals{"min"}{$file} ; }
	$Vals{"upper_wsk"}{$file} = $Vals{"Q3"}{$file} + 1.5*($Vals{"Q3"}{$file}-$Vals{"Q1"}{$file});
	if ($Vals{"upper_wsk"}{$file} > $Vals{"max"}{$file}) { $Vals{"upper_wsk"}{$file} = $Vals{"max"}{$file} ; }
	}

my@fracOfMean = (10,5,2);
open(my$fhOut, ">", "$outdir/$outName\_depthReport.txt") || die "can't create file $outdir/$outName\_depthReport.txt($!)\n";
print $fhOut "$outName.bed :\n\n";
print $fhOut "total length : $lengthBed bp\n\n";
print $fhOut "samples\tmin\tQ1\tmedian\tQ3\tmax\tmean";
foreach my$frac (@fracOfMean) { print $fhOut "\t%"."cov <=(mean depth/$frac)x"; }
if ($threshold) { print $fhOut "\t%"."cov >=".$threshold."x"; }
print $fhOut "\n";
foreach my$file (@{$Files}) {
	print $fhOut ${$smplName}{$file}."\t".$Vals{"min"}{$file}."\t".$Vals{"Q1"}{$file}."\t".$Vals{"median"}{$file}."\t".$Vals{"Q3"}{$file}."\t".$Vals{"max"}{$file};
	my$totMean = ${$totBases}{$file} / $lengthBed;
	print $fhOut "\t".sprintf("%.1f", $totMean);
	my%fracMean;
	foreach my$frac (@fracOfMean) {
		my$fracMean = 0;
		for (my$d=0;$d<=int($totMean / $frac);$d++) {
			if (exists ${$depthCount}{$file}{$d}) { $fracMean += ${$depthCount}{$file}{$d}; }
			}
		$fracMean{$frac} = $fracMean / $lengthBed;
		print $fhOut "\t".100*(sprintf("%.3f", $fracMean{$frac}));
		}
	if ($threshold && exists ${$covBases}{$file}) { print $fhOut "\t".100*(sprintf("%.3f", (${$covBases}{$file} / $lengthBed))); }
	else { print $fhOut "\t0",}
	print $fhOut "\n";
	}
close $fhOut;



my$Rtxt = "#!/usr/bin/env Rscript\n
pdf(\"$outdir/$outName\_depth_boxplot.pdf\", width=11.69, height=4.135)\n";
#Vals <- list( stats=matrix(c(lower_whisker{1},Q1{1},Median{1},Q3{1},upper_whisker{1},lower_whisker{2},Q1{2},Median{2},Q3{2},upper_whisker{2})), ncol=2) , n=c(tot,tot) , out=c(min{1},max{1},min{2},max{2}) , group=c(1,1,2,2), names=c(\"smpl1\",\"smpl2\") )
$Rtxt .= "Vals <- list( stats=matrix(c(";
foreach (@{$Files}) { $Rtxt .= $Vals{"lower_wsk"}{$_}.",".$Vals{"Q1"}{$_}.",".$Vals{"median"}{$_}.",".$Vals{"Q3"}{$_}.",".$Vals{"upper_wsk"}{$_}.","; }
chop $Rtxt;
$Rtxt .= "),ncol=".scalar(@{$Files}).") , n=c(";
foreach (@{$Files}) { $Rtxt .= $lengthBed.","; }
chop $Rtxt;
$Rtxt .= ") , out=c(";
foreach (@{$Files}) { $Rtxt .= $Vals{"min"}{$_}.",".$Vals{"max"}{$_}.","; }
chop $Rtxt;
$Rtxt .= ") , group=c(";
my$i = 1;
foreach (@{$Files}) { $Rtxt .= $i.",".$i.","; $i++; }
chop $Rtxt;
$Rtxt .= ") , names=c(";
foreach (@{$Files}) { $Rtxt .= "\"".${$smplName}{$_}."\","; }
chop $Rtxt;
$Rtxt .= ") )
bxp(Vals, outline=TRUE , show.names = TRUE)
dev.off();\n";

open ($fhOut, ">", "$outdir/plot_temp.R") || die "can't create $outdir/plot_temp.R ($!)\n";
print $fhOut "$Rtxt";
close($fhOut);
system "Rscript $outdir/plot_temp.R";
unlink "$outdir/plot_temp.R";

}
#####
sub medianFromHistog {
my($length,$Depths,$depthCount) = @_;
my$median;
my$countD = 0;
#odd?
if($length%2) {
	foreach my$d (@{$Depths}) {
		$countD += ${$depthCount}{$d};
		if ($countD > $length/2) {
			$median = $d;
			last;
			}
		}
	}
#even
else {
	for my$i (0..$#{$Depths}) {
		$countD += ${$depthCount}{${$Depths}[$i]};
		if ($countD >= $length/2) {
			if ($countD == $length/2) {
				$median = (${$Depths}[$i] + ${$Depths}[$i+1]) / 2;
				last;
				}
			else {
				$median = ${$Depths}[$i];
				last;
				}
			}
		}
	}
return($median);

}



############################
#if (!$bedFile || $jobs{"genePlot"}) { 
#	InterS($bin,$maxPl,$maxGr,$outPlotDir,"all_Genes_Regions",\@Files,$smplIdx,\%{ ${$depthFilePerChr}{"raw"} },$coordBed_r); 
#	}
#if ($bedFile) {		
#	InterS($bin,$maxPl,$maxGr,$outPlotDir,$bedName,\@Files,$smplIdx,\%{ ${$depthFilePerChr}{"raw"} },$Bed_r);
#	}

sub InterS {

my($bin,$maxPl,$maxGr,$outdir,$outName,$Files_r,$smplIdx,$depthFilePerChr,$intervals_r)=@_;
print "\n\n####\n\ndoing interPlot for $outName\n\n####\n\n";

if (scalar(keys%{$intervals_r}) != 0) {

	my@allBins; my$i=0;
	while ( ($i*$bin) <= $maxPl ) {
		push(@allBins,($bin*$i+1));
		$i++;
		}
	my(%nPos,%nCov);
	foreach (@allBins) {
		$nPos{$_}=0;
		$nCov{$_}=0;
		}

	#intersection by chr
	foreach my$chr (keys%{$intervals_r}) {
		#initialyzes %regCov
		my%regCov;
		foreach (@allBins) {
			my$r=0; #idx of reg of $NM
			foreach my$startReg (sort{$a<=>$b}keys%{ ${$intervals_r}{$chr} }) {
				my$pos=$startReg;
				while ($pos <= ${$intervals_r}{$chr}{$startReg} ) {
					push (@{ $regCov{$_}{$r} }, 0);
					$pos++;
					}
				$r++;
				}
			}
		## %allDepth: depth foreach position in intervals from %allInterval : ${$allDepth}{$file}{$loc} = depth
		print "\n\t##intersecting plot target regions with $chr depth File\n";
		my$allDepth_r = DeCovA::covByRegion::intersectCovFile($Files_r,$smplIdx,${$depthFilePerChr}{$chr},$maxGr,\@allBins,\%{ ${$intervals_r}{$chr} },"");

		foreach my$file (@{$Files_r}) {
			print "\t\tanalysing $file\n";
			# -> @{ $regCov{$threshold}{$nReg} } = [nber of covered samples foreach pos]
			foreach my$threshold (@allBins) {
				DeCovA::covByRegion::covByThreshold($threshold,\%{ ${$intervals_r}{$chr} },\%{ ${$allDepth_r}{$file} },\%{ $regCov{$threshold} });
				}
			}

		foreach my$b (@allBins) {
			foreach my$r (keys%{ $regCov{$b} }) {
				foreach my$d (@{ $regCov{$b}{$r} }) {	#print "$b $r $d\n";
					$nPos{$b}++;
					if ($d == scalar@{$Files_r}) { $nCov{$b}++; }
					}
				}
			}
		}
	#plot
	my$histTxt = "depth\t";
	$i=0;
	while ( ($i*$bin) <= $maxPl ) {
		$histTxt .= ">".($bin*$i)."\t";
		$i++;
		}
	chop $histTxt;
	$histTxt .= "\ncov\t";
	foreach (@allBins) {
		$histTxt .=  ($nCov{$_}/$nPos{$_})."\t";
		}
	chop $histTxt; $histTxt .= "\n";
	open(my$fh, ">", "$outdir/$outName\_interPlot.txt") || die "can't create $outdir/$outName\_interPlot.txt ($!)\n";
	print $fh $histTxt;
	close($fh); 

	#R barplot -sum, 1graph / sample
	interPlot($bin,$maxPl,$nPos{$allBins[0]},$outdir,$outName);

	}

}


##########
sub interPlot {
my($bin,$maxPl,$lengthBed,$outdir,$outName) = @_;
print "\n\t## doing inter-samples CovPlot\n";
my$line="";
my$i=0;
while ( ($i*$bin) <= $maxPl ) {
	$line .= "\">".($bin*$i)."\",";
	$i++;
	}
chop $line;
open (CMDR, ">$outdir/hist_temp.R") || die;
print CMDR "#!/usr/bin/env Rscript\n\n" ;
print CMDR "pdf(\"$outdir/$outName\_interPlot.pdf\", width=11.69, height=4.135)
mat <- read.table(\"$outdir/$outName\_interPlot.txt\",header=T, row.names=1)\n";
print CMDR "barplot(t(matrix(mat[\"cov\",])), names.arg=c(";
print CMDR $line;
print CMDR "), ylim=c(0,1), col=\"blue\", axis.lty=1, main=\"intersection of all samples\",xlab=\"Depth\", ylab=\"fraction of $outName ($lengthBed bp)\", cex.lab=1.2, cex.main=1.5, cex.axis=1.2, cex.names=1.2)\n";
print CMDR "dev.off();warnings();\n";
close CMDR;
system "Rscript $outdir/hist_temp.R";
unlink "$outdir/hist_temp.R";
}



########################




1;

