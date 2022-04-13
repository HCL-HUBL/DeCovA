package DeCovA::bedCov;


use strict;
use warnings;


################################

##$bedLines_r = Bed_Cov(\@Files,\%fName,\%sName2,$outdir,$bedFile,$bedName,\%Bed,$pThreshold,$toReseq,\@covFields,$depthFilePerChr,$smplIdx,\%CNV_opt);
##$bedLines_r = Bed_Cov(\@Files,\%fName,\%sName2,$outdir,$id2Bed,$bedName,\%coordBed,$pThreshold,$toReseq,\@covFields,$depthFilePerChr,$smplIdx,\%CNV_opt);

sub Bed_Cov {

	my ($Files_r,$fName_r,$sName_r,$outdir,$bedFile,$bedName,$mergeBed_r,$threshold,$toReseq,$bedReport,$covFields_r,$depthFilePerChr,$smplIdx_r,$CNV_opt_r) = @_;

	## in bedFile, some intervals can overlap each other
	## in $mergeBed_r, intervals have been merged (-> only 1 end per start); #$Bed{$chr}{$start} = $end;
	print "\n\n####\n\nperform bed coverage by interval";
	if ($threshold) { print " (>=$threshold X)"; }
	print "\n\n####\n\n";

	my (@allLines,@bedLines,%Chr,%Start,%End,%Intervals,%Intervals2);
	my $i = 0; my $j = 0;
	open(BED, "$bedFile") || die "can't open file $bedFile\n";
	while (my $line = <BED>) {
		$line =~ s/\s+$//;
		$allLines[$i] = $line;
		if ( ($line =~ /^\w+\t\d+\t\d+/) && ($line !~ /^#/) ) {
			my @tab = split(/\t/,$line);
			if ($tab[1] < $tab[2]) {
				my $chr = $tab[0];
				$chr =~ s/chr//i;
				$bedLines[$j]{"allLine"} = $line;
				$bedLines[$j]{"Chrom"} = $chr;
				$bedLines[$j]{"Start"} = $tab[1]+1;		## -> 1-based
				$bedLines[$j]{"End"} = $tab[2];
				if ($tab[3] && $tab[3] !~ /^\s*$/) { $bedLines[$j]{"Infos"} = $tab[3]; }
				$Chr{$i} = $chr;
				$Start{$i} = $tab[1]+1;		## -> 1-based
				$End{$i} = $tab[2];
				$Intervals{$Chr{$i}}{$Start{$i}}{$End{$i}} = 1;
				$j++;
			}
		}
		$i++;
	}
	close BED;
	## add merged intervals, to compute totDepth
	foreach my $chr (keys%{$mergeBed_r}) {
		foreach my $start (keys%{ ${$mergeBed_r}{$chr} }) {
			$Intervals{$chr}{$start}{${$mergeBed_r}{$chr}{$start}} = 1;
		}
	}
	## intervals2 : array of increasing ends foreach start
	foreach my $chr (keys%Intervals) {
		foreach my $start (keys%{ $Intervals{$chr} }) {
			@{ $Intervals2{$chr}{$start} } = sort{$a<=>$b}(keys%{ $Intervals{$chr}{$start} });
		}
	}

	## initialize vals
	my $notCov = {};
	my $TotByDepth = {};
	my ($min,$max,$sum,$mean,$median,$cov,$TotBases,$TotByCov);
	foreach my $file (@{$Files_r}) {
		${$TotByCov}{$file} = 0;
		${$TotBases}{$file} = 0;
		foreach my $chr (keys%Intervals) {
			foreach my $start (keys%{ $Intervals{$chr} }) {
				foreach my $end (@{ $Intervals2{$chr}{$start} }) {
					${$min}{$file}{$chr}{$start}{$end} = 0;
					${$max}{$file}{$chr}{$start}{$end} = 0;
					${$mean}{$file}{$chr}{$start}{$end} = 0;
					${$sum}{$file}{$chr}{$start}{$end} = 0;
					${$median}{$file}{$chr}{$start}{$end} = 0;
					${$cov}{$file}{$chr}{$start}{$end} = 0;
				}
			}
		}
	}
	##intersect with each $depthFilePerChr
	foreach my $chr (keys%Intervals) {
		my $c = 0;			#idx of $startByChr
		my @Starts = sort{$a<=>$b}(keys%{ $Intervals{$chr} });
		my ($noCovStart,$noCovEnd,$vals);
		foreach my $start (@Starts) {
			foreach my $end (@{ $Intervals2{$chr}{$start} }) {
				foreach (@{$Files_r}) {
					${$noCovStart}{$_}{$start}{$end} = 0;
					${$noCovEnd}{$_}{$start}{$end} = 0;
					@{ ${$vals}{$_}{$start}{$end} } = ();
				}
			}
		}
		open(my $fh, "<", ${$depthFilePerChr}{$chr}) || die "can't read file ".${$depthFilePerChr}{$chr}." ($!)\n";
		while (my $line = <$fh>) {
			chomp $line;
			my @tab = split(/\t/,$line);
			my $pos = $tab[0];
			while ( ($pos > $Intervals2{$chr}{$Starts[$c]}[-1]) && ($c < $#Starts) ) { $c++; }
			my $c2 = $c;			##@Starts iteration within while loop
			while ( ($pos >= $Starts[$c2]) && ($c2 < $#Starts) ) {
				fillVals($Files_r,$smplIdx_r,$mergeBed_r,$threshold,$chr,$pos,$Starts[$c2],$Intervals2{$chr}{$Starts[$c2]},$cov,$notCov,$noCovStart,$noCovEnd,\@tab,$min,$max,$sum,$mean,$median,$TotBases,$TotByCov,$TotByDepth,$vals);
				$c2++;
			}
			if ($pos >= $Starts[$c2]) {
				fillVals($Files_r,$smplIdx_r,$mergeBed_r,$threshold,$chr,$pos,$Starts[$c2],$Intervals2{$chr}{$Starts[$c2]},$cov,$notCov,$noCovStart,$noCovEnd,\@tab,$min,$max,$sum,$mean,$median,$TotBases,$TotByCov,$TotByDepth,$vals);
			}
			if ( $pos > $Intervals2{$chr}{$Starts[-1]}[-1] ) { last; }
		}
		close $fh;
	}

	my $meanTotBases = 0;
	foreach (@{$Files_r}) { $meanTotBases += ${$TotBases}{$_}; }
	if ($meanTotBases) {
		$meanTotBases /= scalar@{$Files_r};
	} else {
		print "no bases sequenced within bed in all samples!\n";
	}

	my $TotLength = 0;
	foreach my $chr (keys%{$mergeBed_r}) {
		foreach my $start (keys%{ ${$mergeBed_r}{$chr} }) {
			$TotLength += (${$mergeBed_r}{$chr}{$start} - $start + 1);
		}
	}

	my $nCol;
	foreach (@allLines) {
		if (($_ =~ /^\w+\t\d+\t\d+/) && ($_ !~ /^#/)) {
			my @Col = split(/\t/,$_);
			$nCol = scalar@Col;
		}
		if ($nCol) { last; }
	}

	open(my $fh1, ">", "$outdir/$bedName.cov.txt") || die "can't create file $outdir/$bedName.cov.txt ($!)\n";
	print $fh1 "##samples:";
	for (my $i=0;$i<($nCol-1);$i++) { print $fh1 "\t"; }
	foreach (@{$Files_r}) { 
		print $fh1 "\t${$sName_r}{$_}";
		for (my $i=1;$i<scalar@{$covFields_r};$i++) { print $fh1 "\t"; }
	}
	print $fh1 "\n##bed col.";
	for (my $i=0;$i<($nCol-1);$i++) { print $fh1 "\t"; }
	foreach (@{$Files_r}) {
		foreach (@{$covFields_r}) {
			if ($_ eq "cov") {
				print $fh1 "\t%>=$threshold"."X";
			} else {
				print $fh1 "\t$_";
			}
		}
	}
	print $fh1 "\n";
	my (%reseq,$fh2);
	if ($toReseq) {
		open($fh2, ">", "$outdir/$bedName.toReseq.txt") || die "can't create file $outdir/$bedName.toReseq.txt ($!)\n";
	}
	$j=0;
	for (my $i=0;$i<scalar@allLines;$i++) {
		print $fh1 $allLines[$i];
		if (exists $Chr{$i}) {
			foreach my $file (@{$Files_r}) {
				$bedLines[$j]{$file}{"mean"} = ${$mean}{$file}{$Chr{$i}}{$Start{$i}}{$End{$i}};
				if (exists ${$CNV_opt_r}{"seuil_cov"}) { $bedLines[$j]{$file}{"cov"} = ${$cov}{$file}{$Chr{$i}}{$Start{$i}}{$End{$i}}; }
				if (${$CNV_opt_r}{"RefDepth"} eq "tot") { $bedLines[$j]{$file}{"tot"} = ${$sum}{$file}{$Chr{$i}}{$Start{$i}}{$End{$i}}; }
				foreach (@{$covFields_r}) {
					if ($_ eq "min") {
						print $fh1 "\t".${$min}{$file}{$Chr{$i}}{$Start{$i}}{$End{$i}};
					} elsif ($_ eq "max") {
						print $fh1 "\t".${$max}{$file}{$Chr{$i}}{$Start{$i}}{$End{$i}};
					} elsif ($_ eq "tot") {
						print $fh1 "\t".${$sum}{$file}{$Chr{$i}}{$Start{$i}}{$End{$i}};
					} elsif ($_ eq "mean") {
						print $fh1 "\t".sprintf("%.1f",${$mean}{$file}{$Chr{$i}}{$Start{$i}}{$End{$i}});
					} elsif ($_ eq "median") {
						print $fh1 "\t".sprintf("%.1f",${$median}{$file}{$Chr{$i}}{$Start{$i}}{$End{$i}});
					} elsif ($_ eq "cov") {
						print $fh1 "\t".100*(sprintf("%.3f", ${$cov}{$file}{$Chr{$i}}{$Start{$i}}{$End{$i}}))."%";
					}
				}
				if ($toReseq && ${$cov}{$file}{$Chr{$i}}{$Start{$i}}{$End{$i}} < $toReseq) {
					push(@{ $reseq{$i} },"${$sName_r}{$file} (".100*(sprintf("%.3f", ${$cov}{$file}{$Chr{$i}}{$Start{$i}}{$End{$i}}))."%)");
					#print STDERR "$toReseq ${$sName_r}{$file} $i $cov{$file}{$Chr{$i}}{$Start{$i}}{$End{$i}}\n";
				}
			}
			$j++;
		}
		print $fh1 "\n";
		if ($toReseq && exists $reseq{$i}) {
			print $fh2 $allLines[$i];
			foreach (@{ $reseq{$i} }) { print $fh2 "\t$_"; }
			print $fh2 "\n";
		}
	}
	close $fh1; 
	if ($toReseq) { close $fh2; }

	if ($bedReport && $threshold) {
		foreach my $file (@{$Files_r}) {
			open($fh1, ">$outdir/cov\_${$sName_r}{$file}/$bedName.${$sName_r}{$file}.notCov.txt") || die "can't create file $outdir/cov\_${$sName_r}{$file}/$bedName.${$sName_r}{$file}.notCov.bed\n";
			print $fh1 "depthAnd Coverage report on $bedName.bed for ${$sName_r}{$file} sample\n\n";
			print $fh1 "total length of bed : $TotLength bp\n\n";
			print $fh1 "total cov >=$threshold"."x : ".100*(sprintf("%.3f", (${$TotByCov}{$file}/$TotLength)))." %\n\n";
			print $fh1 "mean depth : ".sprintf("%.1f", (${$TotBases}{$file}/$TotLength))." x\n\n";
			print $fh1 "domains covered less than $threshold"."x\n";
			for (my $i=0;$i<scalar@allLines;$i++) {
				if (exists $Chr{$i}) {
					if (exists ${$notCov}{$file}{$Chr{$i}}{$Start{$i}}{$End{$i}}) {
						print $fh1 "\n".$allLines[$i]."\n";
						foreach (sort{$a<=>$b}keys%{ ${$notCov}{$file}{$Chr{$i}}{$Start{$i}}{$End{$i}} }) {
							print $fh1 "\t".($_-1)."\t".${$notCov}{$file}{$Chr{$i}}{$Start{$i}}{$End{$i}}{$_}."\n";
						}
					}
				}
			}
			close $fh1;
		}
	}

	return(\@bedLines,$TotLength,$TotByDepth,$TotByCov,$TotBases);

}


## fillVals($Files_r,$smplIdx_r,$threshold,$chr,$pos,$Starts[$c2],$Intervals2{$chr}{$Starts[$c2]},$cov,$notCov,$noCovStart,$noCovEnd,\@tab,$min,$max,$sum,$mean,$median,$TotBases,$TotByCov,$TotByDepth,$vals);
sub fillVals {
	my ($Files_r,$smplIdx_r,$mergeBed_r,$threshold,$chr,$pos,$start,$Ends,$cov,$notCov,$noCovStart,$noCovEnd,$fields,$min,$max,$sum,$mean,$median,$TotBases,$TotByCov,$TotByDepth,$vals) = @_;
	if ($pos == $start) {
		foreach my $end (@{$Ends}) {
			foreach (@{$Files_r}) {
				${$min}{$_}{$chr}{$start}{$end} = ${$fields}[${$smplIdx_r}{$_}];
				${$max}{$_}{$chr}{$start}{$end} = ${$fields}[${$smplIdx_r}{$_}];
			}
		}
	}
	foreach my $end (@{$Ends}) {
		if ($pos <= $end) {
			getVals($Files_r,$smplIdx_r,$threshold,$chr,$pos,$start,$end,$cov,$notCov,$noCovStart,$noCovEnd,$fields,$min,$max,$vals);
		}
		if ($pos == $end) {
			computeVals($Files_r,$mergeBed_r,$threshold,$chr,$pos,$start,$end,$cov,$sum,$mean,$median,$TotBases,$TotByCov,$TotByDepth,$vals);
		}
	}
}

## getVals($Files_r,$smplIdx_r,$threshold,$chr,$pos,$start,$end,$cov,$notCov,$noCovStart,$noCovEnd,$fields,$min,$max,$vals);
sub getVals {
	my ($Files_r,$smplIdx_r,$threshold,$chr,$pos,$start,$end,$cov,$notCov,$noCovStart,$noCovEnd,$fields,$min,$max,$vals) = @_;
	foreach (@{$Files_r}) {
		if ($threshold && ${$fields}[${$smplIdx_r}{$_}] >= $threshold) {
			${$cov}{$_}{$chr}{$start}{$end}++;
		} else {
			if ( $pos == (${$noCovEnd}{$_}{$start}{$end} + 1) ) { 
				${$notCov}{$_}{$chr}{$start}{$end}{${$noCovStart}{$_}{$start}{$end}} = $pos ;
				${$noCovEnd}{$_}{$start}{$end} = $pos; 
			} else { 
				${$notCov}{$_}{$chr}{$start}{$end}{$pos} = $pos;
				${$noCovStart}{$_}{$start}{$end} = $pos;
				${$noCovEnd}{$_}{$start}{$end} = $pos;
			}
		}
		if (${$fields}[${$smplIdx_r}{$_}] < ${$min}{$_}{$chr}{$start}{$end}) {
			${$min}{$_}{$chr}{$start}{$end} = ${$fields}[${$smplIdx_r}{$_}];
		}
		if (${$fields}[${$smplIdx_r}{$_}] > ${$max}{$_}{$chr}{$start}{$end}) {
			${$max}{$_}{$chr}{$start}{$end} = ${$fields}[${$smplIdx_r}{$_}];
		}
		push(@{ ${$vals}{$_}{$start}{$end} }, ${$fields}[${$smplIdx_r}{$_}]);
	}
}

## computeVals($Files_r,$mergeBed_r,$threshold,$chr,$pos,$start,$end,$cov,$sum,$mean,$median,$TotBases,$TotByCov,$TotByDepth,$vals);
sub computeVals {
	my ($Files_r,$mergeBed_r,$threshold,$chr,$pos,$start,$end,$cov,$sum,$mean,$median,$TotBases,$TotByCov,$TotByDepth,$vals) = @_;
	my $len = $end - $start + 1;
	foreach my $f (@{$Files_r}) {
		#cov
		${$cov}{$f}{$chr}{$start}{$end} /= $len;
		#sum
		${$sum}{$f}{$chr}{$start}{$end} = 0;
		foreach (@{ ${$vals}{$f}{$start}{$end} }) {
			${$sum}{$f}{$chr}{$start}{$end} += $_;
		}
		#mean
		${$mean}{$f}{$chr}{$start}{$end} = (${$sum}{$f}{$chr}{$start}{$end} / $len);
		#median
		my @sortDepth=();
		for (my $i=0;$i<($len - scalar@{ ${$vals}{$f}{$start}{$end} });$i++) { push(@sortDepth,0); }
		push(@sortDepth,@{ ${$vals}{$f}{$start}{$end} });
		@sortDepth = sort{$a<=>$b}@sortDepth;
		if($len % 2) { #odd?
			${$median}{$f}{$chr}{$start}{$end} = $sortDepth[int($len/2)];
		} else { #even
			${$median}{$f}{$chr}{$start}{$end} = (($sortDepth[int($len/2)-1]+$sortDepth[int($len/2)])/2 );
		}
	}
	##just for non-overlapping intervals
	if (exists ${$mergeBed_r}{$chr}{$start} && $end == ${$mergeBed_r}{$chr}{$start}) {
		foreach my $f (@{$Files_r}) {
			foreach my $d (@{ ${$vals}{$f}{$start}{$end} }) {
				${$TotBases}{$f} += $d;
				${$TotByDepth}{$f}{$d}++;
				if ($threshold && $d >= $threshold) { ${$TotByCov}{$f}++; }
			}
		}
	}
	foreach my $f (@{$Files_r}) {
		delete ${$vals}{$f}{$start}{$end};
	}
}

################################



1;
