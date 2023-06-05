package DeCovA::bedProcessing;


use strict;
use warnings;



######################
sub readBed {

	my $bedFile = shift;
	my %covBed;			# $covbed{$chr}{$start} = $end;	transform in hash with 1-based coord
	open(my $fh, "<", $bedFile) || die "can't open file $bedFile\n";
	print "\n\treading $bedFile\n";
	while (my $line = <$fh>) {
		if ( ($line =~ /^\w+\t\d+\t\d+/) && ($line !~ /^#/) ) {
			chomp $line;
			my @tab = split(/\t/,$line);
			$tab[0] =~ s/^chr//i;
			# keep longest interval
			if ( exists $covBed{$tab[0]}{($tab[1]+1)} ) {
				if ( $tab[2] > $covBed{$tab[0]}{($tab[1]+1)}) {
					$covBed{$tab[0]}{($tab[1]+1)} = $tab[2];
				}
			} else {
				$covBed{$tab[0]}{($tab[1]+1)} = $tab[2];
			}
		}
	}
	close ($fh);
	unless (%covBed) { die "\n!! no bed formatted line in $bedFile file\n"; }
	return(\%covBed);

}


#########################

sub readMut {

	my ($file,$chromName) = @_;
	my %Mut = ();			# $Mut{$chr}{$start} = $end;
	my $fh;
	if ($file =~ /.gz$/) {
		use IO::Zlib;
		$fh = new IO::Zlib;
		$fh->open($file, "rb")
	} else {
		open($fh, "<", $file) or die "could not read $file ($!)\n";
	}
	while (my $line = <$fh>) {
		unless ($line =~ /^\s*$/ || $line =~ /^#/) {
			$line =~ s/\s+$//;
			my @tab = split(/\t/,$line);
			my $chr = $tab[0];
			$chr  =~ s/chr//i;
			$Mut{$chr}{$tab[1]} = $tab[2];	# already in 1-based coord
			${$chromName}{"mut"}{$chr} = $tab[0];
		}
	}
	close $fh;
	return(\%Mut);

}


######################
#
sub getChrNameInBed {

	my ($bedFile,$chromOrder_r,$chromName_r) = @_;
	# bedFile in %intervals
	my (@tmpOrder,%allChrom);
	open(my $fh, "<", $bedFile) || die "can't open file $bedFile\n";
	while (my $line = <$fh>) {
		if ( ($line =~ /^\w+\t\d+\t\d+/) && ($line !~ /^#/) ) {
			chomp $line;
			my @tab = split(/\t/,$line);
			unless (exists $allChrom{$tab[0]}) {
				$allChrom{$tab[0]} = 1;
				my $chr = $tab[0];
				$chr =~ s/^chr//i;
				push(@tmpOrder,$chr);
				${$chromName_r}{"bed"}{$chr} = $tab[0];
			}
		}
	}
	close ($fh);
	unless (@tmpOrder) { die "\n!! no bed formatted line in $bedFile file\n"; }
	unless (@{$chromOrder_r}) {
		@{$chromOrder_r} = @tmpOrder;
	}

}

######################
#genes and transcipts and exon nbr, and if none, introns or intergenic (closer gene and distance to it)
sub reAnnotBed3 {

	# print"re_annotating with transcripts overlapping bed file\n";
	my ($bedFile,$reAnnot_r,$RefID,$chromOrder_r,$chromName_r,$nonCod) = @_;	#refseq is in 0-based		if (exists ${$IDrf}{$nm})

	my %RefCoord;		# $RefCoord{chrom}{start_nm}{end_nm} = [nm1, nm2, ...]
	foreach my $nm (keys%{$RefID}) {
		if( $nonCod || (exists ${$RefID}{$nm}{"CDS_start"}) ) {
			if (exists ${$RefID}{$nm}{"ex_starts"}) {
				push(@{ $RefCoord{${$RefID}{$nm}{"chr"}}{${$RefID}{$nm}{"ex_starts"}[0]}{${$RefID}{$nm}{"ex_ends"}[-1]} }, $nm);
			}
		}
	}


	## bedFile in %intervals
	my %intervals;
	open(my $fh, "<", $bedFile) || die "can't open file $bedFile\n";
	while (my $line = <$fh>) {
		if (($line =~ /^\w+\t\d+\t\d+/)&&($line !~ /^#/)) {
			chomp $line;
			my @tab = split(/\t/,$line);
			my $chr = $tab[0];
			$chr =~ s/^chr//i;
			$intervals{$chr}{$tab[1]}{$tab[2]}{"info"} = "";
			if ($tab[4]) {
				#for my$i (4..$#tab) { $intervals{$chr}{$tab[1]}{$tab[2]}{"info"} .= "\t".$tab[$i]; }
				$intervals{$chr}{$tab[1]}{$tab[2]}{"info"} = "\t".join("\t",@tab[4..$#tab]);
			}
		}
	}
	close ($fh);
	unless (%intervals) { die "\n!! no bed formatted line in $bedFile file\n"; }

	## intersecting both interval hashes
	foreach my $chr (keys%intervals) {
		if (exists $RefCoord{$chr}) {
			my @bedStarts = sort{$a<=>$b}keys%{ $intervals{$chr} };
			my @refStarts = sort{$a<=>$b}keys%{ $RefCoord{$chr} };
			my $c = 0;	#idx of @refStarts
			my @RefEndC = sort{$b<=>$a}keys%{ $RefCoord{$chr}{$refStarts[$c]} };
			foreach my $bedstart (@bedStarts) {
				while ( ($c < $#refStarts) && ($bedstart > $RefEndC[0]) ) {
					$c++;
					@RefEndC = sort{$b<=>$a}keys%{ $RefCoord{$chr}{$refStarts[$c]} };
				}
				foreach my $bedEnd (sort{$a<=>$b}keys%{ $intervals{$chr}{$bedstart} }) {
					if ($bedEnd < $refStarts[$c]) {				##out of genes
						if ($c == 0) {							##before 1st gene
							foreach my $refEnd (@RefEndC) {
								foreach my $nm (@{ $RefCoord{$chr}{$refStarts[$c]}{$refEnd} }) {
									$intervals{$chr}{$bedstart}{$bedEnd}{"interG"}{$nm} = "--".($refStarts[$c] - $bedEnd);
								}
							}
						} else {								##between 2 genes
							my @RefEnd0 = sort{$b<=>$a}keys%{ $RefCoord{$chr}{$refStarts[($c-1)]} };
							if (($refStarts[$c] - $bedEnd) < ($bedstart - $RefEnd0[0])) {
								foreach my $nm (@{ $RefCoord{$chr}{$refStarts[$c]}{$RefEndC[0]} }) {
									$intervals{$chr}{$bedstart}{$bedEnd}{"interG"}{$nm} = "--".($refStarts[$c] - $bedEnd);
								}
							} else {
								foreach my $nm (@{ $RefCoord{$chr}{$refStarts[($c-1)]}{$RefEnd0[0]} }) {
									$intervals{$chr}{$bedstart}{$bedEnd}{"interG"}{$nm} = "+".($bedstart - $RefEnd0[0]);
								}
							}
						}
					} elsif (($c == $#refStarts) && ($bedstart > $RefEndC[0])) {	##beyond last gene
						foreach my $nm (@{ $RefCoord{$chr}{$refStarts[$c]}{$RefEndC[0]} }) {
							$intervals{$chr}{$bedstart}{$bedEnd}{"interG"}{$nm} = "+".($bedstart - $RefEndC[0]);
						}
					} else {
						my $c2 = $c;
						my @RefEndC2 = @RefEndC;
						while (($c2 < $#refStarts) && ($bedEnd >= $refStarts[$c2])) {
							foreach my $refEnd (@RefEndC2) {
								if ($bedstart <= $refEnd) {
									foreach my $nm (@{ $RefCoord{$chr}{$refStarts[$c2]}{$refEnd} }) {
										complete_intervalBed($chr,$bedstart,$bedEnd,$nm,\%intervals,$RefID);
									}
								} else {
									last;
								}
							}
							$c2++;
							@RefEndC2 = sort{$b<=>$a}keys%{ $RefCoord{$chr}{$refStarts[$c2]} };
						}
						if ($bedEnd >= $refStarts[$c2]) {	##for $c2==$#refStarts
							foreach my $refEnd (@RefEndC2) {
								if ($bedstart <= $refEnd) { 
									foreach my $nm (@{ $RefCoord{$chr}{$refStarts[$c2]}{$refEnd} }) {
										complete_intervalBed($chr,$bedstart,$bedEnd,$nm,\%intervals,$RefID);
									}
								} else {
									last;
								}
							}
						}
					}
				}
			}
		}
	}

	## print annotated bed
	open($fh, ">", ${$reAnnot_r}{"file"}) || die "can't create file ".${$reAnnot_r}{"file"}." ($!)\n";
	foreach my $chr (@{$chromOrder_r}) {
		foreach my $start ( sort{$a<=>$b}(keys%{ $intervals{$chr} }) ) {
			foreach my $end ( sort{$a<=>$b}(keys%{ $intervals{$chr}{$start} }) ) {
				my $line = ${$chromName_r}{"bed"}{$chr}."\t".$start."\t".$end."\t";
				if (exists $intervals{$chr}{$start}{$end}{"interG"}) {
					if (exists ${$reAnnot_r}{"item"}{"o"}) {
						$line .= "interG";
						if (exists ${$reAnnot_r}{"item"}{"g"} || exists ${$reAnnot_r}{"item"}{"t"}) {
							$line .= ":";
							my %genes = ();
							foreach my $nm (keys%{ $intervals{$chr}{$start}{$end}{"interG"} }) {
								push(@{ $genes{${$RefID}{$nm}{"gene"}} }, $nm);
							}
							my @geneTxt = (); my $g = 0;
							foreach my $gene (sort(keys%genes)) {
								if (exists ${$reAnnot_r}{"item"}{"g"}) { $geneTxt[$g] = "$gene"; }
								if (exists ${$reAnnot_r}{"item"}{"g"} && exists ${$reAnnot_r}{"item"}{"t"}) { $geneTxt[$g] .= ":"; }
								if (exists ${$reAnnot_r}{"item"}{"t"}) {
									my @nmTxt = ();
									foreach my $nm (sort@{ $genes{$gene} }) {
										push(@nmTxt , $nm.$intervals{$chr}{$start}{$end}{"interG"}{$nm});
									}
									$geneTxt[$g] .= join(",", @nmTxt);
								}
								$g++;
							}
							$line .=  join(";", @geneTxt);
						}
					} else {
						$line .= ".";
					}
				} elsif (exists $intervals{$chr}{$start}{$end}{"NM"}) {
					my (%genes,%someExon);
					foreach my $nm (keys%{ $intervals{$chr}{$start}{$end}{"NM"} }) {
						push(@{ $genes{${$RefID}{$nm}{"gene"}} }, $nm);
						#if (exists $intervals{$chr}{$start}{$end}{"NM"}{$nm}{"exon"}) { $someExon{${$RefID}{$nm}{"gene"}} = 1; }  # this not to print introns of nms if some other nms with exons, in same gene
					}
					my @geneTxt = (); my $g = 0;
					foreach my $gene (sort(keys%genes)) {
						if (exists ${$reAnnot_r}{"item"}{"g"}) {
							$geneTxt[$g] = "$gene";
							if ( exists ${$reAnnot_r}{"item"}{"t"} || exists ${$reAnnot_r}{"item"}{"i"} || exists ${$reAnnot_r}{"item"}{"e"} ) {
								$geneTxt[$g] .= ":";
							}
						}
						my @nmTxt = (); my $n = 0;
						foreach my $nm (sort@{ $genes{$gene} }) {
							#if (exists $intervals{$chr}{$start}{$end}{"NM"}{$nm}{"exon"} || (exists $intervals{$chr}{$start}{$end}{"NM"}{$nm}{"intron"} && !exists $someExon{$gene}) ) {
								if (exists ${$reAnnot_r}{"item"}{"t"}) {
									$nmTxt[$n] = "$nm";
									if ((exists $intervals{$chr}{$start}{$end}{"NM"}{$nm}{"intron"} && exists ${$reAnnot_r}{"item"}{"i"}) || ($intervals{$chr}{$start}{$end}{"NM"}{$nm}{"exon"} && exists ${$reAnnot_r}{"item"}{"e"}) ) {
										$nmTxt[$n] .= ":";
									}
								}
								if (exists $intervals{$chr}{$start}{$end}{"NM"}{$nm}{"intron"}) {
									if (exists ${$reAnnot_r}{"item"}{"i"}) {		# && !exists $someExon{$gene}) {
										$nmTxt[$n] .= "ivs-".$intervals{$chr}{$start}{$end}{"NM"}{$nm}{"intron"};
									}
								} else {
									if (exists ${$reAnnot_r}{"item"}{"e"}) {
										my @exTxt = (); my$e = 0;
										foreach (@{ $intervals{$chr}{$start}{$end}{"NM"}{$nm}{"exon"} }) {
											$exTxt[$e] = "$_";
											if ($intervals{$chr}{$start}{$end}{"NM"}{$nm}{"NC"}) { $exTxt[$e] .= ".NC"; }
											$e++;
										}
										$nmTxt[$n] .= "exon-".join("-", @exTxt);
									}
								}
								$n++;
							#}
						}
						$geneTxt[$g] .= join(",", @nmTxt);
						$g++;
					}
					$line .= join(";", @geneTxt);
				} else {
					$line .= ".";
				}

				if ($intervals{$chr}{$start}{$end}{"info"}) {
					$line .= $intervals{$chr}{$start}{$end}{"info"};
				}

				print $fh "$line"."\n";
			}
		}
	}
	close ($fh);

}

####
sub complete_intervalBed {

	my ($chr,$bedstart,$bedEnd,$nm,$interval_r,$RefID_r) = @_;
	my $i = 0;
	my $N_exons = scalar@{ ${$RefID_r}{$nm}{"ex_starts"} };
	while ( ($i<($N_exons-1)) && ($bedstart > ${$RefID_r}{$nm}{"ex_ends"}[$i]) ) { $i++; }
	if ($bedEnd < ${$RefID_r}{$nm}{"ex_starts"}[$i]) {
		if (${$RefID_r}{$nm}{"strand"} eq "+") { 
			${$interval_r}{$chr}{$bedstart}{$bedEnd}{"NM"}{$nm}{"intron"} = $i;
		} else { 
			${$interval_r}{$chr}{$bedstart}{$bedEnd}{"NM"}{$nm}{"intron"} = $N_exons-$i;
		}
	} else {
		while ( ($i<($N_exons-1)) && ($bedEnd >= ${$RefID_r}{$nm}{"ex_starts"}[$i]) ) {
			if (${$RefID_r}{$nm}{"strand"} eq "+") {
				push(@{ ${$interval_r}{$chr}{$bedstart}{$bedEnd}{"NM"}{$nm}{"exon"} }, ($i+1));
			} else {
				push(@{ ${$interval_r}{$chr}{$bedstart}{$bedEnd}{"NM"}{$nm}{"exon"} }, ($N_exons-$i));
			}
			if (exists ${$RefID_r}{$nm}{"CDS_start"}) {
				if ($bedstart >= ${$RefID_r}{$nm}{"CDS_end"}) {
					${$interval_r}{$chr}{$bedstart}{$bedEnd}{"NM"}{$nm}{"NC"} = 1;
				}
				if ($bedEnd <= ${$RefID_r}{$nm}{"CDS_start"}) {
					${$interval_r}{$chr}{$bedstart}{$bedEnd}{"NM"}{$nm}{"NC"} = 1;
				}
			} else {
				${$interval_r}{$chr}{$bedstart}{$bedEnd}{"NM"}{$nm}{"NC"} = 1;
			}
			$i++
		}
		if ($bedEnd >= ${$RefID_r}{$nm}{"ex_starts"}[$i]) {
			if (${$RefID_r}{$nm}{"strand"} eq "+") {
				push(@{ ${$interval_r}{$chr}{$bedstart}{$bedEnd}{"NM"}{$nm}{"exon"} }, ($i+1));
			} else {
				push(@{ ${$interval_r}{$chr}{$bedstart}{$bedEnd}{"NM"}{$nm}{"exon"} }, ($N_exons-$i));
			}
			if (exists ${$RefID_r}{$nm}{"CDS_start"}) {
				if ($bedstart >= ${$RefID_r}{$nm}{"CDS_end"}) {
					${$interval_r}{$chr}{$bedstart}{$bedEnd}{"NM"}{$nm}{"NC"} = 1;
				}
				if ($bedEnd <= ${$RefID_r}{$nm}{"CDS_start"}) {
					${$interval_r}{$chr}{$bedstart}{$bedEnd}{"NM"}{$nm}{"NC"} = 1;
				}
			} else {
				${$interval_r}{$chr}{$bedstart}{$bedEnd}{"NM"}{$nm}{"NC"} = 1;
			}
		}
	}

}


######################
sub addLength2Bed {

	my ($bedFile,$newBed,$len) = @_;
	open(my $fhIn, "<", $bedFile) || die "\n!! can't open file $bedFile\n";
	open(my $fhOut, ">", $newBed) || die "\n!! can't create file $newBed\n";
	while (my $line = <$fhIn>) {
		if ( ($line !~ /^\w+\t\d+\t\d+/) || ($line =~ /^#/) ) {
			print $fhOut "$line";
		} else {
			chomp $line;
			my @tab = split(/\t/,$line);
			$tab[1] = $tab[1] - $len;
			$tab[2] = $tab[2] + $len;
			print $fhOut join("\t",@tab)."\n";
			}
		}
	close($fhIn);
	close($fhOut);

}


######################
sub sliceBedIntervals {

	my ($bedFile, $outFile, $splitBed, $rmOverlapBed, $cutBed, $cutB_opt_r, $chromOrder_r, $chromName_r) = @_;
	my $currentBed = $bedFile;
	if ($splitBed) {
		print "\n\tsplitting overlapping intervals from $currentBed\n";
		$splitBed = "$outFile";
		$splitBed =~ s/.bed$//;
		$splitBed .= "_split.bed";
		split_Bed($currentBed, $splitBed, $chromOrder_r, $chromName_r);
		$currentBed = $splitBed;
	} elsif ($rmOverlapBed) {
		print "\n\tremoving overlapping intervals from $currentBed\n";
		$rmOverlapBed = "$outFile";
		$rmOverlapBed =~ s/.bed$//;
		$rmOverlapBed .= "_flat.bed";
		rm_overlap_Bed($currentBed, $rmOverlapBed, $chromOrder_r, $chromName_r);
		$currentBed = $rmOverlapBed;
	}
	if ($cutBed) {
		print "\n\tcutting intervals from $currentBed 
			(above ".${$cutB_opt_r}{"maxL"}." bp ; in ".${$cutB_opt_r}{"cutL"}." bp pieces , if longer than ".${$cutB_opt_r}{"minL"}." bp , else ";
		if (${$cutB_opt_r}{"keepLast"} eq "merge") {
			print "merge 2 last)\n";
		} elsif (${$cutB_opt_r}{"keepLast"} eq "half") {
			print "keep 2 last with half of their sum)\n";
		} else {
			print "through it)\n";
		}
		$cutBed = "$outFile";
		$cutBed =~ s/.bed$//;
		$cutBed .= "_cut".${$cutB_opt_r}{"cutL"}.".bed";
		cut_Bed($currentBed, $cutBed, $cutB_opt_r);
		unless ($currentBed eq $bedFile) { unlink $currentBed; }
		$currentBed = $cutBed;
	}
	return($currentBed);

}


######################
sub split_Bed {

	my ($bedFile,$newBed,$chromOrder_r,$chromName_r) = @_;

	##bedFile in %intervals
	my %intervals;
	open(my $fh, "<", $bedFile) || die "\n!! can't open file $bedFile\n";
	while (my $line = <$fh>) {
		if ( ($line =~ /^\w+\t\d+\t\d+/) && ($line !~ /^#/) ) {
			chomp $line;
			my @tab = split(/\t/,$line);
			my $chr = $tab[0];
			$chr =~ s/^chr//i;
			$intervals{$chr}{$tab[1]}{$tab[2]} = "";
			if ($tab[3]) {
				for my $i (3..$#tab) {
					$intervals{$chr}{$tab[1]}{$tab[2]} .= "\t".$tab[$i];
				}
			}
			#else { $intervals{$chr}{$tab[1]}{$tab[2]} = "\t."; }
		}
	}
	close ($fh);
	unless (%intervals) { die "\n!! no bed formatted line in $bedFile file"; }

	#split intervals with same start
	# $interval2{$chr}[idx]{"start"/"end"/"info"} = ;
	my %interval2;
	foreach my $chr (keys%intervals) {
		my $i = 0;
		foreach my $start ( sort{$a<=>$b}(keys%{ $intervals{$chr} }) ) {
			my $currentStart=$start;
			foreach my $end ( sort{$a<=>$b}(keys%{ $intervals{$chr}{$start} }) ) {
				$interval2{$chr}[$i]{"start"} = $currentStart;
				$interval2{$chr}[$i]{"end"} = $end;
				$interval2{$chr}[$i]{"info"} = $intervals{$chr}{$start}{$end};
				$currentStart=$end;
				$i++;
			}
		}
	}
	%intervals = ();
	#split overlapping intervals
	# $interval3{$chr}{$start}{"end"/"info"} = ;
	my %interval3;
	foreach my $chr (keys%interval2) { 
		my $i = 0;
		while (($i+1) < scalar@{ $interval2{$chr} } ) {
			if ($interval2{$chr}[$i+1]{"start"} < $interval2{$chr}[$i]{"end"}) {
				my $End = $interval2{$chr}[$i]{"end"};
				my $Info = $interval2{$chr}[$i]{"info"};
				my @allPos = ($interval2{$chr}[$i]{"start"},$interval2{$chr}[$i]{"end"}); 
				my $j = $i;
				while (($j+1) < scalar@{ $interval2{$chr} } && $interval2{$chr}[$j+1]{"start"} < $End) {
					push(@allPos,$interval2{$chr}[$j+1]{"start"});
					if ($interval2{$chr}[$j+1]{"end"} < $End) {
						push(@allPos,$interval2{$chr}[$j+1]{"end"});
					} else {
						if ($interval2{$chr}[$j+1]{"end"} > $End) {
							$interval2{$chr}[$j+1]{"start"} = $End;
						}
					}
					$j++;
				}
				my @allPos2 = sort{$a<=>$b}@allPos;
				for my $k (0..($#allPos2-1)) {
					$interval3{$chr}{$allPos2[$k]}{"end"} = $allPos2[$k+1];
					$interval3{$chr}{$allPos2[$k]}{"info"} = $Info;
				}
			} else {
				unless (exists $interval3{$chr}{$interval2{$chr}[$i]{"start"}}{"end"}) {
					$interval3{$chr}{$interval2{$chr}[$i]{"start"}}{"end"} = $interval2{$chr}[$i]{"end"};
					$interval3{$chr}{$interval2{$chr}[$i]{"start"}}{"info"} = $interval2{$chr}[$i]{"info"};
				}
			}
			$i++;
		}
		unless (exists $interval3{$chr}{$interval2{$chr}[$i]{"start"}}{"end"}) {
			$interval3{$chr}{$interval2{$chr}[$i]{"start"}}{"end"} = $interval2{$chr}[$i]{"end"};
			$interval3{$chr}{$interval2{$chr}[$i]{"start"}}{"info"} = $interval2{$chr}[$i]{"info"};
		}
	}
	%interval2 = ();
	#print splitted bed
	open($fh, ">", "$newBed") || die "\n!! can't create file $newBed\n";
	foreach my $chr (@{$chromOrder_r}) {
		if (exists $interval3{$chr}) {
			foreach my $start ( sort{$a<=>$b}(keys%{ $interval3{$chr} }) ) { 
				print $fh ${$chromName_r}{"bed"}{$chr} ."\t".$start."\t".$interval3{$chr}{$start}{"end"}.$interval3{$chr}{$start}{"info"}."\n";
			}
		}
	}
	close($fh);

}


######################
sub rm_overlap_Bed {

	my ($bedFile,$newBed,$chromOrder_r,$chromName_r) = @_;

	##bedFile in %intervals
	my %intervals;
	open(my $fh, "<", $bedFile) || die "\n!! can't open file $bedFile\n";
	while (my $line = <$fh>) {
		if ( ($line =~ /^\w+\t\d+\t\d+/) && ($line !~ /^#/) ) {
			chomp $line;
			my @tab = split(/\t/,$line);
			my $chr = $tab[0];
			$chr =~ s/^chr//i;
			$intervals{$chr}{$tab[1]}{$tab[2]} = "";
			if ($tab[3]) {
				for my $i (3..$#tab) {
					$intervals{$chr}{$tab[1]}{$tab[2]} .= "\t".$tab[$i];
				}
			}
			#else { $intervals{$chr}{$tab[1]}{$tab[2]} = "\t."; }
		}
	}
	close ($fh);
	unless (%intervals) { die "\n!! no bed formatted line in $bedFile file"; }

	#intervals with same start : start = end of second to last; end = end of last
	# $interval2{$chr}[idx]{"start"/"end"/"info"} = ;
	my %interval2;
	foreach my $chr (keys%intervals) {
		my $i = 0;
		foreach my $start ( sort{$a<=>$b}(keys%{ $intervals{$chr} }) ) {
			my @Ends = sort{$a<=>$b}(keys%{ $intervals{$chr}{$start} });
			if (scalar@Ends > 1) {
				$interval2{$chr}[$i]{"start"} = $Ends[-2];
				$interval2{$chr}[$i]{"end"} = $Ends[-1];
				$interval2{$chr}[$i]{"info"} = $intervals{$chr}{$start}{$Ends[-1]};
				$i++;
			} else {
				$interval2{$chr}[$i]{"start"} = $start;
				$interval2{$chr}[$i]{"end"} = $Ends[0];
				$interval2{$chr}[$i]{"info"} = $intervals{$chr}{$start}{$Ends[0]};
				$i++;
			}
		}
	}
	%intervals = ();
	#overlapping intervals : 
	# $interval3{$chr}{$start}{"end"/"info"} = ;
	my %interval3;
	foreach my $chr (keys%interval2) { 
		my $i = 0;
		while ( ($i+1) < scalar@{ $interval2{$chr} } ) {
			if ($interval2{$chr}[$i+1]{"start"} < $interval2{$chr}[$i]{"end"}) {
				if ($interval2{$chr}[$i+1]{"end"} > $interval2{$chr}[$i]{"start"}) {
					if ($interval2{$chr}[$i+1]{"start"} > $interval2{$chr}[$i]{"start"}) {
						$interval3{$chr}{$interval2{$chr}[$i]{"start"}}{"end"} = $interval2{$chr}[$i+1]{"start"};
						$interval3{$chr}{$interval2{$chr}[$i]{"start"}}{"info"} = $interval2{$chr}[$i]{"info"};
					}
					if ($interval2{$chr}[$i+1]{"end"} < $interval2{$chr}[$i]{"end"}) {
						$interval2{$chr}[$i+1]{"start"} = $interval2{$chr}[$i+1]{"end"};
						$interval2{$chr}[$i+1]{"end"} = $interval2{$chr}[$i]{"end"};
					} else {
						$interval2{$chr}[$i+1]{"start"} = $interval2{$chr}[$i]{"end"};
					}
				} else {
					$interval2{$chr}[$i+1]{"start"} = $interval2{$chr}[$i]{"start"};
					$interval2{$chr}[$i+1]{"end"} = $interval2{$chr}[$i]{"end"};
				}
			} else {
				if (!exists $interval3{$chr}{$interval2{$chr}[$i]{"start"}}{"end"} && $interval2{$chr}[$i]{"start"} != $interval2{$chr}[$i]{"end"}) {
					$interval3{$chr}{$interval2{$chr}[$i]{"start"}}{"end"} = $interval2{$chr}[$i]{"end"};
					$interval3{$chr}{$interval2{$chr}[$i]{"start"}}{"info"} = $interval2{$chr}[$i]{"info"};
				}
			}
			$i++;
		}
		if (!exists $interval3{$chr}{$interval2{$chr}[$i]{"start"}}{"end"} && $interval2{$chr}[$i]{"start"} != $interval2{$chr}[$i]{"end"}) {
			$interval3{$chr}{$interval2{$chr}[$i]{"start"}}{"end"} = $interval2{$chr}[$i]{"end"};
			$interval3{$chr}{$interval2{$chr}[$i]{"start"}}{"info"} = $interval2{$chr}[$i]{"info"};
		}
	}
	%interval2 = ();
	#print no_overlap bed
	open($fh, ">", $newBed) || die "\n!! can't create file $newBed\n";
	foreach my $chr (@{$chromOrder_r}) {
		if (exists $interval3{$chr}) {
			foreach my $start ( sort{$a<=>$b}(keys%{ $interval3{$chr} }) ) { 
				print $fh ${$chromName_r}{"bed"}{$chr}."\t".$start."\t".$interval3{$chr}{$start}{"end"}.$interval3{$chr}{$start}{"info"}."\n";
			}
		}
	}
	close($fh);

}



######################
sub cut_Bed {

	my ($bedFile, $newBed, $c_opt) = @_;
	my$cutL = ${$c_opt}{"cutL"};

	open(my $fhOut, ">", "$newBed") || die "\n!! can't create file $newBed\n";
	open(my $fhBed, "<", "$bedFile") || die "\n!! can't open file $bedFile\n";
	while (my $line = <$fhBed>) {
		chomp $line;
		my @tab = split(/\t/,$line);
		my $line2;
		if ( ($tab[2]-$tab[1]) >= ${$c_opt}{"maxL"} ) {
			my $start = $tab[1];
			my $end = $tab[2];
			my $div = int(($end-$start) / $cutL);
			my $mod = ($end-$start) % $cutL;
			my $info = "";
			for (my $i = 3; $i < scalar@tab; $i++) { $info .= $tab[$i]."\t"; }
			chop $info;
			my$i = 0;
			if (${$c_opt}{"keepLast"} && ($mod<${$c_opt}{"minL"})) {
				while ($i < ($div-1)) {
					$line2 = $tab[0]."\t".($start+($i*$cutL))."\t".($start+(($i+1)*$cutL))."\t";
					if ($info) { $line2 .= $info; }
					else { chop $line2; }
					print $fhOut $line2."\n";
					$i++;
				}
				if (${$c_opt}{"keepLast"} eq "merge") {
					$line2 = $tab[0]."\t".($start+($i*$cutL))."\t".$end."\t";
					if ($info) { $line2 .= $info; }
					else { chop $line2; }
				} else {
					my $lastI = $end-($start+($i*$cutL));
					$line2 = $tab[0]."\t".($start+($i*$cutL))."\t".($start+($i*$cutL)+int($lastI/2))."\t";
					if ($info) { $line2 .= $info; }
					else { chop $line2; }
					print $fhOut $line2."\n";
					$line2 = $tab[0]."\t".($start+($i*$cutL)+int($lastI/2))."\t".$end."\t";
					if ($info) { $line2 .= $info; }
					else { chop $line2; }
				}
				print $fhOut $line2."\n";
			} else {
				while ($i < $div) {
					$line2 = $tab[0]."\t".($start+($i*$cutL))."\t".($start+(($i+1)*$cutL))."\t";
					if ($info) { $line2 .= $info; }
					else { chop $line2; }
					print $fhOut $line2."\n";
					$i++;
				}
				if ($mod >= ${$c_opt}{"minL"}) {
					$line2 = $tab[0]."\t".($start+($i*$cutL))."\t".$end."\t";
					if ($info) { $line2 .= $info; }
					else { chop $line2; }
					print $fhOut $line2."\n";
				}
			}
		} else {
			foreach my $i (@tab) { $line2 .= $i."\t"; }
			chop $line2;
			print $fhOut $line2."\n";
		}
	}
	close($fhBed);
	close($fhOut);

}
#########################
#printBed("$tmpDir/original",\%Bed,\@ChromOrder);

sub printBed {

	my ($bedFile, $interval_r, $chromOrder_r) = @_;
	open (my $fh0, ">", "$bedFile\_0Chr.bed")  || die "\n!! cannot create $bedFile\_0Chr.bed\n";
	open (my $fhW, ">", "$bedFile\_wChr.bed")  || die "\n!! cannot create $bedFile\_wChr.bed\n";
	if (@{$chromOrder_r}) {
		foreach my $chr (@{$chromOrder_r}) {
			if (exists ${$interval_r}{$chr}) {
				foreach my $start (sort{$a<=>$b}keys%{ ${$interval_r}{$chr} }) {
					print $fh0 $chr."\t".($start-1)."\t".${$interval_r}{$chr}{$start}."\n";
					print $fhW "chr".$chr."\t".($start-1)."\t".${$interval_r}{$chr}{$start}."\n";
				}
			}
		}
	} else {
		foreach my $chr (sort(keys%{$interval_r})) {
			foreach my $start (sort{$a<=>$b}keys%{ ${$interval_r}{$chr} }) {
				print $fh0 $chr."\t".($start-1)."\t".${$interval_r}{$chr}{$start}."\n";
				print $fhW "chr".$chr."\t".($start-1)."\t".${$interval_r}{$chr}{$start}."\n";
			}
		}
	}
	close $fh0;
	close $fhW;

}

#########################


1;


