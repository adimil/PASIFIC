#!/bin/env perl
use strict; 
use Bio::SeqIO;

# Set params
my ($tmpdir, $name, $p_P1, $p_simpAT, $p_limitAT, $p_findterm) = @ARGV;
my @TMminenergy;
my @ATminenergy;
my $stem3_len;
my $loop34_len;
my $stem4_len;
my $allavg = 0;
my $SLavg = 0;
my $allmul = 1;
my $SLmul = 1;
my $ID;
my $seq;
my $seq_len;

# Open the in fasta file, and the output file
my $seq_file = Bio::SeqIO -> new(-file => $tmpdir.$name."\.fasta", 
								-format => 'Fasta' );
open(my $fh_out,">",$tmpdir.$name.".out") or die "can't open file out: $!";

# Go through sequences
while (my $seq_obj = $seq_file->next_seq())
{
	# get ID and sequence
	$ID = $seq_obj->display_id;
	$seq = $seq_obj->seq;
	$seq =~ tr/acgutT/ACGUUU/;
	
	# If 3' is unknown
	if ($p_findterm)
	{
		#######################################
		#### find terminator - 3' unknown #####
		#######################################		
		@TMminenergy = find_term_noseq($seq);
	}
	else
	{
		# remove polyU
		if ($seq !~ /(?<polyU>(U{3,}(.{0,2}U+){0,3}.{0,5}$))/)
		{
			print $fh_out "$ID\t$seq\tCould not find terminator" . "\t"x54 . "\n";
			next;
		}
		else
		{
			my $noUseq = substr($seq, 0, $-[0]);

			########################################
			#### find terminator - TTS is known ####
			########################################
			@TMminenergy = find_term($noUseq);
		}
	}

	# Check that found a terminator
	if (!@TMminenergy)
	{
		print $fh_out "$ID\t$seq\tcould not find terminator" . "\t"x54;
	}
	else
	{
		$seq_len = length($seq);
		$allavg = 0;
		$SLavg = 0;
		$allmul = 1;
		$SLmul = 1;

		#######################################
		#### look for antiterminator helix ####
		#######################################
		# get lengths of stems
		$TMminenergy[0] =~ /(?<stem3>(\(+(\.*\(+)*))(?<loop34>(\.+))(?<stem4>(\)+(\.*\)+)*))/;
		($stem3_len,$loop34_len,$stem4_len) = (length($+{stem3}),length($+{loop34}),length($+{stem4}));
		my $no4seq = substr($seq, 0, $TMminenergy[4] - $stem4_len +1);

		@ATminenergy = find_antiterminator($no4seq);

		# Check that found possible AT
		if (!@ATminenergy)
		{
			print $fh_out "$ID\t$seq\tNo antiterminator was found" . "\t"x54;
		}
		else
		{
			#################
			#### Find P1 ####
			#################
			my @longestP1;
			
			if($p_P1)
			{
				# Sequnce to look for P1 in - 
				$ATminenergy[0] =~ /\(\.+\)/;
				my $P1seq = substr($seq, 0, $ATminenergy[2]+$-[0]);

				@longestP1 = find_P1($P1seq);
			}

			my $TMconstraint;

			# set constraints of TM fold
			if (!@longestP1)
			{
				@longestP1 = ("No P1 was found",0,0,0,0);
				$TMconstraint = "." x ($TMminenergy[2]-1) . $TMminenergy[0] . "x" x ($seq_len-$TMminenergy[4]);
			}
			else
			{
				$TMconstraint = "." x ($longestP1[2]-1) . $longestP1[0] . "." x ($TMminenergy[2]-$longestP1[4]-1) . $TMminenergy[0] . "x" x ($seq_len-$TMminenergy[4]);
			}
			
			# Set constraints of AT fold
			my $ATconstraint = "." x ($ATminenergy[2]-1) . $ATminenergy[0] . "x" x ($seq_len-$ATminenergy[4]);

			##############################
			#### fold the two options ####
			##############################
			my $folddesc = ">".$tmpdir."TM_".$name."_".$ID;
			my @TMfold = `printf "$folddesc\n$seq\n$TMconstraint" | RNAfold -p -d2 --noLP -C`;
			$folddesc = ">".$tmpdir."AT_".$name."_".$ID;
			my @ATfold = `printf "$folddesc\n$seq\n$ATconstraint" | RNAfold -p -d2 --noLP -C`;
			
			
			################################
			#### Start printing results ####
			################################
			print $fh_out "$ID\t$seq\t";
			
			############################
			# Print terminator details #
			############################
			my $TMseq = substr($seq, $TMminenergy[2]-1, $TMminenergy[3]);

			# terminator,terminator_seq, start, len, end, G34
			print $fh_out $TMminenergy[0];
			print $fh_out "\t";
			print $fh_out $TMseq;
			print $fh_out "\t";
			print $fh_out $TMminenergy[2]; # irrelevant for current classifier
			print $fh_out "\t";
			print $fh_out $TMminenergy[3]; # irrelevant for current classifier
			print $fh_out "\t";
			print $fh_out $TMminenergy[4]; # irrelevant for current classifier
			print $fh_out "\t";
			print $fh_out $TMminenergy[1]; # irrelevant for current classifier
			print $fh_out "\t";
			
			# G34/len
			my $G34_len = $TMminenergy[1]/$TMminenergy[3]; # irrelevant for current classifier
			print $fh_out $G34_len . "\t";
			
			# Gfold34, Gfold34/len
			$TMfold[2]=~ m/(?<GTMfold>(\-\d+(\.\d+)?))/;
			my $GTMfold = $+{GTMfold};
			print $fh_out $GTMfold . "\t";
			my $GTMfold_len = $GTMfold/$seq_len;
			print $fh_out $GTMfold_len . "\t"; # irrelevant for current classifier

			# frequency of mfe structure in ensemble, ensemble diversity
			$TMfold[5]=~ m/(?<FRQ_TM>(\d+(\.\d+)?(e\-\d+)?))\;/;
			print $fh_out $+{FRQ_TM} . "\t"; # irrelevant for current classifier
			$TMfold[5]=~ m/\;.+?(?<DIV_TM>(\d+(\.\d+)?))/;
			print $fh_out $+{DIV_TM} . "\t"; # irrelevant for current classifier

			# Average prob
			avgProb("TM", $TMminenergy[2], $TMminenergy[4]);
			print $fh_out $allavg . "\t"; # irrelevant for current classifier
			print $fh_out $SLavg . "\t"; # irrelevant for current classifier
			print $fh_out $allmul . "\t"; # irrelevant for current classifier
			print $fh_out $SLmul . "\t"; # irrelevant for current classifier
			
			my $TMallavg = $allavg;
			my $TMSLavg = $SLavg;
			my $TMallmul = $allmul;
			my $TMSLmul = $SLmul;
			
			#################################
			# Print anti-terminator details #
			#################################
			my $ATseq = substr($seq, $ATminenergy[2]-1, $ATminenergy[3]);

			# antiterminator,ATseq ,start, len, end, G23
			print $fh_out $ATminenergy[0];
			print $fh_out "\t";
			print $fh_out $ATseq;
			print $fh_out "\t";
			print $fh_out $ATminenergy[2]; # irrelevant for current classifier
			print $fh_out "\t";
			print $fh_out $ATminenergy[3]; # irrelevant for current classifier
			print $fh_out "\t";
			print $fh_out $ATminenergy[4];
			print $fh_out "\t";
			print $fh_out $ATminenergy[1];
			print $fh_out "\t";
			
			# G23/len
			print $fh_out $ATminenergy[1]/$ATminenergy[3];
			print $fh_out "\t";
			
			# Gfold23, Gfold23/len
			$ATfold[2]=~ m/(?<GATfold>(\-\d+(\.\d+)?))/;
			my $GATfold = $+{GATfold};
			print $fh_out "$+{GATfold}\t";
			print $fh_out $+{GATfold}/$seq_len;
			print $fh_out "\t";

			# frequency of mfe structure in ensemble, ensemble diversity
			$ATfold[5]=~ m/(?<FRQ_AT>(\d+(\.\d+)?(e\-\d+)?))\;/;
			print $fh_out "$+{FRQ_AT}\t"; # irrelevant for current classifier
			$ATfold[5]=~ m/\;.+?(?<DIV_AT>(\d+(\.\d+)?))/;
			print $fh_out "$+{DIV_AT}\t"; # irrelevant for current classifier
			
			# Average prob
			avgProb("AT", $ATminenergy[2], $ATminenergy[4]);
			print $fh_out "$allavg\t$SLavg\t$allmul\t$SLmul"; # irrelevant for current classifier
			print $fh_out "\t";
						
			####################
			# Print P1 details #
			####################
			my $P1seq = substr($seq, $longestP1[2]-1, $longestP1[3]);
			if ($P1seq eq "")
			{
				$P1seq = "N/A";
			}

			# P1,P1seq ,start, len, end, GP1, GP1/len
			print $fh_out $longestP1[0];
			print $fh_out "\t";
			print $fh_out $P1seq;
			print $fh_out "\t";
			print $fh_out $longestP1[2]; # irrelevant for current classifier
			print $fh_out "\t";
			print $fh_out $longestP1[3];
			print $fh_out "\t";
			print $fh_out $longestP1[4];  # irrelevant for current classifier
			print $fh_out "\t";
			print $fh_out $longestP1[1];
			print $fh_out "\t";
			print $fh_out $longestP1[1]/$seq_len;
			print $fh_out "\t";

			#################
			# Print other:  #
			#################
			
			# overlap
			print $fh_out $ATminenergy[4]-$TMminenergy[2]-1; # irrelevant for current classifier
			print $fh_out "\t";
			
			# ddG of TM and AT
			print $fh_out $GATfold-$GTMfold; # irrelevant for current classifier
			print $fh_out "\t";
			print $fh_out $GATfold/$GTMfold;
			print $fh_out "\t";
			print $fh_out $ATminenergy[1]-$TMminenergy[1]; # irrelevant for current classifier
			print $fh_out "\t";
			print $fh_out $ATminenergy[1]/$TMminenergy[1];
			print $fh_out "\t";

			###########################
			# print naive fold params #
			###########################
			my @naivefold = `echo $seq | RNAfold -p -d2 --noLP`;
			$naivefold[1]=~ m/(?<Gnaive>(\-?\d+(\.\d+)?))/;
			my $Gnaive = $+{Gnaive};
			print $fh_out "$+{Gnaive}\t"; # irrelevant for current classifier
			$naivefold[4]=~ m/(?<FRQnaive>(\d+(\.\d+)?(e\-\d+)?))\;/;
			print $fh_out "$+{FRQnaive}\t"; # irrelevant for current classifier
			$naivefold[4]=~ m/\;.+?(?<DIVnaive>(\d+(\.\d+)?))/;
			print $fh_out "$+{DIVnaive}\t"; # irrelevant for current classifier

			# naive fold dG
			print $fh_out $GATfold-$Gnaive;
			print $fh_out "\t";
			print $fh_out $GATfold/$Gnaive;
			print $fh_out "\t";
			print $fh_out $GTMfold-$Gnaive;
			print $fh_out "\t";
			print $fh_out $GTMfold/$Gnaive;
			print $fh_out "\t";
			
			# probability ratios
			print $fh_out $allavg/$TMallavg; # irrelevant for current classifier
			print $fh_out "\t";
			print $fh_out $SLavg/$TMSLavg; # irrelevant for current classifier
			print $fh_out "\t";
			print $fh_out $allmul/$TMallmul; # irrelevant for current classifier
			print $fh_out "\t";
			print $fh_out $SLmul/$TMSLmul; # irrelevant for current classifier
			print $fh_out "\t";
			
			#######################################
			# print folds in dot-bracket notation #
			#######################################
			print $fh_out ($TMfold[2] =~ /([\(\)\.]+)/);
			print $fh_out "\t";
			print $fh_out ($ATfold[2] =~ /([\(\)\.]+)/);
			
			
			####################################################################
			# delete ps files - uncomment if do not need visual representation #
			####################################################################
			# `rm -f $tmpdir."TM_".$name."_".$ID."_ss.ps"`;
			# `rm -f $tmpdir."TM_".$name."_".$ID."_dp.ps"`;
			# `rm -f $tmpdir."AT_".$name."_".$ID."_ss.ps"`;
			# `rm -f $tmpdir."AT_".$name."_".$ID."_dp.ps"`;
			####################################################################
		}
	}
	print $fh_out "\n";
}
close($fh_out);

####################################################################
# find_term                                                        #
####################################################################
# Finds a terminator when the TTS is known.
#
# Conditions for terminator:
# - Simple helix
# - Helix end is in the end of seq (allow 2-3 spaces)
# - Helix length is less than 45 but longer than 15
#
# Recieves the sequence without the polyU.
#
####################################################################
sub find_term
{
	my ($noUseq) = @_;
	
	# get all local folding
	my @helices = `echo $noUseq | RNALfold -d2 --noLP -L45`;
	chomp(@helices);

	# set all possible terminators in matrix with columns: fold, energy, start, length, end
	my @resultsMTX;

	for (my $hlxline=0; $hlxline<scalar(@helices)-2; ++$hlxline)
	{
		my @helix = (split(/\s+/, $helices[$hlxline]));

		# take only simple stem-loop
		if ($helix[0] !~ /\).*\(/)
		{
			my $energy=0;
			my $start=0;
			
			if (scalar(@helix) == 3)
			{
			$energy = ($helix[1] =~ /(\-?\d+\.\d+)/)[0];
			$start = $helix[2];
			}	
			if (scalar(@helix) == 4)
			{
			$energy = ($helix[2] =~ /(\-?\d+\.\d+)/)[0];
			$start = $helix[3];
			}
			
			my $len = length($helix[0]);
			my $end = $start + $len - 1;
			
			# Check if terminator:
			
			# end is in the end of seq (allow 2-3 spaces)
			if ($end >= length($noUseq)-3)
			{
				# len is less than 45 but longer than 15
				if ($len >= 15)
				{
					# fold, energy, start, length, end
					my @term = ($helix[0], $energy, $start, $len, $end);
					push(@resultsMTX, \@term);
				}
			}
		}
	}

	if (@resultsMTX)
	{
		return @{(sort { $a->[1] <=> $b->[1] } @resultsMTX)[0]};
	}
	else
	{
		return;
	}
}

####################################################################
# find_term_noseq                                                  #
####################################################################
# Finds a terminator when the TTS is unknown.
#
# Conditions for terminator:
# - Simple helix
# - Helix starts after first 25 bases
# - Helix length is less than 45 but longer than 15
# - Helix ends up to 9 bases before end of sequence
# - loop length is 3-10
# - PolyU:
#	* spacer - 0-1 bases
#	* proximal part - 5 bases, at least 3 U
#	* distal part - 4 bases, no 4C or 4pur
#
# Recieves the whole sequence.
#
####################################################################
sub find_term_noseq
{
	my ($tseq) = @_;
	# get all local folding
	my @helices = `echo $tseq | RNALfold -d2 --noLP -L45`;
	chomp(@helices);

	# set all possible terminator in matrix with columns: fold, energy, start, length, end
	my @resultsMTX;

	for (my $hlxline=0; $hlxline<scalar(@helices)-2; ++$hlxline)
	{
		my @helix = (split(/\s+/, $helices[$hlxline]));
		
		# take only simple stem-loop
		if ($helix[0] !~ /\).*\(/)
		{
			my $energy=0;
			my $start=0;

			if (scalar(@helix) == 3)
			{
			$energy = ($helix[1] =~ /(\-?\d+\.\d+)/)[0];
			$start = $helix[2];
			}	
			if (scalar(@helix) == 4)
			{
			$energy = ($helix[2] =~ /(\-?\d+\.\d+)/)[0];
			$start = $helix[3];
			}
			
			# Check if terminator:
			
			# start after 25 (if min length is 80 and max term lengt is 55)
			if ($start >= 25)
			{
				my $len = length($helix[0]);
				
				# len > 15
				if ($len >= 15)
				{
					my $end = $start + $len - 1;
					
					# end leaves minimum 9b for polyU
					if ($end <= length($tseq)-9)
					{
						$helix[0] =~ /\(\.+\)/;
						my $looplen = ($+[0]-1) - ($-[0]+1);
						
						# loop len is 3-10
						if ($looplen >= 3 && $looplen <= 10)
						{
							# check poly U
							my $polyU = substr($tseq, $end-1, 11);
							
							# 0-1 spacer, 5 prox - at least 3 U, 4 dist - no 4C or 4pur
							if ($polyU =~ /^.{0,2}(?=(U.?){3,}).{5}(?!(C{4}|[AG]{4})).{4}/)
							{
								# and overall 4 U
								if ($polyU =~ /(U.*){4,}/)
								{
									# fold, energy, start, length, end
									my @term = ($helix[0], $energy, $start, $len, $end);
									push(@resultsMTX, \@term);
								}
							}
						}
					}
				}
			}
		}
	}
	
	if (@resultsMTX)
	{
		my @chosenterm = @{(sort { $a->[1] <=> $b->[1] } @resultsMTX)[0]};
		$seq = substr($seq, 0, $chosenterm[4] + 11);
		return @chosenterm;
	}
	else
	{
		return;
	}
}

####################################################################
# find_antiterminator                                              #
####################################################################
# Looks for an antiterminator.
#
# Conditions for antiterminator:
# - Helix starts before stem 3
# - Helix ends few bases after stem 3 starts
# - Helix length is more than 20
# - if p_limitAT is checked - distance between strand #2 and #3 should be 55nt or shorter
# - if p_simpAT is checked - simple helix
#
# Recieves the sequence up to stem 4.
#
####################################################################
sub find_antiterminator
{
	my ($no4seq) = @_;
	my $sizelim = 55+2*$stem3_len;
	my $simp = "\\)\\.*\\(";
	if (!$p_limitAT)
	{
		$sizelim = 1000;
	}
	if (!$p_simpAT)
	{
		$simp = "Q";
	}

	# get all local folding
	my @helices = `echo $no4seq | RNALfold -d2 --noLP -L$sizelim`;
	chomp(@helices);

	# set all possible terminators in matrix with columns: fold, energy, start, length, end
	my @resultsMTX;

	for (my $hlxline=0; $hlxline<scalar(@helices)-2; ++$hlxline)
	{
		my @helix = (split(/\s+/, $helices[$hlxline]));

		my $energy=0;
		my $start=0;
		
		if (scalar(@helix) == 3)
		{
		$energy = ($helix[1] =~ /(\-?\d+\.\d+)/)[0];
		$start = $helix[2];
		}	
		if (scalar(@helix) == 4)
		{
		$energy = ($helix[2] =~ /(\-?\d+\.\d+)/)[0];
		$start = $helix[3];
		}
		
		# starts before 3
		if ($start < $TMminenergy[2])
		{
			my $len = length($helix[0]);
			my $end = $start + $len - 1;
			
			# ends few bases after 3 starts
			if ($end > $TMminenergy[2]+5)
			{
				# longer than 20
				if ($len > 20)
				{
					# if checked - simple
					if ($helix[0] !~ /$simp/)
					{
						# fold, energy, start, length, end
						my @antiterm = ($helix[0], $energy, $start, $len, $end, $energy/$len);
						push(@resultsMTX, \@antiterm);
					}
				}
				
			}
		}
	}

	if (@resultsMTX)
	{
		return @{(sort { $a->[1] <=> $b->[1] } @resultsMTX)[0]};
	}
	else
	{
		return;
	}
}

####################################################################
# find_P1                                                          #
####################################################################
# Looks for an anti-antiterminator (P1).
#
# Conditions for P1:
# - Helix starts before the antiterminator
# - Helix ends after antiterminator starts and up to the terminator
#
# Recieves the sequence up to stem 3.
#
####################################################################
sub find_P1
{
	my ($P1seq) = @_;
	my $P1seq_len = length($P1seq);
	
	# get all local folding
	my @helices = `echo $P1seq | RNALfold -d2 --noLP -L$seq_len`;
	chomp(@helices);

	# set all possible terminators in matrix with columns: fold, energy, start, length, end
	my @resultsMTX;

	for (my $hlxline=0; $hlxline<scalar(@helices)-2; ++$hlxline)
	{
		my @helix = (split(/\s+/, $helices[$hlxline]));

		my $energy=0;
		my $start=0;
		
		if (scalar(@helix) == 3)
		{
		$energy = ($helix[1] =~ /(\-?\d+\.\d+)/)[0];
		$start = $helix[2];
		}	
		if (scalar(@helix) == 4)
		{
		$energy = ($helix[2] =~ /(\-?\d+\.\d+)/)[0];
		$start = $helix[3];
		}
		
		my $len = length($helix[0]);
		my $end = $start + $len - 1;
		
		# Choose P1 that:
		
		# starts before AT
		if ($start < $ATminenergy[2])
		{
			# ends after AT starts and up to the TM
			if (($end > $ATminenergy[2]+5) && ($end < $TMminenergy[2]))
			{
				# fold, energy, start, length, end
				my @P1 = ($helix[0], $energy, $start, $len, $end);
				push(@resultsMTX, \@P1);
			}
		}
	}

	if (@resultsMTX)
	{
		# Choose longest
		return @{(sort { $b->[3] <=> $a->[3] } @resultsMTX)[0]};
	}
	else
	{
		return;
	}
}

####################################################################
# avgProb                                                          #
####################################################################
# Calculates the probabilities from the RNAfold output.
#
# Irrelevant for current classifier.
#
####################################################################
sub avgProb
{
	# get pairs
	my $type = $_[0];
	my $sum = 0;
	my $count = 0;
	my $SLsum = 0;
	my $SLcount = 0;
	my $pair;
	my $Lpair;
	my $Rpair;
	$allmul=1;
	$SLmul=1;
	
	open(my $TM_RNAfile, "<", $tmpdir.$type."_".$name."_".$ID."_ss.ps") or die "can't open file $!";
	
	
	my $line = <$TM_RNAfile>;
	
	while(defined $line)
	{
		if ($line =~ /\/pairs/)
		{
			$line = <$TM_RNAfile>;
			
			while ((defined $line) && ($line =~ /\[(\d+)\s+(\d+)\]/))
			{
				$Lpair = $1;
				$Rpair = $2;

				my $pairfile = $tmpdir.$type."_".$name."_".$ID."_dp.ps";
				$pair = `grep "^$Lpair $Rpair [0-9][0-9]*\.[0-9][0-9]* ubox" $pairfile`;
				$sum += (split(/\s+/, $pair))[2];
				$allmul *= (split(/\s+/, $pair))[2];
				++$count;
				
				if (($Lpair >= $_[1]) && ($Lpair <= $_[2]))
				{
					$SLsum += (split(/\s+/, $pair))[2];
					$SLmul *= (split(/\s+/, $pair))[2];
					++$SLcount;
				}
				
				$line = <$TM_RNAfile>;
			}
		}
		
		$line = <$TM_RNAfile>;
	}
	
	close $TM_RNAfile;
	$allavg = $sum/$count;
	$SLavg = $SLsum/$SLcount;	
}
