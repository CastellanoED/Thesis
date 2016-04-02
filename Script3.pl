#!/usr/bin/perl


### LIBS
use Bio::SeqIO;
use Getopt::Long;
use Math::BigInt;
use POSIX;

$|=1; #elimina el buffer

GetOptions('file=s' => \$file, #Files
'chr=s' => \$chr,
'fold=s' => \$xfold, #Clase funcional 
#'chromatin_state=s' => \$chromatin, 
'fixnum=i' => \$fixnum); #Número de sitios (n)
die("Usage = ???.pl \n") if(!$file);

#Put each sequence in a new array 
(!&readfasta($file,\@seqArray)) && do{die "$!\n";};

#Output file
open (RESULTS, ">>", "$chr:SFS:k.txt") or die "Error: Carpeta de resultados inexsistente"."\n";
#print RESULTS "chr	file	xfold	chromatin	m_dmel	m_dyak	folded_lines	polimorficos	foldedSFS	divergentes	monomorficos	k	kJC\n";

#FUNCTIONAL ANNOTATION DATA (CHROMATIN + CLASS_SITES)
for ($i=0; $i<=1; $i++) { # To take the first 3 sequences (Functional)
	push (@functional,$seqArray[$i]);
}

#POPULATION GENOMIC DATA (DGRP SEQ.)
for ($i=6; $i<=166; $i++) { # To take the rest of sequences (Population with Outgroups)
	push (@population,$seqArray[$i]);
}

unshift (@population,$seqArray[3],$seqArray[5]); # Ahora las dos primeras filas son dsim y dyak! xq hemos eliminado dere y dsec (Depende del formato del multifasta)

#print $population[0]."\n";
#print $population[1]."\n";

	#Sliding windows loop variables
	my $len = length($seqArray[0]);
	my $winsize = length($seqArray[0]);
	die ("Window size must be smaller than total seq. lenght(Total lenght:$len bp)\n") if ($len < $winsize);
	my $numArray = $#population;   # number of aligned sequences;
	
	##Combinatorial number
	my $n = $fixnum;
	#my $nbin2 = $n*($n-1)/2;
	
	#Number of discarded sites to consider
	my $discard = $numArray-$fixnum if $fixnum <= $numArray or die("discarded($fixnum) larger than total($n)");

	#GENE: for(my $i = 0; $i < $len; $i += $winsize)
	GENE: for(my $i = 0; $i < 1; $i += 1) # Lo importante es que entre solo una vez si queremos coger de las bases 8 a la 30
	{
		#my $positions= getPos($i, $functional[0], $xfold, $winsize);  
		my $positions= getPos(8, $functional[0], $xfold, 22);  

				
		my $stop = $i+$winsize;
		
		my $miss = 0;
		my $het  = 0;
		my $gaps = 0;
		my $analyzed = 0;
		my @clean_col = ();
		my @position_pol = ();
		my @position_mon = ();	
		my @sfs_folded_dmel	= ();	
		my @count = ();
		my @sfs_unfolded_dmel = ();	
		my @count2 = ();				
		my $folded70 ="";
		my $unfolded70="";

		POSITION: foreach my $pos (@$positions)
		{
					
			#get the column in the alignment for the current position
			my $col_out ="";		
			my $col_pop ="";		
				
				for (0 .. 1) 
				{
					my $nuc1=uc(substr($population[$_],$pos,1));
					$col_out .= $nuc1;
				}
				
				for (2 .. $numArray) 
				{
					my $nuc2=uc(substr($population[$_],$pos,1));
					$col_pop .= $nuc2;
				}
		
			#count missing and heterozigous sites in position
			#discard sites with many N's and gaps from the analysis
			my $count_miss = ($col_pop =~ tr/N//);
			my $count_het  = ($col_pop =~ tr/(M|R|W|S|Y|K)//);
			my $count_gaps = ($col_out =~ tr/-//);
			$miss += $count_miss;
			$gaps += $count_gaps;
			$het  += $count_het;
											
			my $total_non_valid = $count_miss + $count_het;
			next POSITION if $total_non_valid > $discard;                        
			$col_pop =~ s/(N|M|R|W|S|Y|K)//g ;   #remove N's, gaps or heterozigous sites if any
			next POSITION if length($col_pop) < $fixnum;
			my @col_split = split(//, $col_pop);  #transform scalar into array
			shuffle(\@col_split);
			my @slice = @col_split[0..$fixnum-1]; #get the first $fixnum elements of the shuffled array
			my $string = join("",@slice);
			push (@clean_col, $col_out.$string); 
			
			#filters passed, count column as valid and proceed to analyze	
			$analyzed++; 		
		}	
			
		foreach $_ (@clean_col) #colum w/o Ns or residual heterozygotes // Mira posicion por posicion en el alineamiento
		{
			$A = 0;
			$C = 0;
			$T = 0;
			$G = 0;
			$mon = 0;
			$pol = 0;
			$leng = 0;
			
			$leng = length $_;
			$dyak = substr($_, 0, 1); #Depende del formato del multifasta!!
#print $dyak."\t";
			$dsim = substr($_, 1, 1); #Depende del formato del multifasta!!
#print $dsim."\n";
			$sw_preconsensus = substr($_, 2, $leng); #Sec. polimorficas del DGRP
			$aligm = length $sw_preconsensus;
			
			# Contadores:		
			my $countA = ($sw_preconsensus =~ tr/A//);
			my $countC = ($sw_preconsensus =~ tr/C//);
			my $countT = ($sw_preconsensus =~ tr/T//);
			my $countG = ($sw_preconsensus =~ tr/G//);
			$A += $countA;
			$C += $countC;
			$T += $countT;
			$G += $countG;
			
			# Monomórficos:
			if ($A == $aligm) {
				
					$mon++;
					$FREQ = $aligm;
					$minor = "N";
					$freq = 0;
					$dmel = "A";
				}
				
			elsif ($C == $aligm) {
					
					$mon++;
					$FREQ = $aligm;
					$minor = "N";
					$freq = 0;
					$dmel = "C";
				}
				
			elsif ($T == $aligm) {
					
					$mon++;
					$FREQ = $aligm;
					$minor = "N";
					$freq = 0;
					$dmel = "T";
				}
				
			elsif ($G == $aligm) {
					
					$mon++;
					$FREQ = $aligm;
					$minor = "N";
					$freq = 0;
					$dmel = "G";
				}
							
			# Polimórficos mayoritarios i minoritarios: ------------> Si hay tres alelos segregando el que esta a más baja freqüencia no lo cuenta!!!
			if ($A < $aligm and $A != 0 and $A >= ($C+$T+$G)) 
			{ 
				$pol++;
				$dmel = "A"; 
				$FREQ = $A;
				
				if ($C > ($T+$G)) {	$minor = "C"; $freq = $C; }
				 
				if ($T > ($C+$G)) {	$minor = "T"; $freq = $T; } 
				
				if ($G > ($C+$T)) {	$minor = "G"; $freq = $G; } 
			}
			
			if ($C < $aligm and $C != 0 and $C >= ($A+$T+$G)) 
			{ 
				$pol++;			
				$dmel = "C"; 
				$FREQ = $C;
				
				if ($A > ($T+$G)) {	$minor = "A"; $freq = $A; }
				 
				if ($T > ($A+$G)) {	$minor = "T"; $freq = $T; } 
				
				if ($G > ($A+$T)) {	$minor = "G"; $freq = $G; } 
			}
			
			if ($T < $aligm and $T != 0 and $T >= ($C+$A+$G)) 
			{ 
				$pol++;			
				$dmel = "T"; 
				$FREQ = $T;
				
				if ($C > ($A+$G)) {	$minor = "C"; $freq = $C; }
				 
				if ($A > ($C+$G)) {	$minor = "A"; $freq = $A; } 
				
				if ($G > ($C+$A)) {	$minor = "G"; $freq = $G; } 
			}
			
			if ($G < $aligm and $G != 0 and $G >= ($C+$A+$T)) 
			{ 
				$pol++;
				$dmel = "G"; 
				$FREQ = $G;
				
				if ($C > ($A+$T)) {	$minor = "C"; $freq = $C; }
				 
				if ($A > ($C+$T)) {	$minor = "A"; $freq = $A; } 
				
				if ($T > ($C+$A)) {	$minor = "T"; $freq = $T; } 
			}
		
			my $positions = $dmel."\t".$FREQ."\t".$minor."\t".$freq."\t".$dsim."\t".$dyak."\t".$pol."\t".$mon."\t".$aligm;	
			push (@position_pol, $positions) if ($pol >= 1); #En los casos donde una posicion haya dos alelos segregando al 50% este contador es igual a 2
			push (@position_mon, $positions) if ($mon == 1 and $pol == 0);		
	
		}	
				
		# Número de linias folded y unfolded:
		#my $unfolded_lines = $aligm-1; # menos 1 porq el número de lineas totales significaria que es monomorfico!
		my $folded_lines = ceil($aligm/2); #redondea hacia arriba
		########
		
		foreach (@position_pol) 
		{
			my @col1 = split (/\t/,$_);
			$col = join("\t",@col1);
		}
			
		foreach (@position_mon) 
		{
			my @col2 = split (/\t/,$_);
			$col2 = join("\t",@col2);
		}		
		
		# Total monomórficos:
		$total_mon = $#position_mon+1;
		# Total polimórficos:
		$total_pol = $#position_pol+1;
		# Total nt:
		$m = $total_pol+$total_mon;
		########
		
		#print "$chr-$xfold-$chromatin-$i-$stop-$m\n";
		
		# Polimórfico NO derivado (sin tener en cuenta la divergencia):
		foreach my $nt (@position_pol) 
		{
			my @colm = split (/\t/,$nt);
			push (@sfs_folded_dmel, ":$colm[3]:"); 
			@sfs_folded_dmel = sort {$a cmp $b} @sfs_folded_dmel;
			$sfs_folded_dmel = join("",@sfs_folded_dmel);
		}	
		
		# Calculo del site frequency spectrum folded:	
		for ($r = 1; $r <= "$folded_lines"; $r++)	
		{
			my $count = grep /:$r:/, @sfs_folded_dmel;
			push (@count, $count);
			$folded70 = join(":",@count);
		}
		#print "$folded70\n";
		
		#Gaps polimorfismo:
		$gaps_pol = 0; #Si los gaps se solapan entre outgroups solo cuenta una vez
		$gaps_pol_dsim = 0; 
		$gaps_pol_dyak = 0;
		########
		
		#Tipo de polimorfismo:
		$pol_1 = 0;	#Polimorfismo de alta calidad
		$pol_2 = 0; #Polimorfismo de alta calidad
		$pol_3 = 0;
		$pol_4 = 0;
		$pol_5 = 0;
		$pol_6 = 0;
		$pol_7 = 0;
		$pol_8 = 0;
		$pol_der = 0; #Polimorficos derivados totales	
		########
			
		# Polimórfico derivado (hay tendencia a contar más minor que major y creer que dsim es más paracido a dmel que dmel a dyak):
		foreach my $nt (@position_pol) 
		{
			my @colu = split (/\t/,$nt);
			my $major = $colu[0];
			my $minor = $colu[2];	
			my $major_freq = $colu[1];
			my $minor_freq = $colu[3];	
			my $dsim = $colu[4];
			my $dyak = $colu[5];

#print "$nt\n";				
#print "$major	$minor	$major_freq	$minor_freq	$dyak	$dsim\n";
			
			if ($dyak eq "-") { $gaps_pol_dyak++; }
			if ($dsim eq "-") { $gaps_pol_dsim++; }
			if ($dsim eq "-" or $dyak eq "-") 
			{ 
				$gaps_pol++; 
				
			} else {
				
				$pol_der++;
					
				if ($dsim eq $dyak) #lo más parsimonioso
				{
					if ($major eq $dyak) {$der = $minor_freq; $pol_1++;} #el menor es el derivado, pues el mayor es como dsim y dyak
					
					if ($minor eq $dyak) {$der = $major_freq; $pol_2++;} #el mayor es el derivado, pues el menor es como dsim y dyak 
					
					if ($minor ne $dyak and $major ne $dyak) {$der = $minor_freq; $pol_3++;} #los dos son derivados (sólo consideramos el menor derivado), pues ninguno es como dsim o dyak 
				}
	
				if ($dsim ne $dyak) 
				{ 
					if ($major eq $dsim and $minor ne $dyak) {$der = $minor_freq; $pol_4++;}  #el mayor es como simulans y el menor es nuevo
	
					if ($minor eq $dsim and $major ne $dyak) {$der = $major_freq; $pol_5++;}  #el menor es como simulans y el mayor es nuevo
					
					if ($major eq $dsim and $minor eq $dyak) {$der = $minor_freq; $pol_6++;}  #el mayor es como simulans y el menor como yakuba (pol ancestral)
	
					if ($minor eq $dsim and $major eq $dyak) {$der = $major_freq; $pol_7++;}  #el menor es como simulans y el mayor como yakuba (pol ancestral)
	
					if ($minor ne $dsim and $minor ne $dyak and $major ne $dsim and $major ne $dyak) {$der = $minor_freq; $pol_8++;}  #todos son diferentes (sólo consideramos el menor derivado)
				}
			
			} 
			
			$derived = "$der";
			push (@sfs_unfolded_dmel, ":$derived:");
			@sfs_unfolded_dmel = sort {$a cmp $b} @sfs_unfolded_dmel;
			$sfs_unfolded_dmel = join("",@sfs_unfolded_dmel);
		
		}		
		
		## Calculo del site frequency spectrum unfolded:	
		#for ($t = 1; $t <= "$unfolded_lines"; $t++)	
		#{
			#my $count2 = grep /:$t:/, @sfs_unfolded_dmel;
			#push (@count2, $count2);
			#$unfolded70 = join(":",@count2);
		#}
		#print "$unfolded70\n";
		########		
		
		#Gaps monomorfismo:
		$gaps_mon = 0; #Si los gaps se solapan entre outgroups solo cuenta una vez
		$gaps_mon_dsim = 0; 
		$gaps_mon_dyak = 0;
		########
		
		#Tipo de monomorficos:
		$i0 = 0;	
		$i1 = 0;
		$i2 = 0;
		$i3 = 0;
		$i4 = 0;
		$i5 = 0;
		$i6 = 0;
		$i7 = 0;
		$i8 = 0;
		$i9 = 0;
		########
		
		# Monomórfico derivado (asumimos q dsim es más paracido a dmel que dmel a dyak):
		foreach (@position_mon) 
		{

#print $_."\n";

			my @row = split (/\t/,$_);
			my $dmel_mon = $row[0];
			my $dsim_mon = $row[4];
			my $dyak_mon = $row[5];
							
			if ($dyak_mon eq "-") { $gaps_mon_dyak++; }
			if ($dsim_mon eq "-") { $gaps_mon_dsim++; }
			
			# Divergencia polarizada necesitamos dos outgroups:
			if ($dsim_mon eq "-" or $dyak_mon eq "-") 
			{ 
				$gaps_mon++; 
			} else {
					
					$i0++; #no gaps
					if ($dmel_mon eq $dsim_mon and $dmel_mon eq $dyak_mon) {$i1++;} #monomorfico polarizado (invariable)
					if ($dmel_mon eq $dsim_mon and $dmel_mon ne $dyak_mon) {$i2++;} #div polarizada especifica dyak (o especifica del linaje dsim-dmel)
					if ($dmel_mon ne $dsim_mon and $dmel_mon eq $dyak_mon) {$i3++;} #div polarizada especifica dsim 
					if ($dmel_mon ne $dsim_mon and $dmel_mon ne $dyak_mon and $dsim_mon eq $dyak_mon) {$i4++;} #div polarizada dmel
					if ($dmel_mon ne $dsim_mon and $dmel_mon ne $dyak_mon and $dsim_mon ne $dyak_mon) {$i5++;} #div en todos los outgroups (incluido dmel)
			}
				 
			#Divergencia no polarizada spp1(dmel) Vs spp2 (dsim/dyak):
			if ($dmel_mon eq $dsim_mon) {$i6++;} 						#monomorfico con dsim
			if ($dmel_mon ne $dsim_mon and $dsim_mon ne "-") {$i7++;} 	#divergente con dsim
			if ($dmel_mon eq $dyak_mon) {$i8++;} 						#monomorfico con dyak
			if ($dmel_mon ne $dyak_mon and $dyak_mon ne "-") {$i9++;} 	#divergente con dyak			 
			
		}	
				
		$m_dyak = ($m-($gaps_mon_dyak+$gaps_pol_dyak));
		$m_dsim = ($m-($gaps_mon_dsim+$gaps_pol_dsim));
	
		###############################################################################################################################################################################################
		
		#Corrección JC	
		
			# Divergencia polarizada específica de dyak (o del linaje dsim-dmel)
			#my $m_dyak = $#position_mon+1-$gaps_mon_dyak;
			my $ks_dyak = 0;
			if ($m_dyak == 0) {
				$ks_dyak = -0.001;
			} else {
				$ks_dyak = $i2/$m_dyak; 
				$ks_dyak = sprintf("%.6f",$ks_dyak);
			}
			my $ks_dyak_JC=JC_correction($m_dyak,$ks_dyak);
			$ks_dyak_JC = sprintf("%.6f",$ks_dyak_JC);
			
			# Divergencia polarizada específica de dsim 
			#my $m_dsim = $#position_mon+1-$gaps_mon_dsim;
			my $ks_dsim = 0;
			if ($m_dsim == 0) {
				$ks_dsim = -0.001;
			} else {
				$ks_dsim = $i3/$m_dsim;
				$ks_dsim = sprintf("%.6f",$ks_dsim);
			}
			my $ks_dsim_JC=JC_correction($m_dsim,$ks_dsim);
			$ks_dsim_JC = sprintf("%.6f",$ks_dsim_JC);
			
			# Divergencia polarizada específica de dmel 
			my $m_dmel = $#position_mon+1-$gaps_mon_dyak-$gaps_mon_dsim;
			my $ks_dmel = 0;
			if ($m_dmel == 0) {
				$ks_dmel = -0.001;
			} else {
				$ks_dmel = $i4/$m_dmel;
				$ks_dmel = sprintf("%.6f",$ks_dmel);
			}
			my $ks_dmel_JC=JC_correction($m_dmel,$ks_dmel);
			$ks_dmel_JC = sprintf("%.6f",$ks_dmel_JC);
		
			# Divergencia no polarizada entre dyak-dsim
			#my $m_dyak = $#position_mon+1-$gaps_mon_dyak;
			my $ks_dyakdmel = 0;
			if ($m_dyak == 0) {
				$ks_dyakdmel = -0.001;
			} else {
				$ks_dyakdmel = $i9/$m_dyak; 
				$ks_dyakdmel = sprintf("%.6f",$ks_dyakdmel);
			}
			my $ks_dyakdmel_JC=JC_correction($m_dyak,$ks_dyakdmel);
			$ks_dyakdmel_JC = sprintf("%.6f",$ks_dyakdmel_JC);
			
			
		###############################################################################################################################################################################################
		
		@filenameparts = split ("/",$file);
		@filename = split ("_",$filenameparts[3]);

		print RESULTS "$chr	$filename[2]	$filename[1]	$filename[3]	$filename[4]	$m	$m_dyak	$total_pol	$folded70	$i9	$ks_dyakdmel	$ks_dyakdmel_JC\n" if $m > 0;
		
		###############################################################################################################################################################################################

}
	#print "Finished SFS and divergence for $filename with $xfold and $fixnum\n";



######## ____________FUNCTIONS____________ ########
###################################################

sub getPos {
	my $i           = shift;
	my $recoded_seq = shift; 	
	my $class_site 	= shift;	
	my $winsize     = shift;

	my $slide_recoded = substr( $recoded_seq, $i, $winsize );

	my @positions3 = ();
		
		while ( $slide_recoded =~ m/($class_site)/gi) 
		{
			my $pos = length($`);
			$pos = $pos + $i if $i > 0;			             
			push( @positions3, $pos );
		}
			

	return \@positions3; 
}
sub readfasta($$){
	my ($file,$seqs)=@_;
	(!open(FASTA,$file)) && do{return 0;};	
	my $id="";
	my $seq="";
	while(<FASTA>){
		chomp;	
		if ($_=~/^>(\S+)/){
			push(@$seqs,$seq) if ($id ne "");
			$seq="";
			$id=$1;
		}else{
			$seq.=$_;
		}
#print $id."\n";
	}
#print $id."\n";
	push(@$seqs,$seq) if ($id ne "");
	close(FASTA);
	return 1;
}
sub JC_correction {
	my ($numSites,$divergence)=@_;
	if ($numSites != '0') {
		#my $divergence2correct = $divergence/$numSites;

		#correccion JC
		if (((1-((4/3)*$divergence2correct)) > 0) ){
			#$divergence_JC = (-(3/4)*log(1-((4/3)*$divergence2correct))); 
			$divergence_JC = (-(3/4)*log(1-((4/3)*$divergence))) or $divergence_JC = 0 ; 
			#$divergence_JC = ($divergence_JC*$numSites);
		} else {
			$divergence_JC = ''; 
		}
	} else {
		$divergence_JC='';
	}
	return ($divergence_JC);
}
sub shuffle {
	#fisher yates shuffle
	srand(8);	
    my $array = shift;
    my $i = @$array;
    while ( --$i )
    {
        my $j = int rand( $i+1 );
        @$array[$i,$j] = @$array[$j,$i];
    }
}
