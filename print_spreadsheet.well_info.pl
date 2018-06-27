#!/usr/bin/perl -w

use strict;
use Data::Dumper;

my $data = do shift;

#print Dumper $barcode_mem;

for my $type (sort{$a cmp $b}(keys %$data)){
    for my $P_TAG (sort{$a cmp $b}(keys %{$data->{$type}})){

	for my $R_TAG (1..16){
	    $R_TAG = sprintf "%.1f", $R_TAG/10;
	    $R_TAG =~ s/\.//g;
	    $R_TAG = "R$R_TAG";
	    
	    for my $C_TAG (1..24){
		$C_TAG = sprintf "%.1f", $C_TAG/10;
		$C_TAG =~ s/\.//g;
		$C_TAG = "C$C_TAG";
		
		my $value = $data->{$type}->{$P_TAG}->{$R_TAG}->{$C_TAG};
		my $array = $value->{barcode_called};
		
		if($array){
		    for my$element (@$array){

			#1  Plate name
			#2  Barcode type
			#3  Row
			#4  Column
			#5  Well stat
			#6  BC ID
			#7  #BC reads
			#8  #lox reads
			#9  UPTAG
			#10  UPTAG stat
			#11 DNTAG
			#12 DNTAG stat
			#13 #BC hit
			#14 Min_BC_hit_rate_imputed
			#15 U1 stat
			#16 U2 stat
			#17 D1 stat
			#18 D2 stat
			#19 #reads supported the stat
			#20 Rate of reads supported the stat
			#21 #lox hit
			#22 loxP stat
			#23 lox2272 stat
			#24 #reads supported the stat
			#25 Rate of reads supported the stat
			
			printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
			$P_TAG,$type,$R_TAG,$C_TAG,
			$value->{WELL_STAT},$element->{BC}->{ID},$value->{BC_all},$value->{lox_all},
			$element->{BC}->{UPTAG},$element->{BC}->{UPTAG_stat},
			$element->{BC}->{DNTAG},$element->{BC}->{DNTAG_stat},
			$element->{BC}->{hit},$element->{BC}->{min_BC_hit_rate_impt},
			$element->{BC}->{likely_stat}->{U1},$element->{BC}->{likely_stat}->{U2},$element->{BC}->{likely_stat}->{D1},$element->{BC}->{likely_stat}->{D2},
			$element->{BC}->{likely_stat_count},$element->{BC}->{stat_rate},
			$element->{lox}->{hit},
			$element->{lox}->{likely_stat}->{loxP},$element->{lox}->{likely_stat}->{lox2272},
			$element->{lox}->{likely_stat_count},$element->{lox}->{stat_rate};
		    }
		}else{
		    #1  Plate name
		    #2  Barcode type
		    #3  Row
		    #4  Column
		    #5  Well stat
		    
		    printf "%s\t%s\t%s\t%s\t%s\n",
		    $P_TAG,$type,$R_TAG,$C_TAG,
		    $value->{WELL_STAT};
		}
		
		
	    }
	}
    }
} 
