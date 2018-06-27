#!/usr/bin/perl -w

use strict;
use Data::Dumper;

my $data = do shift;

my $barcode_mem;

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
		
		my $array = $data->{$type}->{$P_TAG}->{$R_TAG}->{$C_TAG}->{barcode_called};
		
		for my $element (@$array){
		    my $UPTAG = $element->{BC}->{UPTAG};
		    my $DNTAG = $element->{BC}->{DNTAG};
		    
		    $barcode_mem->{$type}->{UPTAG}->{$UPTAG}->{"$P_TAG\:$R_TAG\:$C_TAG"} = 1;
		    $barcode_mem->{$type}->{DNTAG}->{$DNTAG}->{"$P_TAG\:$R_TAG\:$C_TAG"} = 1;
		}
	    }
	}
    }
}

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
		
		#####################################
		##### Define the status of well #####
		#####################################
		my $WELL_STAT = '-';
		if($array){
		    
		    $WELL_STAT = 'S';
		    $WELL_STAT = 'M' if($#$array>0);
		    
		    my $id=0;


		    my $error_sig = 0;
		    my $error_mem;

		    for my $i (0..$#$array){
			my $element = $array->[$i];
			$id++;

			$data->{$type}->{$P_TAG}->{$R_TAG}->{$C_TAG}->{barcode_called}->[$i]->{BC}->{ID} = "$id";

			my $UPTAG = $element->{BC}->{UPTAG};
			my $DNTAG = $element->{BC}->{DNTAG};

			my $UPTAG_stat = '+';
			my $DNTAG_stat = '+';					
			
			if(scalar(values %{$barcode_mem->{$type}->{UPTAG}->{$UPTAG}}) > 1){
			    $error_mem->{NS} = 1;
			    $error_sig = 1;
			    $UPTAG_stat = 'NS';
			}
			if(scalar(values %{$barcode_mem->{$type}->{DNTAG}->{$DNTAG}}) > 1){
			    $error_mem->{NS} = 1;
			    $error_sig = 1;
			    $DNTAG_stat = 'NS';
			}
			
			$data->{$type}->{$P_TAG}->{$R_TAG}->{$C_TAG}->{barcode_called}->[$i]->{BC}->{UPTAG_stat} = $UPTAG_stat;
			$data->{$type}->{$P_TAG}->{$R_TAG}->{$C_TAG}->{barcode_called}->[$i]->{BC}->{DNTAG_stat} = $DNTAG_stat;

			for my $read_type ('BC','lox'){
			    while(my($key,$value) = each %{$element->{$read_type}->{likely_stat}}){

				if($value eq '-'){
				    $error_mem->{E} = 1;
				    $error_sig = 1;
				}
				if($value eq 'Ab'){
				    $error_mem->{Ab} = 1;
				    $error_sig = 1;
				}

				
			    }
			}
		    }
		    if($error_sig){
			$WELL_STAT = join ';', ($WELL_STAT,keys %$error_mem);
		    }
		}
		$data->{$type}->{$P_TAG}->{$R_TAG}->{$C_TAG}->{WELL_STAT} = $WELL_STAT;
		##### END
	    
	    }
	}
    }
} 

print Dumper $data;
