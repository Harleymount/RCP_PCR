#!/usr/bin/perl -w

use strict;
use Data::Dumper;

my $data = do shift;


my $count;

while(my($category,$value1) = each %$data){
    while(my($P_TAG,$value2) = each %$value1){
	while(my($R_TAG,$value3) = each %$value2){
	    while(my($C_TAG,$value4) = each %$value3){
		while(my($read,$value) = each %$value4){

		    ########################
		    # Count DB-BC barcodes #
		    ########################
		    if($category eq 'DB-BC'){

			my $UPTAG   = $value->{R1}->{UPTAG}  || 'N_A';
			my $cDNTAG  = $value->{R2}->{cDNTAG} || 'N_A';
			my $DNTAG   = 'N_A';

			if($cDNTAG =~ /^[ATGCN]+$/){
			    $DNTAG = reverse $cDNTAG;
			    $DNTAG =~ tr/atgcATGC/tacgTACG/;
			}

			my $tags = join ':', ($UPTAG,$DNTAG);

			my $stat = join ':', 
			($value->{R1}->{'DBU1-primer'},$value->{R1}->{'cDBU2-primer'},$value->{R2}->{'cDBD1-primer'},$value->{R2}->{'DBD2-primer'},$value->{R1}->{'lox2272'});


			 $count->{$P_TAG}->{$R_TAG}->{$C_TAG}->{$category}->{$tags}->{$stat} ||= 0;
			 $count->{$P_TAG}->{$R_TAG}->{$C_TAG}->{$category}->{$tags}->{$stat}++;
		    }


		    ########################
		    # Count AD-BC barcodes #
		    ########################
		    if($category eq 'AD-BC'){

			my $UPTAG   = $value->{R1}->{UPTAG}  || 'N_A';
			my $cDNTAG  = $value->{R2}->{cDNTAG} || 'N_A';
			my $DNTAG   = 'N_A';

			if($cDNTAG =~ /^[ATGCN]+$/){
			    $DNTAG = reverse $cDNTAG;
			    $DNTAG =~ tr/atgcATGC/tacgTACG/;
			}

			my $tags = join ':', ($UPTAG,$DNTAG);

			my $stat = join ':', 
			($value->{R1}->{'ADU1-primer'},$value->{R1}->{'cADU2-primer'},$value->{R2}->{'cADD1-primer'},$value->{R2}->{'ADD2-primer'},$value->{R2}->{'loxP'});


			 $count->{$P_TAG}->{$R_TAG}->{$C_TAG}->{$category}->{$tags}->{$stat} ||= 0;
			 $count->{$P_TAG}->{$R_TAG}->{$C_TAG}->{$category}->{$tags}->{$stat}++;
		    }
		    
		    ######################
		    # Count DB-lox reads #
		    ######################
		    if($category eq 'DB-lox'){

			my $merged_UPTAG = $value->{'merged-UPTAG'};

			my $stat = join ':', ($value->{R1}->{'cloxP'},$value->{R2}->{'clox2272'},$value->{R1}->{'DBU1-primer'},$value->{R2}->{'DBU2-primer'});

			$count->{$P_TAG}->{$R_TAG}->{$C_TAG}->{$category}->{$merged_UPTAG}->{$stat} ||= 0;
			$count->{$P_TAG}->{$R_TAG}->{$C_TAG}->{$category}->{$merged_UPTAG}->{$stat}++;
		    }

		    ######################
		    # Count AD-lox reads #
		    ######################
		    if($category eq 'AD-lox'){

			my $merged_DNTAG = $value->{'merged-DNTAG'};

			my $stat = join ':', ($value->{R1}->{'cloxP'},$value->{R2}->{'clox2272'},$value->{R1}->{'ADD1-primer'},$value->{R2}->{'ADD2-primer'});

			$count->{$P_TAG}->{$R_TAG}->{$C_TAG}->{$category}->{$merged_DNTAG}->{$stat} ||= 0;
			$count->{$P_TAG}->{$R_TAG}->{$C_TAG}->{$category}->{$merged_DNTAG}->{$stat}++;
		    }




		}
	    }
	}
    }
}

print Dumper $count;
