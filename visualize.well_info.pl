#!/usr/bin/perl -w

use strict;
use Data::Dumper;

my $data = do shift;

#print Dumper $barcode_mem;

for my $type (sort{$a cmp $b}(keys %$data)){
    for my $P_TAG (sort{$a cmp $b}(keys %{$data->{$type}})){

	print "$P_TAG\n";

	my $sum = 0;
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
		
		my $vis = 0;
		$vis = 1 if($value->{WELL_STAT} eq 'S');
		
		print "$vis\t";
		$sum++ if($vis == 1)
	    }
	    print "\n";
	}
	print "\n\# of good BCs\: $sum\n";
	print "\n";
    }
} 
