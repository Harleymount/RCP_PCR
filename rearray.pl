#!/usr/bin/perl -w

use strict;
use Data::Dumper;

my $file = shift; # decent barcode list

# You need to edit here
##############################################################################
my @dest_db_plate_list = ('DB0001GOLD','DB0002GOLD');
my @dest_ad_plate_list = ('AD0001GOLD','AD0002GOLD');
##############################################################################



my $source;

my @ad_array;
my @db_array;

print '##################'."\n";
print '# Strain ID list #'."\n";
print '##################'."\n";

open FILE,$file;
while(<FILE>){
    chomp;
    my @array = split /\t/, $_;
    if($#array>4){
	my ($plate,$row,$column,$match,$uptag,$dntag) = @array;
	my $r = $row;
	my $c = $column;

	my $m_tag = ($match)? 'N2':'N1';
	
	$r =~ s/^[RC]0*//g;
	$c =~ s/^[RC]0*//g;
	
	my $name = "$plate$row$column$m_tag";

	$source->{$plate}->{$c}->{$r}->{name}  = $name;
	$source->{$plate}->{$c}->{$r}->{uptag} = $uptag;
	$source->{$plate}->{$c}->{$r}->{dntag} = $dntag;

	push @ad_array,$name if($name =~ /^AD/);
	push @db_array,$name if($name =~ /^DB/);

	print "$name\t$uptag\t$dntag\t$match\n";
    }
}
close FILE;

print "\n\n";



##### make source list

print '#############################'."\n";
print '# Source/library plate info #'."\n";
print '#############################'."\n";

for my $plate (sort{$a cmp $b} keys %$source){
    print "$plate\n";
    for my $c (1..24){
	for my $r (1..16){

	    my $row    = sprintf "%.1f", $r/10;
	    my $column = sprintf "%.1f", $c/10;

	    $row    =~ s/\.//g;
	    $column =~ s/\.//g;
	    
	    $row    = "R$row";
	    $column = "C$column";

	    my $name = "$plate$row$column".'N0';
	    if($source->{$plate}->{$c}->{$r}->{name}){
		$name = $source->{$plate}->{$c}->{$r}->{name};
	    }
	    print "$c\t$r\t$name\n";
	}
    }
    print "\n";
}

print "\n\n";


##### make destination list

print '##########################'."\n";
print '# Destination plate info #'."\n";
print '##########################'."\n";

my $ad_count = 0;
for my $d_plate (@dest_ad_plate_list){
    print "$d_plate\n";

    my $mem;
    for my $r (1..16){
	for my $c (1..24){
	    my $name = ($ad_array[$ad_count] || 'BLANK');
	    $mem->{$r}->{$c} = $name;
	    $ad_count++;
	}
    }

    for my $c (1..24){
	for my $r (1..16){
	    my $name = $mem->{$r}->{$c};
	    print "$c\t$r\t$name\n";
	}
    }
    print "\n";
}



my $db_count = 0;
for my $d_plate (@dest_db_plate_list){
    print "$d_plate\n";

    my $mem;
    for my $r (1..16){
	for my $c (1..24){
	    my $name = ($db_array[$db_count] || 'BLANK');
	    $mem->{$r}->{$c} = $name;
	    $db_count++;
	}
    }

    for my $c (1..24){
	for my $r (1..16){
	    my $name = $mem->{$r}->{$c};
	    print "$c\t$r\t$name\n";
	}
    }
    print "\n";
}
