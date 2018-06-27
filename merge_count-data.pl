#!/usr/bin/perl -w

use strict;
use Data::Dumper;

my @files = @ARGV;

my $big_data;
for my $file (@files){
    my $data = do $file;
    
    while(my($key1,$value1) = each %$data){
	while(my($key2,$value2) = each %$value1){
	    while(my($key3,$value3) = each %$value2){
		while(my($key4,$value4) = each %$value3){
		    while(my($key5,$value5) = each %$value4){
			while(my($key6,$value6) = each %$value5){

			    $big_data->{$key1}->{$key2}->{$key3}->{$key4}->{$key5}->{$key6} ||=0;
			    $big_data->{$key1}->{$key2}->{$key3}->{$key4}->{$key5}->{$key6}+= $value6;

			    #while(my($key7,$value7) = each %$value6){
				#while(my($key8,$value8) = each %$value7){
				    
				    #$big_data->{$key1}->{$key2}->{$key3}->{$key4}->{$key5}->{$key6}->{$key7}->{$key8} ||=0;
				    #$big_data->{$key1}->{$key2}->{$key3}->{$key4}->{$key5}->{$key6}->{$key7}->{$key8}+= $value8;
				#}
			   #}
			}
		    }
		}
	    }
	}
    }
}


print Dumper $big_data;
