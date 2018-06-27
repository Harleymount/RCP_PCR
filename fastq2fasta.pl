#!/usr/bin/perl -w


use strict;
use Data::Dumper;

my $seq_dir = shift;
my @files     = @ARGV;

$seq_dir =~ s/\/$//g;
my $save_dir = "$seq_dir".'/fragmented_fasta';

for my $file_1 (@files){

    my $file_name = 'N_A';
    if($file_1 =~ /\/([^\/]+)\.fastq/){
	$file_name = $1;
    }    

    open FILE,$file_1;

    my @array;

    my $line  = 0;
    my $count = 0;
    my $split = 0;

    while(<FILE>){
	chomp;
	
	my $i = $line % 4;
	$array[$i] = $_;

	if($i % 4 == 3){

	    my $tag    = $array[0];
	    my $seq    = $array[1];
	    my $qscore = $array[3];
	    
	    my @temp = split /\s/, $tag;
	    $tag = shift @temp;
	    $tag =~ s/^\@//g;
	    $tag =~ tr/\:/\_/;

	    if(!($count%10000)){
		$split++;
		my $path = $save_dir . '/' . $file_name . '_' . $split . '.fna';

		close SAVE;
		open SAVE, ">$path";
		print "Writing\.\. $path\n";
        }
	    $qscore =~ s/[^\#]//g;
	    ##### if(length($qscore)<10){  ###### Threshold for num of "#" in quality score
	    print SAVE "\>$tag\n$seq\n";
        ##### }
        
	    $count++;
	}
	$line++;
    }
    close FILE;
}

