#!/usr/bin/perl -w

use strict;

my $seq_dir = shift; ##### FULL PATH REQIRED!
$seq_dir =~ s/\/$//g;

my @inputs = @ARGV; ##### FULL PATH REQUIRED!

my $num=0;
my $q  =0;
for my $input (@inputs){
    my $name = 'na';
    if($input =~ /\/([^\/]+)\.fna/){
	$name = $1;
    }

    $q++ if(!($num % 30));

    my $savefile = $seq_dir . '/blast/sh.primers_blast/primers_blast'.$q.'.sh';
    open SAVE,">>$savefile";
    
    my $command = 'blastn -task blastn-short -strand plus -db ~/Desktop/output/scripts/vector_database/newconstseq.fa -outfmt 10 -evalue 1e-3 ';
    my $blastsave = $seq_dir . '/blast/out.primers_blast/'.$name.'.blast';
    $command .= '-query '.$input.' -out '.$blastsave;
    print SAVE "$command\n";
    
    $num++;

    close SAVE;
}
