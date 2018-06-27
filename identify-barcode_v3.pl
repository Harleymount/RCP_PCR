#!/usr/bin/perl -w

use strict;
use Data::Dumper;

my $primers_fasta = '/Users/Harley/Desktop/my_data/data/newconstseq.fna';
my $bar2num_file  = '/Users/Harley/Desktop/my_data/bar2num2.txt';
my $tag_assignment = '/Users/Harley/Desktop/my_data/tag_assignment6282017.csv';


my $seq_dir = shift; ##### FULL PATH REQIRED!
$seq_dir =~ s/\/$//g;

my $R1_id = shift;
my $R2_id = shift;

my $seq_DB  = &get_seq($primers_fasta);
my $bar2num = &get_bar2num($bar2num_file);
my $tag_strct = &get_tag($tag_assignment);

my $primers_R1_blast = $seq_dir.'/blast/out.primers_blast/'. $R1_id.'.blast';
my $primers_R2_blast = $seq_dir.'/blast/out.primers_blast/'. $R2_id.'.blast';
my $fasta_R1         = $seq_dir.'/fragmented_fasta/'.  $R1_id.'.fna';
my $fasta_R2         = $seq_dir.'/fragmented_fasta/'.  $R2_id.'.fna';

my $seq_reads;
$seq_reads->{R1} = &get_seq($fasta_R1);
$seq_reads->{R2} = &get_seq($fasta_R2);

#print Dumper $seq_reads;


my @files = ($primers_R1_blast,$primers_R2_blast);

my $data;
for my $file (@files){

    my $read_direction = 'R1';
    if($file =~ /\_R2\_/){
	$read_direction = 'R2';
    }

    open FILE,$file;
    while(<FILE>){
	chomp;
	my @array = split /\,/,$_;
	if($#array>5){
	    my ($read,$template,$str_on_read,$end_on_read,$str,$end,$evalue) =
		@array[0,1,6,7,8,9,10];

	    if($str_on_read<$end_on_read && $str == 1 && $end == $seq_DB->{$template}->{length}){ ########## THINK ABOUT THIS
		my $element;
		$element->{str}    = $str_on_read;
		$element->{end}    = $end_on_read;
		$element->{evalue} = $evalue;
		$data->{$read}->{$read_direction}->{$template} = $element;
	    }
	}
    }
    close FILE;
}

#print Dumper $data;



my $data2;
while(my($read,$value) = each %$data){
    my $R1     = $value->{R1};
    my $R2     = $value->{R2};
    my $seq;
    $seq->{R1} = $seq_reads->{R1}->{$read}->{seq};
    $seq->{R2} = $seq_reads->{R2}->{$read}->{seq};



    # PS1.0 and PS2.0 at the reasonable positions? (end at <40bp)
    # Row and Column priming sites from the reasonable positions? (start at <50bp)
    # Which category? DB-BC / DB-lox / AD-BC / AD-lox?
    # Specific category assignment?
    # Row, Column and Plate tags?

    my ($category,$P1_seq,$P2_seq,$R_seq,$C_seq) = &assign_category($value,$R1,$R2,$seq);


    my $P1_num = ($P1_seq && length($P1_seq) < 12)? &barcode_matching($bar2num,$P1_seq):0;
    my $P2_num = ($P2_seq && length($P2_seq) < 12)? &barcode_matching($bar2num,$P2_seq):0;
    my $R_num  = ($R_seq  && length($R_seq)  < 12)? &barcode_matching($bar2num,$R_seq) :0;
    my $C_num  = ($C_seq  && length($C_seq)  < 12)? &barcode_matching($bar2num,$C_seq) :0;

    #############################
    # Process for each category #
    #############################

    if($category && $P1_num*$P2_num*$R_num*$C_num){
	
	$P1_num = sprintf "%.1f", $P1_num/10;
	$P1_num =~ s/\.//g;
	$P2_num = sprintf "%.1f", $P2_num/10;
	$P2_num =~ s/\.//g;
	$R_num  = sprintf "%.1f", $R_num/10;
	$R_num  =~ s/\.//g;
	$C_num  = sprintf "%.1f", $C_num/10;
	$C_num  =~ s/\.//g;

       	my ($P_TAG,$R_TAG,$C_TAG) = ("P$P1_num\-P$P2_num","R$R_num","C$C_num");
	#print "$category\t$P_TAG\t$R_TAG\t$C_TAG\n";


	my $stat;	
	if($category eq 'DB-BC'){

	    ##############################
	    # Process for DB-BC category #
	    ##############################
	    $stat = &analyze_DB_BC_stat($value,$R1,$R2,$seq);
      
	}elsif($category eq 'AD-BC'){

	    ##############################
	    # Process for AD-BC category #
	    ##############################    
	    $stat = &analyze_AD_BC_stat($value,$R1,$R2,$seq);
	    
	}elsif($category eq 'DB-lox'){

	    ###############################
	    # Process for DB-lox category #
	    ###############################
	    $stat = &analyze_DB_lox_stat($value,$R1,$R2,$seq);
	    
	}elsif($category eq 'AD-lox'){

	    ###############################
	    # Process for AD-lox category #
	    ###############################
	    $stat = &analyze_AD_lox_stat($value,$R1,$R2,$seq);
	    
	}
	my ($type1,$type2) = split /\-/, $category;
	if($tag_strct->{$type1}->{$type2}->{$P_TAG}){
	    my $plate_name = $tag_strct->{$type1}->{$type2}->{$P_TAG};
	    $data2->{$category}->{$plate_name}->{$R_TAG}->{$C_TAG}->{$read} = $stat;
	}
    }
}

print Dumper $data2;

















############################################################
# functions
############################################################

sub analyze_DB_BC_stat(){

    my $value = shift;
    my $R1    = shift;
    my $R2    = shift;
    my $seq   = shift;

    my $stat;
    $stat->{R1}->{'DBU1-primer'}  = 'there';
    $stat->{R1}->{'cDBU2-primer'} = 'absent';
    $stat->{R1}->{'UPTAG'}        = 'absent';
    $stat->{R1}->{'lox2272'}      = 'absent';
    $stat->{R2}->{'DBD2-primer'}  = 'there';
    $stat->{R2}->{'cDBD1-primer'} = 'absent';
    $stat->{R2}->{'cDNTAG'}       = 'absent';
    
    
    ##### Existence check; UPTAG
    if($R1->{'cDBU2-primer'}){
	if($R1->{'DBU1-primer'}->{end} < $R1->{'cDBU2-primer'}->{str}){
	    my $TAG_str = $R1->{'DBU1-primer'}->{end} +1 -1;
	    my $TAG_end = $R1->{'cDBU2-primer'}->{str}-1 -1;
	    
	    $stat->{R1}->{'cDBU2-primer'} = 'there';
	    $stat->{R1}->{'UPTAG'}        = substr $seq->{R1},$TAG_str,$TAG_end-$TAG_str+1;		    
	} 
    }
    
    ##### Existence check; lox2272
    if($R1->{'cDBU2-primer'} && $R1->{'lox2272'}){
	if($R1->{'cDBU2-primer'}->{end} < $R1->{'lox2272'}->{str}){
	    my $offset = $R1->{'lox2272'}->{str} - $R1->{'cDBU2-primer'}->{end} -1;
	    $stat->{R1}->{'lox2272'} = 'there;offset='.$offset;
	} 
    }
    
    ##### Existence check; DNTAG
    if($R2->{'cDBD1-primer'}){
	if($R2->{'DBD2-primer'}->{end} < $R2->{'cDBD1-primer'}->{str}){
	    my $TAG_str = $R2->{'DBD2-primer'}->{end} +1 -1;
	    my $TAG_end = $R2->{'cDBD1-primer'}->{str}-1 -1;
	    
	    $stat->{R2}->{'cDBD1-primer'} = 'there';
	    $stat->{R2}->{'cDNTAG'}       = substr $seq->{R2},$TAG_str,$TAG_end-$TAG_str+1;		    
	} 
    }
    
    ##### Sequence match check
    while(my($read,$value1) = each %$stat){
	while(my($primer,$value2) = each %$value1){
	    my @array = split /\;/, $value2;
	    if($array[0] eq 'there'){
		my $read_seq = substr $seq->{$read}, $value->{$read}->{$primer}->{str}-1, $value->{$read}->{$primer}->{end}-$value->{$read}->{$primer}->{str}+1;
		$array[0] = ($read_seq eq $seq_DB->{$primer}->{seq})? 'match':'no-match';
		
		$stat->{$read}->{$primer} = join ';', @array;
	    }
	}
    }
    
    return $stat;
}



sub analyze_AD_BC_stat(){

    my $value = shift;
    my $R1    = shift;
    my $R2    = shift;
    my $seq   = shift;

    my $stat;
    $stat->{R1}->{'ADU1-primer'}  = 'there';
    $stat->{R1}->{'cADU2-primer'} = 'absent';
    $stat->{R1}->{'UPTAG'}        = 'absent';
    
    $stat->{R2}->{'ADD2-primer'}  = 'there';
    $stat->{R2}->{'cADD1-primer'} = 'absent';
    $stat->{R2}->{'cDNTAG'}       = 'absent';
    $stat->{R2}->{'loxP'}         = 'absent';
    
    ##### Existence check; UPTAG
    if($R1->{'cADU2-primer'}){
	if($R1->{'ADU1-primer'}->{end} < $R1->{'cADU2-primer'}->{str}){
	    my $TAG_str = $R1->{'ADU1-primer'}->{end} +1 -1;
	    my $TAG_end = $R1->{'cADU2-primer'}->{str}-1 -1;
	    
	    $stat->{R1}->{'cADU2-primer'} = 'there';
	    $stat->{R1}->{'UPTAG'}        = substr $seq->{R1},$TAG_str,$TAG_end-$TAG_str+1;		    
	} 
    }
    
    ##### Existence check; DNTAG
    if($R2->{'cADD1-primer'}){
	if($R2->{'ADD2-primer'}->{end} < $R2->{'cADD1-primer'}->{str}){
	    my $TAG_str = $R2->{'ADD2-primer'}->{end} +1 -1;
	    my $TAG_end = $R2->{'cADD1-primer'}->{str}-1 -1;
	    
	    $stat->{R2}->{'cADD1-primer'} = 'there';
	    $stat->{R2}->{'cDNTAG'}       = substr $seq->{R2},$TAG_str,$TAG_end-$TAG_str+1;		    
	} 
    }
    
    ##### Existence check; loxP
    if($R2->{'cADD1-primer'} && $R2->{'loxP'}){
	if($R2->{'cADD1-primer'}->{end} < $R2->{'loxP'}->{str}){
	    my $offset = $R2->{'loxP'}->{str} - $R2->{'cADD1-primer'}->{end} -1;
	    $stat->{R2}->{'loxP'} = 'there;offset='.$offset;
	} 
    }
    
    ##### Sequence match check
    while(my($read,$value1) = each %$stat){
	while(my($primer,$value2) = each %$value1){
	    my @array = split /\;/, $value2;
	    if($array[0] eq 'there'){
		my $read_seq = substr $seq->{$read}, $value->{$read}->{$primer}->{str}-1, $value->{$read}->{$primer}->{end}-$value->{$read}->{$primer}->{str}+1;
		$array[0] = ($read_seq eq $seq_DB->{$primer}->{seq})? 'match':'no-match';
		
		$stat->{$read}->{$primer} = join ';', @array;
	    }
	}
    }
    
    return $stat;
}


sub analyze_DB_lox_stat(){

    my $value = shift;
    my $R1    = shift;
    my $R2    = shift;
    my $seq   = shift;

    my $stat;
    $stat->{R1}->{'DBloxP-primer'}    = 'there';
    $stat->{R1}->{'cloxP'}            = 'absent';
    $stat->{R1}->{'DBU1-primer'}      = 'absent';
    $stat->{R1}->{'UPTAG-frag'}       = 'absent';
    
    $stat->{R2}->{'DBlox2272-primer'} = 'there';
    $stat->{R2}->{'clox2272'}         = 'absent';
    $stat->{R2}->{'DBU2-primer'}      = 'absent';
    $stat->{R2}->{'cUPTAG-frag'}      = 'absent';
    
    ##### Existence check; cloxP
    if($R1->{'cloxP'}){
	if($R1->{'DBloxP-primer'}->{end} < $R1->{'cloxP'}->{str}){
	    my $offset = $R1->{'cloxP'}->{str} - $R1->{'DBloxP-primer'}->{end} -1;
	    $stat->{R1}->{'cloxP'} = 'there;offset='.$offset;
	} 
    }
    
    ##### Existence check; DBU1-primer
    if($R1->{'cloxP'} && $R1->{'DBU1-primer'}){
	if($R1->{'cloxP'}->{end} < $R1->{'DBU1-primer'}->{str}){
	    my $offset = $R1->{'DBU1-primer'}->{str} - $R1->{'cloxP'}->{end} -1;
	    $stat->{R1}->{'DBU1-primer'} = 'there;offset='.$offset;
	    $stat->{R1}->{'UPTAG-frag'}  = substr $seq->{R1},$R1->{'DBU1-primer'}->{end},length($seq->{R1}) - $R1->{'DBU1-primer'}->{end};
	} 
    }
    
    ##### Existence check; clox2272
    if($R2->{'clox2272'}){
	if($R2->{'DBlox2272-primer'}->{end} < $R2->{'clox2272'}->{str}){
	    my $offset = $R2->{'clox2272'}->{str} - $R2->{'DBlox2272-primer'}->{end} -1;
	    $stat->{R2}->{'clox2272'} = 'there;offset='.$offset;
	} 
    }
    
    ##### Existence check; DBU2-primer
    if($R2->{'clox2272'} && $R2->{'DBU2-primer'}){
	if($R2->{'clox2272'}->{end} < $R2->{'DBU2-primer'}->{str}){
	    my $offset = $R2->{'DBU2-primer'}->{str} - $R2->{'clox2272'}->{end} -1;
	    $stat->{R2}->{'DBU2-primer'} = 'there;offset='.$offset;
	    $stat->{R2}->{'cUPTAG-frag'}  = substr $seq->{R2},$R2->{'DBU2-primer'}->{end},length($seq->{R2}) - $R2->{'DBU2-primer'}->{end};
	} 
    }
    
    ##### Sequence match check
    while(my($read,$value1) = each %$stat){
	while(my($primer,$value2) = each %$value1){
	    my @array = split /\;/, $value2;
	    if($array[0] eq 'there'){
		my $read_seq = substr $seq->{$read}, $value->{$read}->{$primer}->{str}-1, $value->{$read}->{$primer}->{end}-$value->{$read}->{$primer}->{str}+1;
		$array[0] = ($read_seq eq $seq_DB->{$primer}->{seq})? 'match':'no-match';
		
		$stat->{$read}->{$primer} = join ';', @array;
	    }
	}
    }
    
    ##### frag match
    my $tag;
    if($stat->{R1}->{'UPTAG-frag'} ne 'absent' && $stat->{R2}->{'cUPTAG-frag'} ne 'absent'){
	my $frag   = $stat->{R1}->{'UPTAG-frag'};
	my $c_frag = $stat->{R2}->{'cUPTAG-frag'};
	$tag = &frag_match($frag,$c_frag);
	
    }
    $stat->{'merged-UPTAG'} = $tag || 'not-found';
    
    return $stat;
}


sub analyze_AD_lox_stat(){

    my $value = shift;
    my $R1    = shift;
    my $R2    = shift;
    my $seq   = shift;

    my $stat;
    $stat->{R1}->{'ADloxP-primer'}    = 'there';
    $stat->{R1}->{'cloxP'}            = 'absent';
    $stat->{R1}->{'ADD1-primer'}      = 'absent';
    $stat->{R1}->{'DNTAG-frag'}       = 'absent';
    
    $stat->{R2}->{'ADlox2272-primer'} = 'there';
    $stat->{R2}->{'clox2272'}         = 'absent';
    $stat->{R2}->{'ADD2-primer'}      = 'absent';
    $stat->{R2}->{'cDNTAG-frag'}      = 'absent';
    
    ##### Existence check; cloxP
    if($R1->{'cloxP'}){
	if($R1->{'ADloxP-primer'}->{end} < $R1->{'cloxP'}->{str}){
	    my $offset = $R1->{'cloxP'}->{str} - $R1->{'ADloxP-primer'}->{end} -1;
	    $stat->{R1}->{'cloxP'} = 'there;offset='.$offset;
	} 
    }
    
    ##### Existence check; ADD1-primer
    if($R1->{'cloxP'} && $R1->{'ADD1-primer'}){
	if($R1->{'cloxP'}->{end} < $R1->{'ADD1-primer'}->{str}){
	    my $offset = $R1->{'ADD1-primer'}->{str} - $R1->{'cloxP'}->{end} -1;
	    $stat->{R1}->{'ADD1-primer'} = 'there;offset='.$offset;
	    $stat->{R1}->{'DNTAG-frag'}  = substr $seq->{R1},$R1->{'ADD1-primer'}->{end},length($seq->{R1}) - $R1->{'ADD1-primer'}->{end};
	} 
    }
    
    ##### Existence check; clox2272
    if($R2->{'clox2272'}){
	if($R2->{'ADlox2272-primer'}->{end} < $R2->{'clox2272'}->{str}){
	    my $offset = $R2->{'clox2272'}->{str} - $R2->{'ADlox2272-primer'}->{end} -1;
	    $stat->{R2}->{'clox2272'} = 'there;offset='.$offset;
	} 
    }
    
    ##### Existence check; ADD2-primer
    if($R2->{'clox2272'} && $R2->{'ADD2-primer'}){
	if($R2->{'clox2272'}->{end} < $R2->{'ADD2-primer'}->{str}){
	    my $offset = $R2->{'ADD2-primer'}->{str} - $R2->{'clox2272'}->{end} -1;
	    $stat->{R2}->{'ADD2-primer'} = 'there;offset='.$offset;
	    $stat->{R2}->{'cDNTAG-frag'}  = substr $seq->{R2},$R2->{'ADD2-primer'}->{end},length($seq->{R2}) - $R2->{'ADD2-primer'}->{end};
	} 
    }
    
    ##### Sequence match check
    while(my($read,$value1) = each %$stat){
	while(my($primer,$value2) = each %$value1){
	    my @array = split /\;/, $value2;
	    if($array[0] eq 'there'){
		my $read_seq = substr $seq->{$read}, $value->{$read}->{$primer}->{str}-1, $value->{$read}->{$primer}->{end}-$value->{$read}->{$primer}->{str}+1;
		$array[0] = ($read_seq eq $seq_DB->{$primer}->{seq})? 'match':'no-match';
		
		$stat->{$read}->{$primer} = join ';', @array;
	    }
	}
    }
    
    ##### frag match
    my $tag;
    if($stat->{R1}->{'DNTAG-frag'} ne 'absent' && $stat->{R2}->{'cDNTAG-frag'} ne 'absent'){
	my $frag   = $stat->{R1}->{'DNTAG-frag'};
	my $c_frag = $stat->{R2}->{'cDNTAG-frag'};
	$tag = &frag_match($frag,$c_frag);
    }
    $stat->{'merged-DNTAG'} = $tag || 'not-found';
    

    return $stat;
}




sub frag_match(){
    my $frag1  = shift;
    my $c_frag = shift;
    my $tag = '';

    my $frag2 = reverse $c_frag;
    $frag2 =~ tr/atgcATGC/tacgTACG/;
    

    #print "$frag1\n$frag2\n";

    my $rel;
    for my $i (0..length($frag1)-1){
	my $n1 = substr $frag1,$i,1;
	for my $j (0..length($frag2)-1){
	    my $n2 = substr $frag2,$j,1;
	    if($n1 eq $n2){
		my $diff = $i-$j;
		$rel->{$diff} ||= 0;
		$rel->{$diff}++;		
	    }
	}
    }

    my @array;
    while(my($diff,$count) = each %$rel){
	my $element;
	$element->{diff}  = $diff;
	$element->{count} = $count;
	push @array,$element;
    }
    
    @array = sort{$b->{count} <=> $a->{count}}@array;

    if($array[0]->{count}>=18){ ##### Think about this threshold
	my $diff  = $array[0]->{diff};
	my $count = $array[0]->{count};

	for my $k1 (0..length($frag2)-1+$diff){
	    my $k2 = $k1-$diff;
	    my $n1 = (0<=$k1 && $k1<length($frag1))? substr $frag1,$k1,1 : 'N';
	    my $n2 = (0<=$k2 && $k2<length($frag2))? substr $frag2,$k2,1 : 'N';

	    my $n = $n1;
	    if($n1 eq 'N'){
		$n = $n2;
	    }elsif($n1 ne 'N' && $n2 ne 'N' && $n1 ne $n2){
		$n = 'N';
	    }
	    $tag .= $n;
	}
	#print "$diff\n$count\n";
    }

    return $tag;
}


sub assign_category(){

    my $value = shift;
    my $R1    = shift;
    my $R2    = shift;
    my $seq   = shift;

    # PS1.0 and PS2.0 at the reasonable positions? (end at <40bp)
    my $goto1 = 0;


    my $P_SEQ1 = 0;
    my $P_SEQ2 = 0;

    if($R1->{'PS1.0-primer'} && $R2->{'PS2.0-primer'}){
	if($R1->{'PS1.0-primer'}->{end} < 40 && $R2->{'PS2.0-primer'}->{end} < 40){
	    $P_SEQ1 = substr $seq->{R1},$R1->{'PS1.0-primer'}->{str}-10,9;
	    $P_SEQ2 = substr $seq->{R2},$R2->{'PS2.0-primer'}->{str}-10,9;
	    $goto1 = 1;
	}
    }	

    # Row and Column priming sites from the reasonable positions? (start at <50bp)
    # Which category? DB-BC / DB-lox / AD-BC / AD-lox?
    my $DB_BC  = 0;
    my $DB_lox = 0;
    my $AD_BC  = 0;
    my $AD_lox = 0;
    my $R_seq;
    my $C_seq;

    if($goto1){
	#DB-BC
	if($R1->{'DBU1-primer'} && $R2->{'DBD2-primer'}){
	    if($R1->{'PS1.0-primer'}->{end} < $R1->{'DBU1-primer'}->{str} && $R1->{'DBU1-primer'}->{str} < 50 &&
	       $R2->{'PS2.0-primer'}->{end} < $R2->{'DBD2-primer'}->{str} && $R2->{'DBD2-primer'}->{str} < 50){

		$R_seq->{'DB-BC'} = substr $seq->{R1},$R1->{'PS1.0-primer'}->{end},$R1->{'DBU1-primer'}->{str}-$R1->{'PS1.0-primer'}->{end}+1-2;
		$C_seq->{'DB-BC'} = substr $seq->{R2},$R2->{'PS2.0-primer'}->{end},$R2->{'DBD2-primer'}->{str}-$R2->{'PS2.0-primer'}->{end}+1-2;
		$DB_BC = 1;
	    }
	}
	#DB-lox
	if($R1->{'DBloxP-primer'} && $R2->{'DBlox2272-primer'}){
	    if($R1->{'PS1.0-primer'}->{end} < $R1->{'DBloxP-primer'}->{str} && $R1->{'DBloxP-primer'}->{str} < 50 &&
	       $R2->{'PS2.0-primer'}->{end} < $R2->{'DBlox2272-primer'}->{str} && $R2->{'DBlox2272-primer'}->{str} < 50){

		$R_seq->{'DB-lox'} = substr $seq->{R1},$R1->{'PS1.0-primer'}->{end},$R1->{'DBloxP-primer'}->{str}-$R1->{'PS1.0-primer'}->{end}+1-2;
		$C_seq->{'DB-lox'} = substr $seq->{R2},$R2->{'PS2.0-primer'}->{end},$R2->{'DBlox2272-primer'}->{str}-$R2->{'PS2.0-primer'}->{end}+1-2;
		$DB_lox = 1;
	    }
	}
	#AD-BC
	if($R1->{'ADU1-primer'} && $R2->{'ADD2-primer'}){
	    if($R1->{'PS1.0-primer'}->{end} < $R1->{'ADU1-primer'}->{str} && $R1->{'ADU1-primer'}->{str} < 50 &&
	       $R2->{'PS2.0-primer'}->{end} < $R2->{'ADD2-primer'}->{str} && $R2->{'ADD2-primer'}->{str} < 50){

		$R_seq->{'AD-BC'} = substr $seq->{R1},$R1->{'PS1.0-primer'}->{end},$R1->{'ADU1-primer'}->{str}-$R1->{'PS1.0-primer'}->{end}+1-2;
		$C_seq->{'AD-BC'} = substr $seq->{R2},$R2->{'PS2.0-primer'}->{end},$R2->{'ADD2-primer'}->{str}-$R2->{'PS2.0-primer'}->{end}+1-2;
		$AD_BC = 1;
	    }
	}
	#AD-lox
	if($R1->{'ADloxP-primer'} && $R2->{'ADlox2272-primer'}){
	    if($R1->{'PS1.0-primer'}->{end} < $R1->{'ADloxP-primer'}->{str} && $R1->{'ADloxP-primer'}->{str} < 50 &&
	       $R2->{'PS2.0-primer'}->{end} < $R2->{'ADlox2272-primer'}->{str} && $R2->{'ADlox2272-primer'}->{str} < 50){

		$R_seq->{'AD-lox'} = substr $seq->{R1},$R1->{'PS1.0-primer'}->{end},$R1->{'ADloxP-primer'}->{str}-$R1->{'PS1.0-primer'}->{end}+1-2;
		$C_seq->{'AD-lox'} = substr $seq->{R2},$R2->{'PS2.0-primer'}->{end},$R2->{'ADlox2272-primer'}->{str}-$R2->{'PS2.0-primer'}->{end}+1-2;
		$AD_lox = 1;
	    }
	}
    }
    
    # Specific category assignment?
    my $category = 0;
    my $R_SEQ    = 0;
    my $C_SEQ    = 0;

    if($DB_BC + $DB_lox + $AD_BC + $AD_lox == 1){
	if($DB_BC){
	    $category = 'DB-BC';
	}elsif($DB_lox){
	    $category = 'DB-lox';
	}elsif($AD_BC){
	    $category = 'AD-BC';
	}elsif($AD_lox){
	    $category = 'AD-lox';
	}

	$R_SEQ = $R_seq->{$category};
	$C_SEQ = $C_seq->{$category};
    }


    return ($category,$P_SEQ1,$P_SEQ2,$R_SEQ,$C_SEQ);
}



sub barcode_matching(){
    my $bar2num = shift;
    my $seq     = shift;

    my @bar_count;
    for my $key (keys %$bar2num){

	my $count;
	for(my $k=0; $k<length($key)-1; $k++){
	    my $k_nuc = substr $key,$k,2;
	    for(my $i=0; $i<length($seq)-1; $i++){
		my $i_nuc = substr $seq,$i,2;
		
		my $rel = $k-$i;
		if($k_nuc eq $i_nuc){
		    $count->{$rel} ||= 0;
		    $count->{$rel}  += 1;
		}else{
		    $count->{$rel} ||= 0;
		    $count->{$rel}  += 0;
		}
	    }
	}

	my @array = sort{$b <=> $a} values %$count;
	my $max_count = shift @array;
	my $element;
	$element->{bar} = $key;
	$element->{max} = $max_count;
	push @bar_count, $element;
    }
    my $val;
    @bar_count = sort{$b->{max} <=> $a->{max}} @bar_count;
    if($bar_count[0]->{max} == $bar_count[1]->{max}){
	$val = 0;
    }elsif($bar_count[0]->{max}>5){ ########## THINK ABOUT THIS THRESHOLD
	my $called_bar = $bar_count[0]->{bar};
	$val = $bar2num->{$called_bar};
    }else{
	$val = 0;
    }

    return $val;
}






sub complement(){
    my $seq = shift;
    $seq = reverse $seq;
    $seq =~ tr/ATGCatgc/TACGtacg/;
    return $seq;
}

sub get_bar2num(){
    my $file = shift;
    my $bar2num;
    open FILE,$file;
    while(<FILE>){
	chomp;
	my @array = split /\,/,$_;
	if($#array){
	    $bar2num->{$array[1]} = $array[0];
	}
    }
    close FILE;

    return $bar2num;
}


sub get_tag(){
    my $file = shift;
    my $tag_strct;
    open FILE,$file;
    while(<FILE>){
	chomp;
	my @array = split /\,/,$_;
	if($#array == 3){
	    $tag_strct->{$array[1]}->{$array[2]}->{$array[3]} = $array[0];
	}
    }
    close FILE;

    return $tag_strct;
}

sub get_seq(){
    my @files = @_;

    my $seq;
    for my $file (@files){

	my $tag = 'N_A';
	open FILE, $file;
	while(<FILE>){
	    chomp;
	    if($_ =~ /\>(.+)$/){
		$tag = $1;
	    }else{
		$seq->{$tag}->{seq} .= $_;
		$seq->{$tag}->{seq} =~ s/[^a-zA-Z]//g;
	    }
	}
	close FILE;
	
	while(my($tag,$value) = each %$seq){
	    $seq->{$tag}->{length} = length $seq->{$tag}->{seq}; 
	}
    }

    return $seq;
}
