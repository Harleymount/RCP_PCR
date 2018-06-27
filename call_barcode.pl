#!/usr/bin/perl -w

use strict;
use Data::Dumper;

my $count = do shift;


my $data2;

for my $P_TAG (sort{$a cmp $b}(keys %$count)){

    my $type = 'N_A';
    $type = 'DB' if($P_TAG =~ /^DB/);
    $type = 'AD' if($P_TAG =~ /^AD/);

    #print "$P_TAG\n";
    for my $R_TAG (1..16){
	$R_TAG = sprintf "%.1f", $R_TAG/10;
	$R_TAG =~ s/\.//g;
	$R_TAG = "R$R_TAG";
	
	for my $C_TAG (1..24){
	    $C_TAG = sprintf "%.1f", $C_TAG/10;
	    $C_TAG =~ s/\.//g;
	    $C_TAG = "C$C_TAG";
	    
	    #print "$type\t$P_TAG\t$R_TAG\t$C_TAG\n";		

	    my $rate_strct = &call_barcode($count->{$P_TAG}->{$R_TAG}->{$C_TAG},$type);
	    
	    my $BC_all  = $rate_strct->{BC_all};
	    my $lox_all = $rate_strct->{lox_all};
	    
	    #print "$BC_all\n";
	    #print "$lox_all\n";
	    
	    #print Dumper $rate_strct;
	    
	    $data2->{$type}->{$P_TAG}->{$R_TAG}->{$C_TAG}->{BC_all}  = $BC_all;
	    $data2->{$type}->{$P_TAG}->{$R_TAG}->{$C_TAG}->{lox_all} = $lox_all;
	    
	    for my $element (@{$rate_strct->{array}}){
		my $evaluate = &evaluate($element,$type);
		
		push @{$data2->{$type}->{$P_TAG}->{$R_TAG}->{$C_TAG}->{barcode_called}}, $evaluate if($evaluate);
		
		#print Dumper $evaluate if($evaluate);
	    }
	    
	   
	}
    }
}



print Dumper $data2;




################################################
# function 
################################################

sub evaluate(){
    my $element = shift;
    my $type    = shift;

    my $evaluate;

    my @primer_order = ('U1','U2','D1','D2');
    my @lox_order    = ('loxP','lox2272');


    if($element->{BC}->{'min_BC_hit_rate_impt'} > 0.1){ ##### THINK ABOUT THIS THRESHOLD

	my $primer_stat;
	my $lox_stat;
	
	my @primer_array = split /\:/, $element->{BC}->{'likely_stat'};  @primer_array = (@primer_array,'N_A','N_A','N_A','N_A');
	my @lox_array    = split /\:/, $element->{lox}->{'likely_stat'}; @lox_array    = (@lox_array,   'N_A','N_A','N_A','N_A');

	for my $i (0..3){
	    if($element->{BC}->{'stat_rate'} > 0.5 && $primer_array[$i] =~ /^match/){
		$primer_stat->{$primer_order[$i]} = '+';
	    }elsif($element->{BC}->{'stat_rate'} > 0.5 && $primer_array[$i] !~ /^match/){
		$primer_stat->{$primer_order[$i]} = '-';
	    }else{
		$primer_stat->{$primer_order[$i]} = 'Ab';
	    }
	}

	for my $i (0..1){
	    if($element->{lox}->{'stat_rate'} > 0.5 && $lox_array[$i] =~ /^match/){
		$lox_stat->{$lox_order[$i]} = '+';
	    }elsif($element->{lox}->{'stat_rate'} > 0.5 && $lox_array[$i] !~ /^match/){
		$lox_stat->{$lox_order[$i]} = '-';
	    }else{
		$lox_stat->{$lox_order[$i]} = 'Ab';
	    }
	}

	$evaluate = $element;
	$evaluate->{BC}->{likely_stat}  = $primer_stat;
	$evaluate->{lox}->{likely_stat} = $lox_stat;
    }

    return $evaluate;
}


sub call_barcode(){
    my $data = shift;
    my $type = shift;

    my $count_all;

    my $BC_all  = 0;
    my $lox_all = 0;

    ####################################
    # Count useful BC tag and its stat #
    ####################################
    my @count_array;
    while(my($tags,$value1) = each %{$data->{"$type\-BC"}}){
        if($tags =~ /^[ATGCN]+\:[ATGCN]+$/){

            my $element;
            $element->{tags} = $tags;

            while(my($stat,$count) = each %$value1){
		$element->{count}->{All}   += $count;
		$element->{count}->{$stat} += $count;
		$BC_all+=$count;
            }

            push @count_array, $element;
        }
    }
    @count_array = sort{$b->{count}->{All} <=> $a->{count}->{All}} @count_array;


    #######################################
    # Count useful lox reads and its stat #
    #######################################
    while(my($tag,$value1) = each %{$data->{"$type\-lox"}}){
	if($tag =~ /^[ATGCN]+$/){
	    while(my($stat,$count) = each %$value1){
		$lox_all+=$count;
	    }
	}
    }

    $count_all->{BC_all}  = $BC_all;
    $count_all->{lox_all} = $lox_all;

    #print "$BC_all\n";
    #print "$lox_all\n";


    ########################
    ##### impute count #####
    ########################
    for my $i (0..$#count_array){
	for my $j ($i+1..$#count_array){

	    my $tags_i = $count_array[$i]->{tags};
	    my $tags_j = $count_array[$j]->{tags};
	    my $all_i  = $count_array[$i]->{count}->{All};
	    my $all_j  = $count_array[$j]->{count}->{All};

	    if($all_i > 0 && $all_i > $all_j){

		my $match = &match($tags_i,$tags_j);

		if($match->{rate}>0.9){
		    for my $stat (keys %{$count_array[$j]->{count}}){
			$count_array[$i]->{count}->{$stat} ||= 0;
			$count_array[$i]->{count}->{$stat}  += $count_array[$j]->{count}->{$stat};
			$count_array[$j]->{count}->{$stat}   = 0;
		    }
		    
		}    
	    }
	}
    }
    @count_array = sort{$b->{count}->{All} <=> $a->{count}->{All}} @count_array;

    for my $element (@count_array){
	if($element->{count}->{All}>0){
	    my ($UPTAG,$DNTAG) = split /\:/, $element->{tags};

	    my $lox;
	    while(my($tag,$value1) = each %{$data->{"$type\-lox"}}){
		if($tag =~ /^[ATGCN]+$/){
		    my $match = 0;
		    $match = &match($UPTAG,$tag) if($type eq 'DB');
		    $match = &match($DNTAG,$tag) if($type eq 'AD');

		    while(my($stat,$count) = each %$value1){
			if($match->{rate}>0.8){
			    $lox->{count}->{All}   += $count;
			    $lox->{count}->{$stat} += $count;
			}
		    }	    
		}
	    }	    

	    # likely BC stat
	    my @stat_BC;
	    while(my($stat,$count) = each %{$element->{count}}){
		if($stat ne 'All'){
		    my $i;
		    $i->{stat}  = $stat;
		    $i->{count} = $count;
		    push @stat_BC,$i;
		}
	    }
	    @stat_BC = sort{$b->{count} <=> $a->{count}} @stat_BC;
	    my $likely_BC_stat       = $stat_BC[0]->{stat}  || 'N_A';
	    my $likely_BC_stat_count = $stat_BC[0]->{count} || 0;


	    # likely lox stat
	    my @stat_lox;
	    while(my($stat,$count) = each %{$lox->{count}}){
		if($stat ne 'All'){
		    my $i;
		    $i->{stat}  = $stat;
		    $i->{count} = $count;
		    push @stat_lox,$i;
		}
	    }
	    @stat_lox = sort{$b->{count} <=> $a->{count}} @stat_lox;
	    my $likely_lox_stat       = $stat_lox[0]->{stat}  || 'N_A';
	    my $likely_lox_stat_count = $stat_lox[0]->{count} || 0;


	    #############################
	    # counts need for screening #
	    #############################
	    
	    my $count_strct;
	    $count_strct->{BC}->{UPTAG}              = $UPTAG || 'N_A';
	    $count_strct->{BC}->{DNTAG}              = $DNTAG || 'N_A';
	    $count_strct->{BC}->{hit}                = $element->{count}->{All} || 0;
	    $count_strct->{BC}->{likely_stat}        = $likely_BC_stat;
	    $count_strct->{BC}->{likely_stat_count}  = $likely_BC_stat_count;
	    $count_strct->{lox}->{hit}               = $lox->{count}->{All}     || 0;
	    $count_strct->{lox}->{likely_stat}       = $likely_lox_stat;
	    $count_strct->{lox}->{likely_stat_count} = $likely_lox_stat_count;

	    push @{$count_all->{array}}, $count_strct;

	}
    }

    my $rate_strct = &cal_rate($count_all,$type);

    return $rate_strct;
}


sub cal_rate(){
    my $count_all = shift;
    my $type      = shift;
    
    my $BC_all  = $count_all->{BC_all};
    my $lox_all = $count_all->{lox_all};

    if($count_all->{array}){
	my @array = @{$count_all->{array}};


	my $hit_BC_all = 0;
	
	for(my $i = $#array; $i>=0; $i--){
	    
	    my $untreated_BC_all = 0;
	    for(my $j = $i-1; $j>=0; $j--){
		$untreated_BC_all += $array[$j]->{BC}->{hit};
	    }
	    
	    my ($BC_hit_rate, $BC_hit_rate_error) = &ErrorProp_div($array[$i]->{BC}->{hit}, $BC_all-$untreated_BC_all-$hit_BC_all);
	    my $min_BC_hit_rate   = $BC_hit_rate   - $BC_hit_rate_error;
	    
	    $array[$i]->{BC}->{'min_BC_hit_rate_impt'} = $min_BC_hit_rate;
	    
	    my ($BC_stat_rate, $BC_stat_rate_error)  = &ErrorProp_div($array[$i]->{BC}->{likely_stat_count},$array[$i]->{BC}->{hit});
	    my ($lox_stat_rate,$lox_stat_rate_error) = &ErrorProp_div($array[$i]->{lox}->{likely_stat_count},$array[$i]->{lox}->{hit});
	    
	    my $min_BC_stat_rate  = $BC_stat_rate  - $BC_stat_rate_error;
	    my $min_lox_stat_rate = $lox_stat_rate - $lox_stat_rate_error;
	    
	    $array[$i]->{BC}->{stat_rate}  = $min_BC_stat_rate;
	    $array[$i]->{lox}->{stat_rate} = $min_lox_stat_rate;
	    
	}

	$count_all->{array} = \@array;
    }

    return $count_all;
}


sub ErrorProp_div(){
    my $value1 = shift || 0;
    my $value2 = shift || 0;

    my $value = 0;
    my $error = 0;

    if($value1 && $value2){
	$value = $value1/$value2;
	my $d_value1 = $value1**(1/2); 
	my $d_value2 = $value2**(1/2);
	
	$error = ($d_value1/$value1) ** 2 + ($d_value2/$value2) ** 2;
	$error **= 1/2;
	$error  *= $value;
    }

    return ($value,$error);
}



sub match(){
    my $tags_i = shift;
    my $tags_j = shift;

    my $match;

    my $min_len = length($tags_i);
    if(length($tags_j)<$min_len){
	$min_len = length($tags_j);
    }

    my $num = 0;
    for my $k (0..$min_len-1){
	my $ni = substr $tags_i,$k,1;
	my $nj = substr $tags_j,$k,1;
	$num++ if($ni eq $nj);
    }

    $match->{num}  = $num;
    $match->{rate} = $num/$min_len;

    #print "$tags_i\t$tags_j\n";
    #print Dumper $match;


    return $match;
}
