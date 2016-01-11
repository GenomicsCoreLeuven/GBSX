#!/usr/bin/perl

####################################################################################################################
# This is GBSX v1.0. A toolkit for experimental design and demultiplexing genotyping by sequencing experiments.    #
#                                                                                                                  #
#  Copyright (C) 2014 KU Leuven                                                                                    #
#                                                                                                                  #
#  This file is part of GBSX.                                                                                      #
#                                                                                                                  #
#  GBSX is free software: you can redistribute it and/or modify                                                    #
#  it under the terms of the GNU General Public License as published by                                            #
#  the Free Software Foundation, either version 3 of the License, or                                               #
#  (at your option) any later version.                                                                             #
#                                                                                                                  #
#  GBSX is distributed in the hope that it will be useful,                                                         #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of                                                  #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                   #
#  GNU General Public License for more details.                                                                    #
#                                                                                                                  #
#  You should have received a copy of the GNU General Public License                                               #
#  along with GBSX.  If not, see <http://www.gnu.org/licenses/>.                                                   #
####################################################################################################################

########################################################################################################################################################################
use Getopt::Std;
use Bio::SeqIO;
use Data::Dumper;
use Bio::PrimarySeq;
use Bio::DB::Fasta;
use Bio::Restriction::EnzymeCollection;
use Bio::Restriction::Analysis;
use Bio::Restriction::Enzyme;



########################################################################################################################################################################
my $program_information = "GBSX v1.0";
my $license = "GBSX v1.0  Copyright (C) 2014 KU Leuven\nThis program comes with ABSOLUTELY NO WARRANTY; for details type `perl GBSX_digest_v1.0.pl -w'.\nThis is free software, and you are welcome to redistribute it under certain conditions; type `perl GBSX_digest_v1.0.pl -c' for details.\n";

$Getopt::Std::STANDARD_HELP_VERSION = 1;
my $help_message="\n$license\n####################################################################################\nusage: perl GBSX_digest_v1.0.pl -d \'digest sequence\' -l \'read length\' -f \'file of reference fasta file location(s)\'\n\noptional parameters: \n\t-e \'enzyme name to use\' (default: Enzyme)\n\t-g \'genome name to use in bed file name\' (default: genome)\n\t-n \'minimum size fragments to include\' (default: 100)\n\t-m \'maximum size fragments to use\' (default: 1000)\n\n\t-E \'second enzyme name to use\' (default: Enzyme2)\n\t-D \'digest sequence for a second enzyme\' (default: not declared)\n\t-R \'digest sequence for a third enzyme\' (default: not declared)\n";

########################################################################################################################################################################
########################################################################################################################################################################
#input parameters
sub HELP_MESSAGE {die "$help_message\n";}
sub VERSION_MESSAGE {print "$program_information\n";}

my %opts = ();
getopts('e:d:g:l:f:m:n:E:D:R:cwh');

#die with help message if -h 
if ($opt_h){ die "$help_message\n"; }
#die with warrenty message if -w 
if ($opt_w){ die "GBSX is distributed in the hope that it will be useful,\nbut WITHOUT ANY WARRANTY; without even the implied warranty of\nMERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\nGNU General Public License for more details.\n"; }
#die with conditions message if -c 
if ($opt_c){ die "GBSX is free software: you can redistribute it and/or modify\nit under the terms of the GNU General Public License as published by\nthe Free Software Foundation, either version 3 of the License, or\n(at your option) any later version.\n"; }


#enzyme name (optional)
my $enzyme_name = $opt_e || "Enzyme";

#digest sequence: example "RCATG^Y";
my $digest_seq = $opt_d || "not_defined"; 	if ($digest_seq ne "not_defined"){$digest_seq=uc($digest_seq);}

#genome name to use for bedfile (optional)
my $bed_genome_name = $opt_g || "genome";

#define read length: assumes PE, only report back in bed first and last number of bps corresponding to this
my $read_length = $opt_l || "not_defined";

#file of reference fasta locations
my $fasta_file_locations = $opt_f || "not_defined";

#max/min fragment sizes (optional, otherwise use defaults)
my $max_size = $opt_m || 1000;
my $min_size = $opt_n || 100; 

#enzyme2 name (optional, otherwise if second digest just default to Enzyme2)
my $enzyme2_name = $opt_E || "Enzyme2";
	#make sure starts Upper case
	$enzyme2_name = ucfirst($enzyme2_name);
	#double check this is not the same as Enzyme1
	if ($enzyme2_name eq $enzyme_name){die "Please use different enzyme names for enzyme 1 and enzyme 2\n"; }

#digest sequence for a second enzyme: example "RCATG^Y"; (optional, otherwise no second digest)
my $digest2_seq = $opt_D || "not_defined"; 	if ($digest2_seq ne "not_defined"){$digest2_seq=uc($digest2_seq);}
	#double check this is not the same as Enzyme1
	if ($digest2_seq eq $digest_seq){die "Please use different enzyme cut sequences for enzyme 1 and enzyme 2\n"; }

#digest sequence for a third enzyme: example "RCATG^Y"; (optional, otherwise no third digest)
my $digest3_seq = $opt_R || "not_defined"; 	
	if ($digest3_seq ne "not_defined"){$digest3_seq=uc($digest3_seq);
		#double check this is not the same as Enzyme1 or Enzyme3
		if ( ($digest3_seq eq $digest2_seq) || ($digest3_seq eq $digest_seq) ){die "Please use different enzyme cut sequences for enzyme 1, enzyme 2, and enzyme 3\n"; }
		#we call this Enzyme3, so make sure nothing else uses this name
		if ( ($enzyme2_name eq "Enzyme3") || ($enzyme_name eq "Enzyme3") ){ die "Enzyme3 is a reseved name, please use a different enzyme name\n"; }
	}

#die if required parameters are missing
if ( ($digest_seq eq "not_defined") || ($read_length eq "not_defined") || ($fasta_file_locations eq "not_defined") ){
	my $missing=" ";
	if ($digest_seq eq "not_defined"){$missing = $missing."-d (digest sequence) ";}
	if ($read_length eq "not_defined"){$missing = $missing."-l (read length) ";}
	if ($fasta_file_locations eq "not_defined"){$missing = $missing."-f (reference fasta locations) ";}
	die "missing parameter(s):$missing\n";
}

########################################################################################################################################################################
########################################################################################################################################################################
#first basic steps to get ready for the analysis (mostly defining variables)

#using input information, define outfiles:
my $outfile;	my $main_bed;
if ($digest2_seq ne "not_defined"){
	$outfile = $bed_genome_name.".".$enzyme_name.".".$enzyme2_name.".".$read_length."nt.digest_results";
	$main_bed = $bed_genome_name.".".$enzyme_name.".".$enzyme2_name.".".$read_length."nt.digest.bed";
}else{
	$outfile = $bed_genome_name.".".$enzyme_name.".".$read_length."nt.digest_results";
	$main_bed = $bed_genome_name.".".$enzyme_name.".".$read_length."nt.digest.bed";
}

#get reference fasta files
my @ref_chr;
open (ref_files,$fasta_file_locations)||die "could not open file $fasta_file_locations\n";
while(<ref_files>){
    #update 1.1: change chop to chomp
	#chop;
    chomp;
       	$_=~s/\n//;$_=~s/\r//;
	if ($_ ne ""){ 
		push(@ref_chr,$_); 
		#current script only works with one sequence per file, add a check for this
		my $grep_results = `grep \"^\>\" $_`;
		my @grep_array = split(/\n/,$grep_results); my $number_sequences_in_file = @grep_array;
		if ($number_sequences_in_file > 1){die "please use only one sequence per fasta file, $_ had $number_sequences_in_file sequences\n";}
	}
}
close(ref_files)||die "could not close file $fasta_file_locations\n";

#define enzyme to use
my $rest_enzyme = new Bio::Restriction::Enzyme(-enzyme => $enzyme_name, -seq => $digest_seq);

#optional: define a second enzyme
my $rest_enzyme2; if ($digest2_seq ne "not_defined"){ $rest_enzyme2 = new Bio::Restriction::Enzyme(-enzyme => $enzyme2_name, -seq => $digest2_seq); }

#optional: define a third enzyme
my $rest_enzyme3; if ($digest3_seq ne "not_defined"){ $rest_enzyme3 = new Bio::Restriction::Enzyme(-enzyme => "Enzyme3", -seq => $digest3_seq); }

#define ouput count variables
my $total_cuts=0; my $total_fragments=0; my $total_fragments_between_min_and_max = 0;
my $bin_size=100;#hard coded for now
my $temp_size=$min_size;
	until (($temp_size+$bin_size)>$max_size){
		my $temp_size_plus_bin = $temp_size+$bin_size;
		my $variable_name = "total_fragments_between_".$temp_size."_and_".$temp_size_plus_bin;
		$$variable_name=0;
		$temp_size=$temp_size+$bin_size;
	}
my $total_fragments_less_min = 0;
my %fragments_in_range_per_chr; my %fragments_in_range_per_chr_enzyme1_end; my %fragments_in_range_per_chr_enzyme2_end; my %fragments_in_range_per_chr_enzymeboth_end;

#define ouput count variables, optional for a second enzyme
my $total_cuts_enzyme2=0; my $fragments_both_ends_enzyme1 = 0; my $fragments_both_ends_enzyme2 = 0; my $fragments_end_from_each_enzyme = 0;

#define ouput count variables, optional for a third enzyme
my $fragments_in_range_contain_third_enzyme = 0;

#define ouput count variables, distribution of distances
my $one_kb = 0; my $ten_kb = 0; my $hundred_kb=0; my $one_Mb = 0; my $ten_Mb = 0;	my $more_ten_kb=0;

########################################################################################################################################################################
########################################################################################################################################################################
#per chromosome, perform digest and add to totals

#open bed file to print to
open(outbed,">$main_bed")||die "could not open out bedfile $main_bed\n";

foreach my $chr (@ref_chr){

	#just get chr name and remove directory structure
	my @chr_name_parts =split(/\//,$chr);my $number_parts=@chr_name_parts;
	my $chr_name=@chr_name_parts[$number_parts-1];
	$chr_name=~s/.fa$//i;$chr_name=~s/.fasta$//i;
	$fragments_in_range_per_chr{$chr_name}=0;
	$fragments_in_range_per_chr_enzyme1_end{$chr_name}=0;
	$fragments_in_range_per_chr_enzyme2_end{$chr_name}=0;
	$fragments_in_range_per_chr_enzymeboth_end{$chr_name}=0;

	#get fasta file for this chr
	my $seqio_object1 = Bio::SeqIO->new(-file => $chr, -format => "fasta");
	my $seq = $seqio_object1->next_seq;
	my $fastadb = Bio::DB::Fasta->new($chr);

	#digest
	my $rest_analysis=Bio::Restriction::Analysis->new(-seq=>$seq, -enzymes=>$rest_enzyme);

	#just get cut numbers first
	my $chr_cuts=$rest_analysis->cuts_by_enzyme($enzyme_name);
	$total_cuts=$total_cuts+$chr_cuts;

	#optional digest and cut numbers for a second enzyme
	my $rest_analysis2;my $chr_cuts2; 
	if ($digest2_seq ne "not_defined"){
		$rest_analysis2=Bio::Restriction::Analysis->new(-seq=>$seq, -enzymes=>$rest_enzyme2);
		$chr_cuts2=$rest_analysis2->cuts_by_enzyme($enzyme2_name);
		$total_cuts_enzyme2=$total_cuts_enzyme2+$chr_cuts2;
	}

	#get locations of cuts
	my @locations1 = $rest_analysis->positions($enzyme_name);
	#cut locations if second enzyme
	my %hash_locations1; my %hash_locations2;#use position as keys
	my @locations2;my @all_locations;
	if ($digest2_seq ne "not_defined"){
		@locations2=$rest_analysis2->positions($enzyme2_name);
		@all_locations = (@locations1,@locations2);
		@all_locations=sort {$a <=> $b} @all_locations;

		foreach my $position1 (@locations1){ $hash_locations1{$position1}=1; }
		foreach my $position2 (@locations2){ $hash_locations2{$position2}=1; }
	}else{
		@all_locations=@locations1;
	}

	#go through all cut locations and check fragment sizes
	my @chr_with_directories = split(/\//,$chr);  $chr_name=pop @chr_with_directories;  $chr_name=~s/.fa//;
	my $previous_location=0;	my $previous_fragment_in_range_end=0;

	#to analyze every fragment, add the end of the chromosome as a position
	my $chr_end_position = $seq->length;
	push(@all_locations,$chr_end_position);

	foreach my $location (@all_locations){
		$total_fragments++;
		if ($location<$previous_location){die "error in order, $location < $previous_location";}
		my $fragment_length = $location-$previous_location;
		if (($fragment_length<=$max_size)&&($fragment_length>=$min_size)){
			$total_fragments_between_min_and_max++;
			$fragments_in_range_per_chr{$chr_name}=$fragments_in_range_per_chr{$chr_name}+1;

			#for distribution of distance between sequenced fragments
			if ($previous_fragment_in_range_end!=0){
				my $distance_this_fragment_to_previous=$previous_location-$previous_fragment_in_range_end;
				if ( $distance_this_fragment_to_previous<=1000 ){
					$one_kb++;
				}elsif( ($distance_this_fragment_to_previous>1000) && ($distance_this_fragment_to_previous<=10000)){
					$ten_kb++;
				}elsif( ($distance_this_fragment_to_previous>10000) && ($distance_this_fragment_to_previous<=100000)){
					$hundred_kb++;
				}elsif( ($distance_this_fragment_to_previous>100000) && ($distance_this_fragment_to_previous<=1000000)){
					$one_Mb++;
				}elsif( ($distance_this_fragment_to_previous>1000000) && ($distance_this_fragment_to_previous<=10000000)){
					$ten_Mb++;
				}elsif( ($distance_this_fragment_to_previous>10000000) ){
					$more_ten_kb++;
				}else{
					die "error in distances\n";
				}
			}
				$previous_fragment_in_range_end=$location;
			

			#if second enzyme, are ends from one or both enzymes
			if ($digest2_seq ne "not_defined"){
				if ( (exists($hash_locations1{$previous_location})) && (exists($hash_locations1{$location})) ){
					$fragments_both_ends_enzyme1++;	#both ends from enzyme1
					$fragments_in_range_per_chr_enzyme1_end{$chr_name}=$fragments_in_range_per_chr_enzyme1_end{$chr_name}+1;
				}elsif ( (exists($hash_locations2{$previous_location})) && (exists($hash_locations2{$location})) ){
					$fragments_both_ends_enzyme2++;	#both ends from enzyme2
					$fragments_in_range_per_chr_enzyme2_end{$chr_name}=$fragments_in_range_per_chr_enzyme2_end{$chr_name}+1;
				}elsif ( (exists($hash_locations1{$previous_location})) && (exists($hash_locations2{$location})) ){
					$fragments_end_from_each_enzyme++;	#one end from each enzyme
					$fragments_in_range_per_chr_enzymeboth_end{$chr_name}=$fragments_in_range_per_chr_enzymeboth_end{$chr_name}+1;
				}elsif ( (exists($hash_locations2{$previous_location})) && (exists($hash_locations1{$location})) ){
					$fragments_end_from_each_enzyme++;	#one end from each enzyme
					$fragments_in_range_per_chr_enzymeboth_end{$chr_name}=$fragments_in_range_per_chr_enzymeboth_end{$chr_name}+1;
				}
			}

			#counts on bin sizes
			my $temp_size=$min_size;
			until (($temp_size+$bin_size)>$max_size){
				my $temp_size_plus_bin = $temp_size+$bin_size;
				if ( ($fragment_length>=$temp_size)&&($fragment_length<$temp_size_plus_bin) ){
					my $variable_name = "total_fragments_between_".$temp_size."_and_".$temp_size_plus_bin;
					$$variable_name=$$variable_name+1;
				}
				$temp_size=$temp_size+$bin_size;
			}

			#print to bed the potentially sequenced bases
			if ($fragment_length <= (2*$read_length)){
				print outbed "$chr_name\t$previous_location\t$location\n";
			}else{
				#print first read
				$previous_location_plus=$previous_location+$read_length;
				print outbed "$chr_name\t$previous_location\t$previous_location_plus\n";
				#print second read
				$location_minus=$location-$read_length;
				print outbed "$chr_name\t$location_minus\t$location\n";
			}
			
			#if third enzyme requested, check if there is a digest site in this fragment
			if ($digest3_seq ne "not_defined"){
				my $fragment_seq = $fastadb->seq($chr_name, ($previous_location+1) => ($location));
				my $frag_seq_obj = Bio::Seq->new( -seq => $fragment_seq );
				my $rest_analysis3=Bio::Restriction::Analysis->new(-seq=>$frag_seq_obj, -enzymes=>$rest_enzyme3);
				my $frag_cuts3=$rest_analysis3->cuts_by_enzyme("Enzyme3");
				if ($frag_cuts3>0){ $fragments_in_range_contain_third_enzyme++ }
			}

		}elsif ($fragment_length<$min_size){
			$total_fragments_less_min++;
		}
		$previous_location=$location;
	}

}

#close bedfile
close(outbed)||die "could not open out bedfile $main_bed\n";

########################################################################################################################################################################
########################################################################################################################################################################
#print to output file general digest info
open(out,">$outfile")||die "could not open outfile $outfile\n";

#print input parameters
print out "input parameters:\n\tenzyme name used: $enzyme_name\n\tenzyme cut sequence: $digest_seq\n";
if ($digest2_seq ne "not_defined"){ print out "\tsecond enzyme name used: $enzyme2_name\n\tsecond enzyme cut sequence: $digest2_seq\n"; }
if ($digest3_seq ne "not_defined"){ print out "\tthird enzyme cut sequence: $digest3_seq\n"; }
print out "\tfor bedfile, genome name used: $bed_genome_name\n\tread length provided: $read_length\nminimum fragment size: $min_size\nmaximum fragment size: $max_size\n\nreference sequence(s) analyzed:\n";
foreach my $ref_to_print_to_output (@ref_chr){ print out "\t$ref_to_print_to_output\n"; }
print out "\n";

#print number of cuts
print out "total cuts for enzyme ($enzyme_name): $total_cuts\n";
if ($digest2_seq ne "not_defined"){ print out "total cuts for second enzyme ($enzyme2_name): $total_cuts_enzyme2\n"; }

#print number of fragments
print out "\nlooked at $total_fragments fragments in total\n";
print out "$total_fragments_between_min_and_max fragments <= $max_size nt and >= $min_size nt\n";
if ($digest3_seq ne "not_defined"){ print out "\t-containing one or more third digest site ($digest3_seq): $fragments_in_range_contain_third_enzyme\n"; }
if ($digest2_seq ne "not_defined"){
	print out "\t-with both ends from $enzyme_name: $fragments_both_ends_enzyme1\n";
	print out "\t-with both ends from $enzyme2_name: $fragments_both_ends_enzyme2\n";
	print out "\t-with an end from each enzyme: $fragments_end_from_each_enzyme\n";
	print out "\tfragments in range per reference sequence: (both ends from $enzyme_name / both ends from $enzyme2_name / an end from each enzyme)\n";
}else{
	print out "\tfragments in range per reference sequence:\n";
}

#print fragment chromosome distribution
my @chr_names_for_count_print = keys %fragments_in_range_per_chr;
foreach my $chr_name_for_count_print (@chr_names_for_count_print){
	if ($digest2_seq ne "not_defined"){
		print out "\t\t$chr_name_for_count_print had $fragments_in_range_per_chr{$chr_name_for_count_print} fragments ($fragments_in_range_per_chr_enzyme1_end{$chr_name_for_count_print}/$fragments_in_range_per_chr_enzyme2_end{$chr_name_for_count_print}/$fragments_in_range_per_chr_enzymeboth_end{$chr_name_for_count_print})\n";
	}else{
		print out "\t\t$chr_name_for_count_print had $fragments_in_range_per_chr{$chr_name_for_count_print} fragments\n";
	}
}

#print fragment distance distribution
print out "\tdistances between fragments in range:\n";
print out "\t\t<=1kb $one_kb fragment pairs\n\t\t1kb-10kb $ten_kb fragment pairs\n\t\t10kb-100kb $hundred_kb fragment pairs\n\t\t100kb-1Mb $one_Mb fragment pairs\n";
print out "\t\t1Mb-10Mb $ten_Mb fragment pairs\n\t\t>10Mb $more_ten_kb fragment pairs\n";


#print fragment size distribution
print out "\tfragments less than $min_size nt: $total_fragments_less_min \n\n";
my $temp_size=$min_size;
until (($temp_size+$bin_size)>$max_size){
	my $temp_size_plus_bin = $temp_size+$bin_size;
	my $variable_name = "total_fragments_between_".$temp_size."_and_".$temp_size_plus_bin;
	$temp_size_plus_bin--;
	print out "\tfragments $temp_size-$temp_size_plus_bin nt: $$variable_name\n";
	$temp_size=$temp_size+$bin_size;
}
print out "\n";

close(out)||die "could not close outfile $outfile\n";
