#!/bin/bash

#this code only works on mac OSX


#when running this code will make the directory structure in the pre-made output directory folder
#in the pre-made output folder you will need to add plate_adapter fwed and rev fas 
#also need to have all of the perl scripts in this directory in a folder called scripts

#need to add a user input for the tag assignment csv file



for arg in "$@"; do
  shift
  case "$arg" in
    "--output_directory") set -- "$@" "-r" ;;
    "--read1") set -- "$@" "-x" ;;
    "--read2")   set -- "$@" "-y" ;;
	"--tags")   set -- "$@" "-z" ;; 
    *)        set -- "$@" "$arg"
  esac
done


while getopts r:x:y:z: option
do
 case "${option}"
 in
 r) output=$OPTARG;;
 x) r1=${OPTARG};;
 y) r2=${OPTARG};;
 z) tags=${OPTARG};; 
 esac
done


#----------------Get RCP programs from github

cd ~;
git clone https://github.com/Harleymount/Scientific-Code.git;
cd ~/Scientific-Code;
tar -xzf RCP_PCR_V2.tar.gz;
echo 'RCP Programs downloaded!';
cd ~/Desktop;
#RCP scripts found in ~/Scientific-Code/RCP-PCR/


#---------Make directory structure in output folder 
mkdir output;
cd output;
mkdir fragmented_fasta;
mkdir QC_and_BC;
cd QC_and_BC;
mkdir out.identification;
mkdir out.count;
mkdir sh.identification;
mkdir sh.count;
cd ..;
mkdir blast;
cd blast;
mkdir out.primers_blast;
mkdir sh.primers_blast;
cd ..;
mkdir scripts;
cp -r  ~/Scientific-Code/ scripts;
echo 'Directory structure created!';

#-------------remove PhiX reads and other junk -------------------------------------------
flexbar -r $r1 -p $r2 -b scripts/plate_adapter_fwd.fa -b2 scripts/plate_adapter_rev.fa -bk -u 8 -n 8;
#--------------------------------------------------------------------------------------
echo 'Reads filtered!';


#-------------------------------------------rename output files for simplicity
cp flexbarOut_barcode_plate_adapter_forward-plate_adapter_reverse_1.fastq RCP_R1.fastq;
cp flexbarOut_barcode_plate_adapter_forward-plate_adapter_reverse_2.fastq RCP_R2.fastq;
rm flexbarOut_barcode_plate_adapter_forward-plate_adapter_reverse_1.fastq;
rm flexbarOut_barcode_plate_adapter_forward-plate_adapter_reverse_2.fastq;


#--------------Make fasta files in chunks 
perl scripts/fastq2fasta.pl . ./RCP_R1.fastq ./RCP_R2.fastq;
ECHO "FASTQ split to FASTA complete!";
#-------------------------------------

#-----MAKE A BLAST DATABASE FOR VECTOR CONSTANT SEQUENCE
perl scripts/vector_database/make_const-seq_fasta.pl scripts/vector_database/newconstseq.txt >  scripts/vector_database/newconstseq.fa;
makeblastdb -in ~/Desktop/output/scripts/vector_database/newconstseq.fa -parse_seqids -dbtype nucl;
#-----------------------------------------------------------



#--------------make shell scripts to BLAST FASTA against reference
perl scripts/primers_blast_wrapper.pl . ./fragmented_fasta/*.fna;
#-------------------------------------


echo "BLAST'ing reads, please be patient :)";
#------------Run the blast shell scripts
chmod +x ~/Desktop/output/blast/sh.primers_blast/*.sh;


for f in ~/Desktop/output/blast/sh.primers_blast/*.sh; do  # or wget-*.sh instead of *.sh
  bash "$f" -H   || break # if needed 
done
echo 'Done!';
#-----------------------------------------------




#now infer barcode position using BLAST coordinates
#this script makes the shell script files 
perl scripts/barcode-identification_wrapper.pl . fragmented_fasta/*.fna;
#had to modify the barcode identification perl script, so also modifying the wrapper to point ot the new script V5





#this makes many shell scripts in the sh.identification directory that need to be run 
chmod +x QC_and_BC/sh.identification/*.sh
#---------this will run teh identify barcode V4 code, that needs to be made directory agnostic before proceesing 
cp $tags ./tag_assignment.csv

echo 'Identifying barcodes in reads, please be patient :)';
for f in QC_and_BC/sh.identification/*.sh; do  # or wget-*.sh instead of *.sh
  bash "$f" -H   || break # if needed 
done
echo "Done!";
#note that these files should be ~15.6 mb, if they arent check the tag_assignment file its probably got bad line breaks.
#open with bbedit, re-save with unix line breaks




#-----------After running the shell scripts for identification we now want to count barcode calls


perl ./scripts/read-counting_wrapper.pl . ./QC_and_BC/out.identification/*.dmp;


#------------ Now that we have made many .sh scripts we want to run them like the others 
chmod +x ~/Desktop/output/QC_and_BC/sh.count/*.sh;
echo 'Counting barcode calls in reads, please be patient :)';
for f in ~/Desktop/output/QC_and_BC/sh.count/*.sh; do  # or wget-*.sh instead of *.sh
  bash "$f" -H   || break # if needed 
done



#----------Now that we have counted barcodes for all of our clones we wnt to call the proper barcode for each based on counts 

#----first we merge count data 
mkdir Data;
perl ~/Desktop/output/scripts/merge_count-data.pl ~/Desktop/output/QC_and_BC/out.count/*.dmp > ~/Desktop/output/Data/merged_counts.dmp;

#--------Now that we have merged counts we can call barcodes 
echo "Calling barcodes!";
perl ~/Desktop/output/scripts/call_barcode.pl ~/Desktop/output/Data/merged_counts.dmp > ~/Desktop/output/Data/barcode_calls.dmp;
echo "Done!";

#----------now score wells 
echo 'Scoring calls!';
perl ~/Desktop/output/scripts/score_well.pl ~/Desktop/output/Data/barcode_calls.dmp > ~/Desktop/output/Data/well_scores.dmp;
echo 'Done!';



#-----------and lastly write barcode calls to spreadhseet
echo 'Writing Output as barcodes_calls_dataset.tsv in Data directory';
perl ~/Desktop/output/scripts/print_spreadsheet.well_info.pl ~/Desktop/output/Data/well_scores.dmp > ~/Desktop/output/Data/barcode_calls_dataset.tsv;
echo 'Analysis Complete';




#nopte that I am modifying all of the scripts that need to be modified for this shell script and I will store them in output scripts directory
#I will update github at teh end with the V2 RCP scripts 



