#!/bin/env/perl -w

use strict;
use warnings;

use POSIX;
# use List::Util;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqIO::fastq;
use Bio::Seq::Quality;
use List::Util qw(min max);
use Getopt::Long;
# use Pod::Usage;
use Carp qw(carp cluck);

my $VERSION="v0.0"; # First public release on GitHub

#### This script is run once for each barcode, from the qsub script job-phiXXer.v0.0.sh

sub usage {
	print("Usage: $0 -barcodeNum=VALUE -raw_reads_folder=VALUE -output_folder=VALUE -in_pCC1_trailer_leader=VALUE -CANU_MIN_COVERAGE=VALUE -my_RUN_metadata_tsv=VALUE -my_DEBUG_mode=VALUE\n");
	print "$0 requires 7 command-line arguments listed above\n";
}

print "$0 : $VERSION\n";
my $numArgs = $#ARGV + 1;

## Declare the Command-line arguments required as input-parameters:
my $barcodeNum;
my $in_folder_raw_reads;
my $out_folder_ALL;
my $in_pCC1_trailer_leader;
my $ASSEMBLY_MIN_COVERAGE;
my $my_RUN_metadata_tsv;
my $my_DEBUG_mode;

GetOptions("barcodeNum=i" => \$barcodeNum,
			"raw_reads_folder=s" => \$in_folder_raw_reads,
            "output_folder=s"   => \$out_folder_ALL,
            "in_pCC1_trailer_leader=s"   => \$in_pCC1_trailer_leader,
            "CANU_MIN_COVERAGE=i"   => \$ASSEMBLY_MIN_COVERAGE,
            "my_RUN_metadata_tsv=s"  => \$my_RUN_metadata_tsv,
            "my_DEBUG_mode=s"  => \$my_DEBUG_mode)
or die usage;

die usage unless $numArgs == 7;
if (!(-d $in_folder_raw_reads)) { die "Could not find in_folder_raw_reads = $in_folder_raw_reads\n";	}
if (!(-d $out_folder_ALL)) { die "Could not find out_folder_ALL = $out_folder_ALL\n";	}
if (!(-e $in_pCC1_trailer_leader)) { die "Could not find file in_pCC1_trailer = $in_pCC1_trailer_leader\n";	}
if (!(-e $my_RUN_metadata_tsv)) { die "Could not find file my_RUN_metadata_tsv = $my_RUN_metadata_tsv\n";	}

if (!($my_DEBUG_mode eq "ON")) {	$my_DEBUG_mode = "OFF"; } # Default should be "OFF;
print "my_DEBUG_mode = $my_DEBUG_mode\n";

my $barcode_padded_num = get_padded_barcode($barcodeNum);
my $barcode_padded_string = "barcode".$barcode_padded_num;
# print "DEBUG barcode_padded_string = $barcode_padded_string\n";

print "Opening the my_RUN_metadata_tsv file  : $my_RUN_metadata_tsv\n";
open (METADATA, "<", $my_RUN_metadata_tsv) or die $!;
my $metaNum = 0;
my $my_RUN_NAME;
my %hash_barcode_to_sampleWellName;
while (my $line = <METADATA>) {
	#ONT_barcode	clone_ID
	#1	150-C18
	chomp($line);
	if (!($line eq "")) { # Skip empty lines
		$line =~ s/[\x0A\x0D]//g; # To handle the WINDOWS-CRLF end-of-line EOL character that might sneak in through XL or a windows user
		$metaNum++;
		if (!($line =~ /^#/)) { # Ignore comment lines
			my ($ont_barcode, $clone_id) = split("\t", $line);
			chomp($ont_barcode);
			chomp($clone_id);
			if ($ont_barcode =~ m/^barcode/) {
				$ont_barcode =~ s/barcode//g;
			}
			# Chomp leading zeroes
			$ont_barcode =~ s/^0//g;
# 			print "DEBUG: Chomped leading zeroes for ont_barcode = $ont_barcode\n";
			$hash_barcode_to_sampleWellName{$ont_barcode} = $clone_id;
		}	
	}
}
close(METADATA);

my $clone_ID = "";
if (exists $hash_barcode_to_sampleWellName{$barcodeNum} ) { $clone_ID = $hash_barcode_to_sampleWellName{$barcodeNum} ;	}
else { die "Problems with the metadata tsv file: Could not find clone_ID for barcodeNum=$barcodeNum\n"; }

my $logfile = $clone_ID.".logfile.out";
open (LOGFILE, ">", $logfile) or die $!;

my $logdumps = $clone_ID.".logfile.dump"; # store temporary stderr /stdouts of various bash commands

my $thisBarcode_raw_reads_folder=$in_folder_raw_reads."/".$barcode_padded_string;
print LOGFILE "Using thisBarcode_raw_reads_folder = $thisBarcode_raw_reads_folder\n";
if (!(-d $thisBarcode_raw_reads_folder)) { 
	my $error_msg = "Could not find thisBarcode_raw_reads_folder = $thisBarcode_raw_reads_folder\n";
	print LOGFILE $error_msg;
	die $error_msg;	
}

my $datetime = `date`;
print"###############################################################\n";
print "Starting phiXXer for barcodeNum = $barcodeNum; clone_ID = $clone_ID\n\n"; 
print "Start datetime = $datetime\n";
print "Received: in_folder_raw_reads = $in_folder_raw_reads\n";
print "Received: out_folder_ALL = $out_folder_ALL\n";
print"####### Step 00: Creating all the required folders\n";

########################################################\n";
print LOGFILE "###############################################################\n";
print LOGFILE "Starting phiXXer for barcodeNum = $barcodeNum; clone_ID = $clone_ID\n\n"; 
print LOGFILE "Start datetime = $datetime\n";
print LOGFILE "Received: in_folder_raw_reads = $in_folder_raw_reads\n";
print LOGFILE "Received: out_folder_ALL = $out_folder_ALL\n";
print LOGFILE "####### Step 00: Creating all the required folders\n";

########################################################\n";

my $bash_cmd;

my $tmp_barcodes_processing_dir = $out_folder_ALL."/"."tmp.barcodes_processing_dir"; # Collect all 96 barcodes intermediate files in one dir
create_folder($tmp_barcodes_processing_dir);

my $thisBarcode_outputs_folder=$tmp_barcodes_processing_dir."/".$clone_ID;
if (!(-d $thisBarcode_outputs_folder)) {
	$bash_cmd="mkdir $thisBarcode_outputs_folder";
	execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "continue");
}

# Output folder for each step will be created
## Barcode specific directories
my $thisBarcode_outdir01 = $thisBarcode_outputs_folder."/"."out01.combined_raw_reads";
my $thisBarcode_outdir02 = $thisBarcode_outputs_folder."/"."out02.shredded_fastqs"; 
my $thisBarcode_outdir03 = $thisBarcode_outputs_folder."/"."out03.assembled_fosmid"; 
my $thisBarcode_outdir04 = $thisBarcode_outputs_folder."/"."out04.tmp_files"; 
create_folder($thisBarcode_outdir01);
create_folder($thisBarcode_outdir02);
create_folder($thisBarcode_outdir03);
create_folder($thisBarcode_outdir04);

## collect final outputs from all barcodes in one super-folder
my $dir_canu_assemblies_raw = $out_folder_ALL."/"."00.canu_assemblies"; 
my $dir_canu_assemblies_higestCov = $out_folder_ALL."/"."01.canu_contigs_highestCov"; 
my $dir_final_insert_pcc1_trimmed = $out_folder_ALL."/"."02.final_insert.pcc1_trimmed"; 
my $dir_final_logs = $out_folder_ALL."/"."03.final_logs"; 

create_folder($dir_canu_assemblies_raw);
create_folder($dir_canu_assemblies_higestCov);
create_folder($dir_final_insert_pcc1_trimmed);
create_folder($dir_final_logs);

## Collect internediate files for debugging
my $outdir03_tmp_logs_etc = $out_folder_ALL."/"."tmp.logs_etc"; 
create_folder($outdir03_tmp_logs_etc);

print LOGFILE "####### Step 00: DONE.\n\n";
print "####### Step 00: DONE.\n\n";


#########################################################################################
### STEP 01: Combine raw reads into one fastq file, discard junk reads shorter than 150 bp, assign"goodnames"
#########################################################################################
print LOGFILE "####### Step 01: Combine raw reads into one fastq file\n";
print "####### Step 01: Combine raw reads into one fastq file\n";
my $thisBarcode_raw_fastq = $thisBarcode_outdir01."/".$clone_ID.".raw_reads_all.fastq"; # This file will be created
my $thisBarcode_raw_fastq2tab = $thisBarcode_outdir01."/".$clone_ID.".raw_reads_all.fastq2tab.tsv"; # This file will be created

## Loop through the inout raw reads folder and combine naopore output fastqs into one
if ( -e $thisBarcode_raw_fastq) {
	# leftover from previous run,should be deleted
	print LOGFILE "Deleting the leftovers from previous run, file = $thisBarcode_raw_fastq\n";
	$bash_cmd = "rm -f $thisBarcode_raw_fastq";
	execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "die");
}
$bash_cmd = "touch $thisBarcode_raw_fastq"; # initialize file
execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "die");

opendir(DIR_HANDLE, "$thisBarcode_raw_reads_folder") or die "Could not open $thisBarcode_raw_reads_folder\n\n";
my @raw_files = readdir(DIR_HANDLE);
closedir(DIR_HANDLE);
foreach my $file (@raw_files) {
	# skip . and ..
	next if($file =~ /^\.$/);
	next if($file =~ /^\.\.$/);
	if ($file =~/.gz$/) {
		my $filepath = $thisBarcode_raw_reads_folder."/".$file;
		$bash_cmd = "zcat $filepath >> $thisBarcode_raw_fastq";
		execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "NO", do_on_fail => "die");
	}
} 

### Assign "goodNames" to readNames. Use fx2tab for easy parsing of fastq records
my $tmp_fq_goodnames = $thisBarcode_raw_fastq."tmp.fastq"; # Create a temporary fq
open (TMP_GOODNAMES_FQ, ">", $tmp_fq_goodnames) or die $!;

my $tmp_fq_goodnames_log = $thisBarcode_raw_fastq.".goodnames-log.tsv"; # Log the read-name changes in case needed later
open (TMP_GOODNAMES_LOG, ">", $tmp_fq_goodnames_log) or die $!;

my $num_fqReads_input = 0;
my $num_fqReads_output = 0;

$bash_cmd = "seqkit fx2tab $thisBarcode_raw_fastq  > $thisBarcode_raw_fastq2tab";
execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "die");

open (FQ2TAB, "<", $thisBarcode_raw_fastq2tab) or die $!;
while (my $fastq_read = <FQ2TAB>) {
	$num_fqReads_input++;
	chomp($fastq_read);
	my ($fq_id, $fq_seq, $fq_qual) = split("\t", $fastq_read);
	my $fq_read_len = length($fq_seq);
	if ($fq_read_len >= 150) { # Get rid of nanopore reads which are shorter than Illumina reads, sometimes even just 1bp!!??
		my $thisReadName = $clone_ID."_readNum".$num_fqReads_input."_len".$fq_read_len; # e.g. barcode01_readNum1_len13456
		my $thisRead_fasta = $thisReadName.".fasta";
		print TMP_GOODNAMES_FQ "@", $thisReadName, "\n", $fq_seq, "\n", "+", "\n", $fq_qual,"\n";
		print TMP_GOODNAMES_LOG $thisReadName, "\t", $fq_id, "\n"; # log the name change
		$num_fqReads_output++;
	}	
}
close(FQ2TAB);
close(TMP_GOODNAMES_FQ);
close(TMP_GOODNAMES_LOG);
print LOGFILE "QC_INFO:clone_ID = $clone_ID;  num_fqReads_input = $num_fqReads_input;	After discarding reads shorter than 150 bp, num_fqReads_output = $num_fqReads_output\n";
### Finish the "goodnames" business and make fx2tab for next steps
	$bash_cmd = "mv $tmp_fq_goodnames $thisBarcode_raw_fastq";
	execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "die");

	$bash_cmd = "seqkit fx2tab $thisBarcode_raw_fastq  > $thisBarcode_raw_fastq2tab";
	execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "die");

print LOGFILE "####### Step 01: DONE.\n\n";
print "####### Step 01: DONE.\n\n";

#########################################################################################
### STEP 02: Shred each read in the fastq file by removing the pCC1 matches
#########################################################################################
print LOGFILE "####### Step 02: Shred each read and remove pCC1 using nucmer matches.\n";
print "####### Step 02: Shred each read and remove pCC1 using nucmer matches.\n";
# my $thisBarcode_shredded_fastq =$thisBarcode_outdir02."/".$clone_ID.".shredded_all.fastq"; # This file will be created
my $thisBarcode_shredded_fastA=$thisBarcode_outdir02."/".$clone_ID.".shredded_all.fasta"; # This file will be created
# open (SHREDDED_FASTQ, ">", $thisBarcode_shredded_fastq) or die $!;
open (SHREDDED_FASTA, ">", $thisBarcode_shredded_fastA) or die $!;
open (FQ2TAB, "<", $thisBarcode_raw_fastq2tab) or die $!;

$num_fqReads_input = 0; # Reset counts, as the inout has been reset to the "goodnames" fastX
my $num_fqReads_output_shredded = 0;
my $num_fq_with_pCC1_hit = 0;
my $num_fq_with_pCC1_strand_flip_ERROR = 0;

while (my $fastq_read = <FQ2TAB>) {
	$num_fqReads_input++;
	chomp($fastq_read);
	my ($fq_id, $fq_seq, $fq_qual) = split("\t", $fastq_read);
	my $fq_read_len = length($fq_seq);
	if ($fq_read_len >= 150) { # Get rid of nanopore reads which are shorter than Illumina reads, sometimes even just 1bp!!??
		my $thisReadName = $clone_ID."_readNum".$num_fqReads_input."_len".$fq_read_len; # e.g. barcode01_readNum1_len13456
		my $thisRead_fasta = $thisReadName.".fasta";
		open (READS_FASTA, ">", $thisRead_fasta) or die $!;
		print READS_FASTA ">", $thisReadName, "\n", $fq_seq, "\n"; 
		close(READS_FASTA);
		
		# Run nucmer against pCC1
		my $nucmer_prefix = $thisReadName."_vs_pcc1";
		my $nucmer2bed_file = $nucmer_prefix.".qclTH.coords2bed.bed"; # will be created by run_nucmer_to_bed_steps

		my $delta_file = $nucmer_prefix.".delta"; # was created by run_nucmer_to_bed_steps
		my $coords_file = $nucmer_prefix.".qclTH.coords"; # was be created by run_nucmer_to_bed_steps

		my $num_delta_lines = run_nucmer_to_bed_steps(nucmer_ref => $thisRead_fasta, 
			nucmer_query => $in_pCC1_trailer_leader, 
			nucmer_prefix => $nucmer_prefix,  
			nucmer2bed_file => $nucmer2bed_file); 

		if ($num_delta_lines == 2) {
			# no PCC1 match, dump the read as it is
			print SHREDDED_FASTA ">", $thisReadName, "\n", $fq_seq, "\n"; #, "+", "\n", $fq_qual,"\n";	
		}
		else { # Alignments were found, start shredding
			my $strand = getStrand_from_coords_file($coords_file);
# 			print LOGFILE "DEBUG: nucmer: strand = $strand in coords_file = $coords_file\n";
			if ($strand eq "ERROR_FLIPPED_STRAND") {
				# Document this chimeric read
				print LOGFILE "WARNING: ERROR_FLIPPED_STRAND in a raw read: coords_file = $coords_file\n";
# 				print "WARNING: ERROR_FLIPPED_STRAND in a raw read: coords_file = $coords_file\n";
				## Dump the output to logfile
# 				my $coords_header = "[strand]". "\t". "[S1]". "\t"."[E1]". "\t"."[S2]". "\t"."[E2]". "\t"."[LEN 1]". "\t"."[LEN 2]". "\t"."[% IDY]". "\t"."[LEN R]". "\t"."[LEN Q]". "\t"."[COV R]". "\t"."[COV Q]". "\t"."[TAGS]";
				my $coords_header = "[S1]". "\t"."[E1]". "\t"."[S2]". "\t"."[E2]". "\t"."[LEN 1]". "\t"."[LEN 2]". "\t"."[% IDY]". "\t"."[LEN R]". "\t"."[LEN Q]". "\t"."[COV R]". "\t"."[COV Q]". "\t"."[TAGS]";
				print LOGFILE $coords_header, "\n";
				dump_to_LOGFILE($coords_file);
			}
			my $bed_frags_to_get_file = $thisReadName.".frags2get.bed";
			run_bed_substract($thisRead_fasta, $nucmer2bed_file, $bed_frags_to_get_file);
			dump_to_LOGFILE($bed_frags_to_get_file);
			my $fasta_shreds_as_string = get_fasta_fromfrags_bed_AS_A_STRING($bed_frags_to_get_file, $thisRead_fasta); # , $out_split_fasta);
			print SHREDDED_FASTA $fasta_shreds_as_string, "\n"; # , $fq_seq, "\n"; #, "+", "\n", $fq_qual,"\n";		

			$bash_cmd = "rm -f $bed_frags_to_get_file"; # $outdir03_tmp_logs_etc/";
			execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "NO", do_on_fail => "nothing");

		}
		##### Clean-up : DELETE the nucmer etc files
		$bash_cmd = "rm -f $delta_file"; # $outdir03_tmp_logs_etc/";
		execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "NO", do_on_fail => "nothing");
		$bash_cmd = "rm -f $coords_file"; # $outdir03_tmp_logs_etc/";
		execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "NO", do_on_fail => "nothing");
		$bash_cmd = "rm -f $thisRead_fasta"; # $outdir03_tmp_logs_etc/";
		execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "NO", do_on_fail => "nothing");
		$bash_cmd = "rm -f $nucmer2bed_file"; # $outdir03_tmp_logs_etc/";
		execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "NO", do_on_fail => "nothing");
	}		
}
close(FQ2TAB);

# close(TMP_GOODNAMES_FQ);
close(TMP_GOODNAMES_LOG);
print LOGFILE "QC_INFO:clone_ID = $clone_ID;  num_fqReads_input = $num_fqReads_input;	num_fqReads_output_shredded = $num_fqReads_output_shredded\n"; 


print LOGFILE "####### Step 02: DONE.\n\n";
print "####### Step 02: DONE.\n\n";


#########################################################################################
### STEP 03: Canu assembly of shredded reads
#########################################################################################
print LOGFILE "####### Step 03: Canu assembly of shredded reads.\n";
print "####### Step 03: Canu assembly of shredded reads.\n";
my $canu_prefix = $clone_ID."_canu";
print LOGFILE "DEBUG: Using canu_prefix = $canu_prefix\n";

my $thisBarcode_assembly = $thisBarcode_outdir03."/".$canu_prefix.".contigs.fasta"; # This file will be created, default canu output
my $canu_logs = $thisBarcode_outdir03."/".$canu_prefix.".logs";
# $bash_cmd = "canu useGrid=remote gridEngineMemoryOption=\"-l mem_free=MEMORY\" stopOnLowCoverage=$ASSEMBLY_MIN_COVERAGE minInputCoverage=$ASSEMBLY_MIN_COVERAGE -p $canu_prefix -d $thisBarcode_outdir03 genomeSize=40k -nanopore $thisBarcode_shredded_fastA &> $canu_logs"; 
$bash_cmd = "canu useGrid=true gridEngineMemoryOption=\"-l mem_free=MEMORY\" stopOnLowCoverage=$ASSEMBLY_MIN_COVERAGE minInputCoverage=$ASSEMBLY_MIN_COVERAGE -p $canu_prefix -d $thisBarcode_outdir03 genomeSize=40k -nanopore $thisBarcode_shredded_fastA &> $canu_logs"; 
execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "die");


## CHECK: DEBUG If assemblies are produced
# But first, Store the canu logs to investigate in case the assembly of this barcode fails
$bash_cmd = "cp $canu_logs $dir_canu_assemblies_raw/";
execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "die");

if (!( -e $thisBarcode_assembly)) {
	print LOGFILE "\nQC_INFO : CANU assembly failed for clone_ID = $clone_ID; Could not find file thisBarcode_assembly = $thisBarcode_assembly\n";
	print LOGFILE "\nERROR : CANU assembly failed for clone_ID = $clone_ID; Could not find file thisBarcode_assembly = $thisBarcode_assembly\n";
	exit; # Nothing to do if nothing gets assembled
}	
else { # ELSE FOR: if (!( -e $thisBarcode_assembly))
	# Process the output contigs
	print LOGFILE "\nQC_INFO : CANU assembly successful for clone_ID = $clone_ID\n";

	# Move out the canu assembly
	$bash_cmd = "cp $thisBarcode_assembly $dir_canu_assemblies_raw/";
	execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "die");
	
	## TODO???: add_cloneID_to_fastaID($out_canu_full_path, $clone_ID);

	$bash_cmd = "seqkit stats $thisBarcode_assembly &> $logdumps";
	execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "continue");
	dump_to_LOGFILE($logdumps);

	$bash_cmd = "grep '>' $thisBarcode_assembly &> $logdumps";
	execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "continue");
	dump_to_LOGFILE($logdumps);
	
}

print LOGFILE "####### Step 03: DONE.\n\n";
print "####### Step 03: DONE.\n\n";


#########################################################################################
### STEP 04: Select the highest coverage CANU contig
###		 - Discard low coverage contigs in each barcode
#########################################################################################
print LOGFILE "####### Step 04: Select the highest coverage CANU contig\n";
print "####### Step 04: Select the highest coverage CANU contig\n";

my $numContigs_in_barcode = getNumContigs_in_fasta($thisBarcode_assembly);
print LOGFILE "QC_INFO: For clone_ID = $clone_ID : Number of CANU contigs = $numContigs_in_barcode\n";

# 4a : Parse canu headers to select the contig with highest coverage
my $input_genomeObject = Bio::SeqIO->new(-file => "$thisBarcode_assembly", -format => 'Fasta') or die $!; 
my $final_contig_with_highestCov__name;
my $final_contig_with_highestCov__header;
my $final_contig_with_highestCov__seq;
my $highestCoverage = 0; # default initialization

my $out_canu_best_fasta = $dir_canu_assemblies_higestCov."/".$clone_ID.".fasta"; # Final output
open (CANU_BEST_FASTA_by_COV, ">", $out_canu_best_fasta) or die $!;
# 	print LOGFILE "\nDEBUG: Opened CANU_BEST_FASTA_by_COV out_canu_best_fasta = $out_canu_best_fasta\n";

my $numContigs = 0;
while ( my $contigObj = $input_genomeObject->next_seq) {	
	$numContigs++;
#	print "\nDEBUG> Inside while Loop : Processing numContigs = $numContigs.\n";
	my $contig_id = $contigObj->id();
	my $contig_description= $contigObj->desc();
	my $contig_seq = $contigObj->seq();
	my $contig_len = length($contig_seq);

	my $contig_name_to_use = $contig_id;
	my $contig_readCoverage_fromHeader = get_contig_readCoverage_fromHeader($contig_description);
# 	print "\nDEBUG> COVERAGE: contig_readCoverage_fromHeader = $contig_readCoverage_fromHeader; highestCoverage = $highestCoverage\n";
# 	if (!($contig_id =~ /^bc/)) {	$contig_name_to_use = $barcode_padded_string."_".$contig_id;	} # make sure the barcode prefix is captured
# 	if (!($contig_id =~ /^$clone_ID/)) {	$contig_name_to_use = $clone_ID."_".$contig_id;	} # make sure the barcode prefix is captured
	if ($contig_readCoverage_fromHeader > $highestCoverage) {
		$highestCoverage = $contig_readCoverage_fromHeader;
		$final_contig_with_highestCov__name = $contig_name_to_use;
		$final_contig_with_highestCov__name = $clone_ID; ## Keep the simplest name
		my $final_contig_info_dump = $contig_id." ".$contig_description; ## Capture all the old info and dump to new contig_description
		$final_contig_with_highestCov__header = ">".$clone_ID." ".$final_contig_info_dump;
		$final_contig_with_highestCov__seq = $contig_seq;
	}
}
print LOGFILE "DEBUG: Selected Best Canu contig: final_contig_with_highestCov__header = $final_contig_with_highestCov__header\n";
print CANU_BEST_FASTA_by_COV $final_contig_with_highestCov__header, "\n", $final_contig_with_highestCov__seq, "\n";
close(CANU_BEST_FASTA_by_COV);
# 	print LOGFILE "Printed $final_contig_with_highestCov__header to CANU_BEST_FASTA_by_COV = $out_canu_best_fasta\n";
print LOGFILE "####### Step 04: DONE.\n\n";
print "####### Step 04: DONE.\n\n";


#########################################################################################
### STEP 05: Split the assembled contig at-  and eliminate- the fosmid backbone = in_pCC1_trailer_leader
### 	Based on some code from : /Users/asinha/Data/neb-collab-projects/lea_fosmid_screening/perl.shred-and-re-assemble.pl.v3
#########################################################################################
print LOGFILE "####### Step 05: Trim out the pCC1 fosmid backbone and generate consensus insert for clone_ID = $clone_ID\n";
print "####### Step 05: Trim out the pCC1 fosmid backbone and generate consensus insert for clone_ID = $clone_ID\n";

my $nucmer_query_fasta = $out_canu_best_fasta; # output from previous step

my $nucmer_prefix = $clone_ID.".trim_pcc1";

# Declare the potential output files for this step	
# my $out_consensus_fasta_noPCC =  $clone_ID.".out_consensus.fasta"; 
# my $out_consensus_fasta_noPCC =  $dir_final_insert_pcc1_trimmed."/".$clone_ID.".out_consensus.fasta"; # Final output file e.g. NEB-260-M1.out_consensus.fasta
my $out_consensus_fasta_noPCC =  $dir_final_insert_pcc1_trimmed."/".$clone_ID.".fasta"; # Final output file: e.g. NEB-260-M1.fasta
my $delta_file = $nucmer_prefix.".delta";
my $coords_file = $nucmer_prefix.".qclTH.coords";
my $coords_header = "[strand]". "\t". "[S1]". "\t"."[E1]". "\t"."[S2]". "\t"."[E2]". "\t"."[LEN 1]". "\t"."[LEN 2]". "\t"."[% IDY]". "\t"."[LEN R]". "\t"."[LEN Q]". "\t"."[COV R]". "\t"."[COV Q]". "\t"."[TAGS]";
my $nucmer2bed_file = $coords_file."2bed.bed";
my $bed_frags_to_get_file = $nucmer_prefix.".frags_to_get.bed";
my $out_split_fasta = $nucmer_prefix.".split-for-consensus.fasta"; # 

# my $out_emma_msa = $nucmer_prefix.".aln"; #
# my $out_emma_dnd = $nucmer_prefix.".dnd"; #
# my $out_emma_stderrs = $nucmer_prefix.".stderr"; # capture emma and cons stderr which creates verbiage in log files

# Run nucmer
$bash_cmd = "nucmer -p $nucmer_prefix $in_pCC1_trailer_leader $nucmer_query_fasta";
execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "die");


my $num_delta_lines = getNumLines($delta_file);
if ($num_delta_lines == 2)	{	
	# No matches were found to the pCC1 fosmid trailer sequence, e.g. a contaminant from host E. coli. 
	# But the USER might still want this sequence
	# Dump it out with a warning
	print LOGFILE "QC_INFO: nucmer: No nucmer alignments between pcc1_fosmid and the assembled contig in $nucmer_query_fasta. Proceed with caution\n";		
	print LOGFILE "WARNING: nucmer: No nucmer alignments between pcc1_fosmid and the assembled contig in $nucmer_query_fasta. Proceed with caution\n";		

	$bash_cmd = "cp $nucmer_query_fasta $out_consensus_fasta_noPCC";
	execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "die");

	my $fasta_id_toprint = $clone_ID;
	my $fasta_desc_toprint = "readCoverageCanu=$highestCoverage;";
	$fasta_desc_toprint = $fasta_desc_toprint." INFO: No alignments found to pCC1 fosmid backbone.";
	rename_fasta_single($out_consensus_fasta_noPCC, $fasta_id_toprint, $fasta_desc_toprint);
}
elsif ($num_delta_lines > 2 ) {
	# ($num_delta_lines > 2) only when alignments are found. The first two lines are program specific headers
	#	$bash_cmd = "show-coords -rclTH $delta_file > $coords_file"; # -r => Sorted by reference = pCC1_trailer
	$bash_cmd = "show-coords -qclTH $delta_file > $coords_file"; # -q => Sorted by reference = pCC1_trailer
	execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "die");

	my $strand = getStrand_from_coords_file($coords_file);
	print LOGFILE "DEBUG: nucmer: strand = $strand in coords_file = $coords_file\n";
	if ($strand eq "ERROR_FLIPPED_STRAND") {
		print LOGFILE "ERROR: ERROR_FLIPPED_STRAND in coords_file = $coords_file\n";
		print "ERROR: ERROR_FLIPPED_STRAND in coords_file = $coords_file\n";
		## Dump the output to logfile
			print LOGFILE $coords_header, "\n";
			dump_to_LOGFILE($coords_file);

		$bash_cmd = "cp $nucmer_query_fasta $out_consensus_fasta_noPCC";
		execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "die");

		my $fasta_id_toprint = $clone_ID;
		my $fasta_desc_toprint = "readCoverageCanu=$highestCoverage;";
		$fasta_desc_toprint = $fasta_desc_toprint." WARNING: Possible chimera - proceed with caution. pCC1 fosmid backbone aligned in both forward and reverse directions.";
		rename_fasta_single($out_consensus_fasta_noPCC, $fasta_id_toprint, $fasta_desc_toprint);			
		exit; ## TODO MOPRE STUFF BEFORE EXITing??		
	}
	elsif ($strand eq "minus") {
		print LOGFILE "DEBUG: nucmer: Re-doing nucmer + show-coords on the reverse-complemented contig\n";
		# 0. reverse-complement the query contig
		write_ReverseComplemented_Fasta($nucmer_query_fasta);		
		# 1. redo the nucmer and show-coords			
		$bash_cmd = "nucmer -p $nucmer_prefix $in_pCC1_trailer_leader $nucmer_query_fasta";
		execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "die");
	
		$bash_cmd = "show-coords -qclTH $delta_file > $coords_file"; # -q => Sorted by reference = pCC1_trailer
		execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "die");
	}
	## Dump the coords_file output to logfile
		print LOGFILE $coords_header, "\n";
		dump_to_LOGFILE($coords_file);

	convert_coords_to_bed($coords_file, $nucmer2bed_file, "QUERY");	
	run_bed_substract($nucmer_query_fasta, $nucmer2bed_file, $bed_frags_to_get_file);
	dump_to_LOGFILE($bed_frags_to_get_file);

	## Get fasta frags from bed
	get_fasta_fromfrags_bed($bed_frags_to_get_file, $nucmer_query_fasta, $out_split_fasta);
	my $num_contigs_post_pcc1_trim = getNumContigs_in_fasta($out_split_fasta);
	print LOGFILE "QC_INFO: nucmer: num_contigs_post_pcc1_trim = $num_contigs_post_pcc1_trim; $nucmer_query_fasta\n";
	if ($num_contigs_post_pcc1_trim == 1) {
		# The pcc1 fosmid backbone is most likley at the terminus of assembled contig and was successfully removed.
		print LOGFILE "QC_INFO: nucmer: num_contigs_post_pcc1_trim = $num_contigs_post_pcc1_trim. The pCC1 fosmid backbome was most likely at the terminii\n";		

		$bash_cmd = "cp $out_split_fasta $out_consensus_fasta_noPCC";
		execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "die");
	
		my $fasta_id_toprint = $clone_ID;
		my $fasta_desc_toprint = "readCoverageCanu=$highestCoverage;";
		$fasta_desc_toprint = $fasta_desc_toprint." NOTE: num_contigs_post_pcc1_trim == 1; pCC1 backbone was most probably at the terminii.";
		rename_fasta_single($out_consensus_fasta_noPCC, $fasta_id_toprint, $fasta_desc_toprint);
	}
	else { ## ($num_contigs_post_pcc1_trim is > 1) 
		## Get merged:consensus using EMBOS:merger utiliy
		my $numContigs_after_merger = emboss_merger(input => $out_split_fasta, output => $out_consensus_fasta_noPCC);
# 		print "DEBUG call to emboss_merger: ARGS => out_split_fasta = $out_split_fasta, out_consensus_fasta_noPCC = $out_consensus_fasta_noPCC\n\n";
# 		emboss_merger(input => $out_split_fasta, output => $out_consensus_fasta_noPCC);
# 		my $numContigs_after_merger = getNumContigs_in_fasta($out_consensus_fasta_noPCC);
		print LOGFILE "numContigs_after_merger = $numContigs_after_merger\n";
			$bash_cmd = "seqkit stats $out_consensus_fasta_noPCC &> $logdumps";
			execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "continue");
			dump_to_LOGFILE($logdumps);
		# Check if merged:consensus is collapsed to one contig. raise WARNING otherwise
		if ($numContigs_after_merger == 1) {
			print LOGFILE "QC_INFO: nucmer: Step 05 : PASS: After pcc1_fosmid_trimming and merger, numContigs_after_merger = $numContigs_after_merger\n";	
			my $fasta_id_toprint = $clone_ID;
			my $fasta_desc_toprint = "readCoverageCanu=$highestCoverage;";
			$fasta_desc_toprint = $fasta_desc_toprint." NOTE: trim_pcc1_fosmid = Success;";
			rename_fasta_single($out_consensus_fasta_noPCC, $fasta_id_toprint, $fasta_desc_toprint);
		}
		else {
			# Leave the split contigs in the output file, the user might still be able to extract useful info
			print LOGFILE "QC_INFO: nucmer: Step 05 : FAIL: After pcc1_fosmid_trimming and merger, numContigs_after_merger = $numContigs_after_merger. (Should have been 1!). STOPPING!!!\n";	
			print LOGFILE "WARNING: nucmer: Step 05 : FAIL: After pcc1_fosmid_trimming and merger, numContigs_after_merger = $numContigs_after_merger. (Should have been 1!). STOPPING!!!\n";
			my $out_consensus_fasta_FAILED_CONSENSUS = $out_consensus_fasta_noPCC."failed_consensus_insert.fasta";
			$bash_cmd = "mv $out_consensus_fasta_noPCC $out_consensus_fasta_FAILED_CONSENSUS";
			execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "die");
			exit;
		}
	} # 		else { ## ($num_contigs_post_pcc1_trim is > 1) 
}	
# Cleanup nucmer files
if ( -e $delta_file) {	`mv $delta_file $outdir03_tmp_logs_etc/`;	}
if ( -e $coords_file) {	`mv $coords_file $outdir03_tmp_logs_etc/`;	}
if ( -e $nucmer2bed_file) {	`mv $nucmer2bed_file $outdir03_tmp_logs_etc/`;	}
if ( -e $bed_frags_to_get_file) {	`mv $bed_frags_to_get_file $outdir03_tmp_logs_etc/`;	}
if ( -e $out_split_fasta) {	`mv $out_split_fasta $outdir03_tmp_logs_etc/`;	}
# if ( -e $out_emma_msa) {	`mv $out_emma_msa $outdir03_tmp_logs_etc/`;	}
# if ( -e $out_emma_dnd) {	`mv $out_emma_dnd $outdir03_tmp_logs_etc/`;	}
# if ( -e $out_emma_stderrs) {	`mv $out_emma_stderrs $outdir03_tmp_logs_etc/`;	}
	# Dont move around the out_consensus_fasta_noPCC file yet, will be needed in the next trimming step
		
print LOGFILE "####### Step 05: DONE.\n\n";
print "####### Step 05: DONE.\n\n";

### Cleanup all intermediate files
# if (!($my_DEBUG_mode eq "ON")) {
	print LOGFILE "####### Step xx: Cleaning Up : \n";
	print "####### Step xx: Cleaning Up : \n";
	$bash_cmd = "rm -fr $thisBarcode_outputs_folder";
	execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "continue");
	
	$bash_cmd = "rm -fr $outdir03_tmp_logs_etc";
	execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "continue");
# }

# Delete the logdumps
$bash_cmd = "rm -f $logdumps";
execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "continue");

$datetime = `date`;
print LOGFILE "END processing this_barcode $clone_ID: datetime = $datetime\n";
print "END processing this_barcode $clone_ID: datetime = $datetime\n";

if (-e $logfile) {
	$bash_cmd = "mv $logfile $dir_final_logs/";
	print LOGFILE "bash_cmd = $bash_cmd\n";
	print LOGFILE "####### Step xx: DONE.\n\n";
	print "####### Step xx: DONE.\n\n";
	close(LOGFILE);
	system($bash_cmd) && die "DIED at bash_cmd = $bash_cmd;\n$!"; # A successful execution of bash_cmd returns ZERO
}

print "####### PHIXXER : DONE.\n\n";
print LOGFILE "####### PHIXXER : DONE.\n\n";



########################################################## END MAIN ##########################################################

########################################################## SUBS BELOW ##########################################################

##########################################################
sub get_padded_barcode {
	# invoked as: get_padded_barcode($barcodeNum);
	my $num = $_[0];
	my $to_return = $num;
	if ($num < 10) {
		$to_return = "0".$num;
	}
	return $to_return;
}

##########################################################
sub execute_and_log_bash_cmd {
    my %params = @_;
    my $bash_cmd = $params{bash_cmd};
    my $print_to_log = $params{print_to_log};
    my $do_on_fail = $params{do_on_fail};

	print LOGFILE "bash_cmd = $bash_cmd\n" if ($print_to_log eq "YES");

	my $bash_cmd_returnVal = system($bash_cmd);

	if ($bash_cmd_returnVal == 0) {		# A successful execution of bash_cmd returns ZERO
		print LOGFILE "bash_cmd = SUCCESS.\n\n" if ($print_to_log eq "YES");
	}	
	else {
		if ($do_on_fail eq "continue") {
			print LOGFILE "bash_cmd = FAILED. Continuing anyway...\n\n" if ($print_to_log eq "YES");	
		}	
		elsif ($do_on_fail eq "die") { 
			print LOGFILE "bash_cmd = FAILED. DIE-ing.\n" if ($print_to_log eq "YES");
			print "bash_cmd = FAILED. DIE-ing at bash_cmd = $bash_cmd.\n\n" if ($print_to_log eq "YES");
			die $!;
		}	
		elsif ($do_on_fail eq "carp") { 
			print LOGFILE "bash_cmd = FAILED. CARP-ing.\n" if ($print_to_log eq "YES");
			print "bash_cmd = FAILED. CARP-ing at bash_cmd = $bash_cmd.\n\n" if ($print_to_log eq "YES");
			carp $!;
		}	
	}
}   
    
##########################################################
sub create_folder {
	# invoked as: create_folder($thisBarcode_outdir04);
	my $folder_to_create = $_[0];
	if (!(-d $folder_to_create)) {
		$bash_cmd="mkdir $folder_to_create";
#  		print LOGFILE "bash_cmd: $bash_cmd\n";
# 		system($bash_cmd); # && carp "FAILED at bash_cmd = $bash_cmd;\n$!"; # A successful execution of bash_cmd returns ZERO
# 		print LOGFILE "bash_cmd = SUCCESS.\n\n";
		execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "NO", do_on_fail => "continue");
	}	
}

##########################################################
sub getNumLines {
	# get number of lines in a file : Substitute for bash "wc -l"
	my $in_file = $_[0];
	my $num_lines = 0;
	chomp($in_file);
	open (TMP, "<", $in_file) or die "Could not open file $in_file: $!";
	while (my $line = <TMP>) {
		$num_lines++;
	}
	close(TMP);
# 	my $num_lines = `wc -l $in_file`;
#	$num_lines =~ s/  / /g;
# 	print "in_func=$in_file\n";

	return $num_lines;
}

##########################################################
sub getStrand_from_coords_file {
	# Invokded as:
	# my $strand = getStrand_from_coords_file($coords_file);
	my $coords_file = $_[0];
	my $strand = "";
	my $prev_strand ="";
	my $line_num = 0;
	open (TMP_IN, "<", $coords_file) or die $!;
	while (my $line = <TMP_IN>) {
		chomp($line);
		my @x = split("\t", $line);
		#[S1]	[E1]	[S2]	[E2]	[LEN 1]	[LEN 2]	[%	IDY]	[LEN R]	[LEN Q]	[COV R]	[COV Q]	[TAGS]";
		#	1	365	19736	20100	365	365	100.00	365	197300	100.00	0.18	my_pCC1_leader	bc27_tig00000001
		#	1	365	49544	49908	365	365	100.00	365	197300	100.00	0.18	my_pCC1_leader	bc27_tig00000001
		$line_num++;

		my $start = $x[2] - 1;
		my $end = $x[3] - 1;	
		if ($start < $end) {
			if ($line_num == 1) {
				$prev_strand = "plus";
				$strand = "plus";
			}
			else {
				$prev_strand = $strand;
				$strand = "plus";
			}
		
		}
		elsif ($start > $end) {
			if ($line_num == 1) {
				$prev_strand = "minus";
				$strand = "minus";
			}
			else {
				$prev_strand = $strand;
				$strand = "minus";
			}
		}
# 		my $coords_header = "[strand]". "\t". "[S1]". "\t"."[E1]". "\t"."[S2]". "\t"."[E2]". "\t"."[LEN 1]". "\t"."[LEN 2]". "\t"."[% IDY]". "\t"."[LEN R]". "\t"."[LEN Q]". "\t"."[COV R]". "\t"."[COV Q]". "\t"."[TAGS]";
# 		print LOGFILE "\n", $coords_header, "\n", $strand, "\t", $line, "\n";
		if (!($strand eq $prev_strand)) {
			# If the phi29 amplicons are not all going in the same dirction e.g. multiple primers were used, 
			# weird artefacts are produced. DIE if strand-switching is encountered.
# 			my $ERROR_MSG="\nERROR: Inside sub:getStrand_from_coords_file: Strand changed direction from prev_strand=$prev_strand to strand=$strand\: in coords_file=$coords_file; line_num = $line_num;\ncoords_line = $line\n";
# 			print LOGFILE $ERROR_MSG; # goes to job log
# 			print $ERROR_MSG; # goes to job log
			$strand = "ERROR_FLIPPED_STRAND";
			return $strand;
			exit;
		}
	}
	print LOGFILE "\n";
	return $strand;
}	

############################################################
sub run_nucmer_to_bed_steps {
# 	invoked as: run_nucmer_to_bed_steps($nucmer_ref $nucmer_query  $nucmer_prefix $nucmer2bed_file);
	my %params = @_;
	my %valid = map { $_ => 1 } qw(nucmer_ref nucmer_query  nucmer_prefix nucmer2bed_file);
	for my $param (sort keys %params) {	die "Invalid field '$param'" if not $valid{$param};	}	
	my $nucmer_ref = $params{nucmer_ref};
	my $nucmer_query = $params{nucmer_query};
	my $nucmer_prefix = $params{nucmer_prefix};
	my $nucmer2bed_file = $params{nucmer2bed_file};
# 	$params{mode} //= ''; 	
 	die "Parameter 'nucmer_ref' is required" if not $nucmer_ref;
 	die "Parameter 'nucmer_query' is required" if not $nucmer_query;
 	die "Parameter 'nucmer_prefix' is required" if not $nucmer_prefix;
 	die "Parameter 'nucmer2bed_file' is required" if not $nucmer2bed_file;

# 	$params{nucmer_ref} //= '';
# 	$params{nucmer_query} //= '';
# 	$params{nucmer_prefix} //= '';
# 	$params{nucmer2bed_out} //= '';
	
	$bash_cmd = "nucmer -p $nucmer_prefix $nucmer_ref $nucmer_query";
	execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "NO", do_on_fail => "die");

	my $delta_file = $nucmer_prefix.".delta"; # was created by nucmer
	my $coords_file = $nucmer_prefix.".qclTH.coords"; # will be created by show-ccords
# 	my $coords_header = "[strand]". "\t". "[S1]". "\t"."[E1]". "\t"."[S2]". "\t"."[E2]". "\t"."[LEN 1]". "\t"."[LEN 2]". "\t"."[% IDY]". "\t"."[LEN R]". "\t"."[LEN Q]". "\t"."[COV R]". "\t"."[COV Q]". "\t"."[TAGS]";
# 	my $nucmer2bed_file = $coords_file."2bed.bed";


	my $num_delta_lines = getNumLines($delta_file);
	if ($num_delta_lines == 2)	{	
		return $num_delta_lines;
	}
	elsif ($num_delta_lines > 2 ) {
# 		$bash_cmd = "show-coords -qclTH $delta_file > $coords_file"; # -q => Sorted by reference = pCC1_trailer
		$bash_cmd = "show-coords -rclTH $delta_file > $coords_file"; # -q => Sorted by reference = READ/CONTIG
		execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "NO", do_on_fail => "die");
		convert_coords_to_bed($coords_file, $nucmer2bed_file, "REF");
	}
	return $num_delta_lines;	
}

############################################################
sub getNumContigs_in_fasta {
	# Invocation: getNumContigs_in_fasta($thisBarcode_assembly);
	my $in_fasta = $_[0];
 	my $in_fasta_object = Bio::SeqIO->new(-file => "$in_fasta", -format => 'Fasta') or die $!;
 	my $numContigs = 0; 
 	while ( my $contigObj = $in_fasta_object->next_seq) { $numContigs++; }
 	return $numContigs;
}
############################################################
sub get_contig_readCoverage_fromHeader {
	# Invocation: get_contig_readCoverage_fromHeader($contig_description);
	my $in_header = $_[0];
	# Parse: len=46286 reads=159 class=contig suggestRepeat=no suggestBubble=no suggestCircular=yes trim=4728-41540
	my @x = split(" ", $in_header);
	my $numReads = "todo";
	foreach (@x) {
		if ($_ =~ m/^reads=/) {		($numReads = $_) =~ s/reads=//g;	}
	}
 	return $numReads;
}

############################################################
sub write_ReverseComplemented_Fasta {
	# New version : input contig/file gets re-written with its reverse complement
	# Invocatio:
	# write_ReverseComplemented_Fasta($contig_filename);
	my $contig_filename = $_[0];
	# Read the fasta and s
	my $input_genomeObject = Bio::SeqIO->new(-file => "$contig_filename", -format => 'Fasta') or die $!; 
	my $numContigs = 0;
	my $out_name;
	my $out_desc;
	my $out_seq;
	while ( my $contigObj = $input_genomeObject->next_seq) {	
		$numContigs++;
		$out_name = $contigObj->id;
		$out_desc = $contigObj->desc;

		$contigObj = $contigObj->revcom; # 		Reverse complement the object
		$out_seq = $contigObj->seq(); # Get the Reverse complemented sequence
	}
	if ($numContigs > 1 ) {
		print LOGFILE "ERROR: Expected only one contig in write_ReverseComplemented_Fasta($contig_filename). Received numContigs = $numContigs\n";
		print LOGFILE "ERROR: EXIT-ing now.\n";
		exit;
	}
	else {
		my $tmp_rc = $contig_filename."tmprc.fasta";
		open(TMP_FA, ">", $tmp_rc) or die $!;
		print TMP_FA ">", $out_name, " ", $out_desc, " ; reverse-complemented for pCC1 alignment", "\n";
		print TMP_FA $out_seq, "\n";
		close(TMP_FA);
		my $bash_cmd = "mv $tmp_rc $contig_filename";
# 		system($bash_cmd) && print LOGFILE "SUCCESS.\n\n";
		execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "die");
	}
}

##########################################################
sub convert_coords_to_bed	{
	# Invoked as:	
	# convert_coords_to_bed($coords_file, $nucmer2bed_file);
	my $coords_file = $_[0]; # e.g. nucmer_bc27_tig00000001_vs_pcc1leader.qclTH.coords
	my $nucmer2bed_file = $_[1]; # 
	my $chromosome_to_use = $_[2]; # Ref or query columns??
	# ONLY retain the bed lines with full-length match to the pCC1-leader
	
	open (TMP_IN, "<", $coords_file) or die $!;
	open (TMP_OUT, ">", $nucmer2bed_file) or die $!;
	while (my $line = <TMP_IN>) {
		chomp($line);
		my @x = split("\t", $line);
		
# 			1	365	19736	20100	365	365	100.00	365	197300	100.00	0.18	my_pCC1_leader	bc27_tig00000001
# 			1	365	49544	49908	365	365	100.00	365	197300	100.00	0.18	my_pCC1_leader	bc27_tig00000001

# 			2	152	6249	6099	151	151	100.00	338	8100	44.67	1.86	NEB131-A2_readNum36_len338	pCC1.trailer_leader_catted
		my $start; # = $x[2] - 1;
		my $end; # = $x[3] - 1;
		my $chr; # = $x[12];
		my $frag_length; # = $end - $start + 1;
		if ($chromosome_to_use eq "QUERY") {
			$start = $x[2] - 1;
			$end = $x[3] - 1;
			$chr = $x[12];
		}
		elsif ($chromosome_to_use eq "REF") {
			$start = $x[0] - 1;
			$end = $x[1] - 1;
			$chr = $x[11];
		}
		else {
			die "In sub:convert_coords_to_bed, invalid option $chromosome_to_use for chromosome_to_use\n";
		}
		$frag_length = $end - $start + 1;
		if ($frag_length <0 ) {
			print LOGFILE "WARNING: In convert_coords_to_bed: frag_length = $frag_length is less than zero. Should have been positive if the strands were properly matched\n";
			print "WARNING: In convert_coords_to_bed: frag_length = $frag_length is less than zero. Should have been positive if the strands were properly matched\n";
		}
		else {
			print TMP_OUT $chr, "\t", $start, "\t", $end, "\n";
		}	
	}
	close(TMP_IN);
	close(TMP_OUT);
}

##########################################################
sub run_bed_substract {
	# Invokded as:
	# run_bed_substract($nucmer_query_fasta, $nucmer2bed_file, $bed_frags_to_get_file);

	my $nucmer_query_fasta = $_[0];
	my $nucmer2bed_file = $_[1];
	my $bed_frags_to_get_file = $_[2];

	# Step 1 : Create a bed-file representing the entire contig
# 	my $queryBed = $outdir03_tmp_logs_etc."/".$nucmer_query_fasta."2bed.bed"; # Easy cleanup if file is created in outdir03_tmp_logs_etc
	my $queryBed = $nucmer_query_fasta."2bed.bed"; # Easy cleanup if file is created in outdir03_tmp_logs_etc
	open(TMP_BED, ">", $queryBed) or die $!;
 	my $in_genomeObject = Bio::SeqIO->new(-file => "$nucmer_query_fasta", -format => 'Fasta') or die $!; 
 	while ( my $contigObj = $in_genomeObject->next_seq) {
		# There should be only one contig in this fasta file
		my $contig_id = $contigObj->id();
		my $contig_seq = $contigObj->seq();
 		my $contig_len_0indexed = length($contig_seq) - 1;
 		my $bedLine =  $contig_id."\t"."0"."\t".$contig_len_0indexed."\n"; 
 		print TMP_BED $bedLine;
#  		print LOGFILE "DEBUG: Inside sub:run_bed_substract: queryBed = ", $bedLine;
 	}
 	close(TMP_BED);
 	
 	# Step 2 : run bedtools subtract
 	my $bash_cmd = "bedtools subtract -a $queryBed -b $nucmer2bed_file > $bed_frags_to_get_file";
#  	print LOGFILE "bash_cmd = $bash_cmd\n";
#  	system($bash_cmd) && print LOGFILE "bash_cmd = SUCCESS.\n\n";
	execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "die");

 	$bash_cmd = "rm -f $queryBed"; # okay to lose this file
#  	print LOGFILE "bash_cmd = $bash_cmd\n";
#  	system($bash_cmd) && print LOGFILE "bash_cmd = SUCCESS.\n\n";
	execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "die");
}		
##########################################################
		
sub get_fasta_fromfrags_bed {
	# Invocation: get_fasta_fromfrags_bed($bed_frags_to_get_file, $contig_filename, $out_shredded_fasta);
	my $bed_frags_to_get_file = $_[0];
	my $contig_filename = $_[1];
	my $out_shredded_fasta = $_[2];
	
	my $contig_id;
	my $contig_seq;
 	my $contig_len;
# 	print LOGFILE "DEBUG: Inside sub:get_fasta_fromfrags_bed; Processing bed_frags_to_get_file = $bed_frags_to_get_file\n";

 	my $genomeObject = Bio::SeqIO->new(-file => "$contig_filename", -format => 'Fasta') or die $!; 
 	while ( my $contigObj = $genomeObject->next_seq) {
		# There should be only one contig in this fasta file
		$contig_id = $contigObj->id();
		$contig_seq = $contigObj->seq();
 		$contig_len = length($contig_seq);
 	}	
	
	open (BED, "<", $bed_frags_to_get_file) or die $!;
	open (OUT_FASTA, ">", $out_shredded_fasta) or die $!;
	my $num_bedlines = 0;
	while (my $line = <BED>) {
		chomp($line);
		$num_bedlines++;
# 		print LOGFILE "DEBUG: Inside sub:get_fasta_fromfrags_bed; BED : $line\n";
		#	bc01_tig00000001        0       24804
		#	bc01_tig00000001        24805   61612
		#	bc01_tig00000001        61613   98420
		#	bc01_tig00000001        98421   134287
		my ($bedcontig, $start, $end, $fragName)  = split("\t", $line);
# 		print LOGFILE "DEBUG: bedcontig=$bedcontig; start=$start; end=$end, fragName=$fragName\n";
		my $fragLen = $end - $start + 1;
# 		print LOGFILE "DEBUG: Inside sub:get_fasta_fromfrags_bed; fragLen = $fragLen\n";
		$fragName = $bedcontig."_frag_".$num_bedlines."_len_".$fragLen;
		my $fragSeq = substr($contig_seq, $start, $fragLen); # substr EXPR,OFFSET,LENGTH
# 		print LOGFILE "DEBUG: Inside sub:get_fasta_fromfrags_bed; fragSeq = $fragSeq\n";
		## For pCC1_leader that is usually found at the very 5-prime end, skip the first frag if it is too small
		if ($fragLen > 99) { # Hard-coded value: use fasta in MSA only if it is of sufficient length
			print OUT_FASTA ">", $fragName, "\n", $fragSeq, "\n";
		}
		else {
			print LOGFILE "DEBUG: Inside sub:get_fasta_fromfrags_bed: SKIPPPING frag with fragLen = $fragLen (< 99)\n";
		}
	}
	close(OUT_FASTA);
	close(BED);
}

##########################################################
sub get_fasta_fromfrags_bed_AS_A_STRING {
	# Invocation: get_fasta_fromfrags_bed_AS_STRING($bed_frags_to_get_file, $contig_filename)
	my $bed_frags_to_get_file = $_[0];
	my $contig_filename = $_[1];
	
	my $contig_id;
	my $contig_seq;
 	my $contig_len;
 	my $return_megafasta_string = "";
# 	print LOGFILE "DEBUG: Inside sub:get_fasta_fromfrags_bed; Processing bed_frags_to_get_file = $bed_frags_to_get_file\n";

 	my $genomeObject = Bio::SeqIO->new(-file => "$contig_filename", -format => 'Fasta') or die $!; 
 	while ( my $contigObj = $genomeObject->next_seq) {
		# There should be only one contig in this fasta file
		$contig_id = $contigObj->id();
		$contig_seq = $contigObj->seq();
 		$contig_len = length($contig_seq);
 	}	
	
	open (BED, "<", $bed_frags_to_get_file) or die $!;
# 	open (OUT_FASTA, ">", $out_shredded_fasta) or die $!;
	my $num_bedlines = 0;
	while (my $line = <BED>) {
		chomp($line);
		$num_bedlines++;
# 		print LOGFILE "DEBUG: Inside sub:get_fasta_fromfrags_bed; BED : $line\n";
		#	bc01_tig00000001        0       24804
		#	bc01_tig00000001        24805   61612
		#	bc01_tig00000001        61613   98420
		#	bc01_tig00000001        98421   134287
		my ($bedcontig, $start, $end, $fragName)  = split("\t", $line);
# 		print LOGFILE "DEBUG: bedcontig=$bedcontig; start=$start; end=$end, fragName=$fragName\n";
		my $fragLen = $end - $start + 1;
# 		print LOGFILE "DEBUG: Inside sub:get_fasta_fromfrags_bed; fragLen = $fragLen\n";
		$fragName = $bedcontig."_frag_".$num_bedlines."_len_".$fragLen;
		my $fragSeq = substr($contig_seq, $start, $fragLen); # substr EXPR,OFFSET,LENGTH
# 		print LOGFILE "DEBUG: Inside sub:get_fasta_fromfrags_bed; fragSeq = $fragSeq\n";
		## For pCC1_leader that is usually found at the very 5-prime end, skip the first frag if it is too small
		if ($fragLen > 150) { # Hard-coded value: use fasta in MSA only if it is of sufficient length
# 			print OUT_FASTA ">", $fragName, "\n", $fragSeq, "\n";
			$return_megafasta_string = $return_megafasta_string.">".$fragName."\n".$fragSeq."\n";
		}
		else {
			print LOGFILE "DEBUG: Inside sub:get_fasta_fromfrags_bed: SKIPPPING frag with fragLen = $fragLen (< 150 bp)\n";
		}
	}
	close(OUT_FASTA);
	close(BED);
	return $return_megafasta_string;
}

##########################################################
sub get_seqkitStats_from_fasta {
	# Invocation: get_fasta_fromfrags_bed($bed_frags_to_get_file, $contig_filename, $out_shredded_fasta);
	my $in_fasta = $_[0];
	my $numContigs = $_[1];
	my $tmp_stats_file = $in_fasta.".seqkit_stats";
	my $stats_to_return;
	my $bash_cmd = "seqkit stats $in_fasta > $tmp_stats_file";
# 	system($bash_cmd) && print LOGFILE "bash_cmd = SUCCESS.\n\n";
	execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "die");
	open (STATS, "<", $tmp_stats_file) or die $!;
	while (my $line = <STATS>) {
		# file                 format  type  num_seqs  sum_len  min_len  avg_len  max_len
		# barcode01.consensus  FASTA   DNA          1   36,808   36,808   36,808   36,808
		if (!($line =~ m/num_seqs/)) {
			chomp($line);
			my ($filename, $format, $type, $num_seqs, $sum_len, $min_len, $avg_len, $max_len) = split(/\s+/, $line);
# 			print LOGS "DEBUG $filename; $format; $type; $num_seqs; $sum_len; $min_len; $avg_len; $max_len;\n";
			print LOGFILE "DEBUG $filename; $format; $type; $num_seqs; $sum_len; $min_len; $avg_len; $max_len;\n";
			if ($numContigs == 1) {
				$stats_to_return = "num_seqs = 1; total_length = $sum_len; filename=$filename.\n";
			}
			else {
				$stats_to_return = "num_seqs = $num_seqs; total_length = $sum_len; min_len=$min_len; avg_len = $avg_len; max_len= $max_len; filename=$filename.\n";			
			}
		}
	}
	close(STATS);
	`rm -f $tmp_stats_file`;
	return $stats_to_return;
}	

###############################################################
sub dump_to_LOGFILE {
	my $dump = $_[0];
	open (DUMP, "<", $dump) or die $!;
	print LOGFILE "\n";
	while (my $dumpline = <DUMP>) {
		print LOGFILE $dumpline;
	}
	print LOGFILE "\n";
	close(DUMP);
# 	`rm -f $dump`;
}

###############################################################
sub rename_fasta_single {
	# Invocation: rename_fasta_single($in_fasta_file, $in_new_id, $in_new_desc);
	my $in_fasta_file = $_[0];
	my $in_new_id = $_[1];
	my $in_new_desc = $_[2];
	
 	my $genomeObject = Bio::SeqIO->new(-file => "$in_fasta_file", -format => 'Fasta') or die $!; 
 	my $contig_id;
 	my $contig_seq;
 	my $contig_len;
 	my $numContigs = 0;
 	while ( my $contigObj = $genomeObject->next_seq) {
		# There should be only one contig in this fasta file
		$numContigs++;
		$contig_id = $contigObj->id();
		$contig_seq = $contigObj->seq();
 		$contig_len = length($contig_seq);
 	}	
	if ($numContigs > 1 ) {
		print LOGFILE "ERROR: Expected only one contig in rename_fasta_single($in_fasta_file). Received numContigs = $numContigs\n";
		print LOGFILE "ERROR: EXIT-ing now.\n";
		exit;
	}
	else {
		my $tmp_fa = $in_fasta_file."tmp_renamed.fasta";
		open(TMP_FA, ">", $tmp_fa) or die $!;
		print TMP_FA ">", $in_new_id, " len=$contig_len;", " ", $in_new_desc, "\n";
		print TMP_FA $contig_seq, "\n";
		close(TMP_FA);
		my $bash_cmd = "mv $tmp_fa $in_fasta_file";
		system($bash_cmd) && print LOGFILE "SUCCESS.\n\n";
	}
}

###############################################################
sub rename_fasta_V2 {

	my %params = @_;
	my %valid = map { $_ => 1 } qw(mode newName newDesc newPrefix newSuffix);
	for my $param (sort keys %params) {	die "Invalid field '$param'" if not $valid{$param};	}	
	my $mode = $params{mode};
# 	$params{mode} //= ''; 	
 	die "Parameter 'mode' is required" if not $mode;

	$params{newName} //= '';
	$params{newDesc} //= '';
	$params{newPrefix} //= '';
	$params{newSuffix} //= '';
	
	if ($params{mode} eq "replace") {
		my $newName = $params{newName};
		my $newDesc = $params{newDesc};
		print "mode:$mode\n";
		print "newName = $params{newName}\n";
		print "newDesc = $params{newDesc}\n";
# 		print "newPrefix = $newPrefix\n";
# 		print "newSuffix = $newSuffix\n";
	}     
	elsif ($mode eq "add_prefix") {
		my $newPrefix = $params{newPrefix};
		my $newSuffix = $params{newSuffix};
		print "mode:$mode\n";
# 		print "newName = $newName\n";
# 		print "newDesc = $newDesc\n";
		print "newPrefix = $newPrefix\n";
		print "newSuffix = $newSuffix\n";		
	}
	else {
		die "Invalid mode $mode in sub:rename_fasta\n";
	}
}

###############################################################
# local $SIG{__DIE__} = sub {
#     my ($die_message) = @_;
#     print LOGFILE $die_message, "\n";
#     print $die_message, "\n";
# };

###############################################################
sub emboss_merger {
# Invoked as :  my $numContigs_after_merger = emboss_merger(input => $out_split_fasta, output => $out_consensus_fasta_noPCC);
	print LOGFILE "emboss_merger INVOKED.\n\n";
	print "emboss_merger INVOKED.\n\n";

	my %params = @_;
	my %valid = map { $_ => 1 } qw(input output);
	for my $param (sort keys %params) {	die "Invalid field '$param'" if not $valid{$param};	}	
	my $in_fasta = $params{input};
	my $out_fasta = $params{output};
 	die "Parameter 'in_fasta' is required\n" if not $in_fasta;
 	die "Parameter 'out_fasta' is required\n" if not $out_fasta;

	my $bash_cmd = "";
 	my $numContigs = 0;
	
	my $tmp_fasta1 = "tmp_fasta1.".$in_fasta; # suffix is essential to keep the filename unique
	my $tmp_fasta2 = "tmp_fasta2.".$in_fasta; # suffix is essential to keep the filename unique

 	my $tmp_merged_fasta = "tmp_merged.".$in_fasta;
 	my $tmp_merged_aln = "tmp_merged.".$in_fasta.".aln";
#  	my $tmp_merged_fasta_Object = Bio::SeqIO->new(-file => ">$tmp_merged_fasta", -format => 'Fasta') or die $!; 

	# Process the input fasta file
 	my $genomeObject = Bio::SeqIO->new(-file => "$in_fasta", -format => 'Fasta') or die $!; 
 	while ( my $contigObj = $genomeObject->next_seq) {
		$numContigs++;
		if ($numContigs == 1) {
			# first contig, write it to a tmp file for use in next step
		 	my $tmp_fasta1_Object = Bio::SeqIO->new(-file => ">$tmp_fasta1", -format => 'Fasta') or die $!; 
			$tmp_fasta1_Object->write_seq($contigObj);		
		}
		else {
			# Start pairwise merge
 			my $tmp_fasta2_Object = Bio::SeqIO->new(-file => ">$tmp_fasta2", -format => 'Fasta') or die $!; 
			$tmp_fasta2_Object->write_seq($contigObj);
			
			$bash_cmd = "merger $tmp_fasta1 $tmp_fasta2 -outfile $tmp_merged_aln -outseq $tmp_merged_fasta";
			execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "die");

			my $numContigs_after_merger = getNumContigs_in_fasta($tmp_merged_fasta);
			if ($numContigs_after_merger > 1) {
				print LOGFILE "MERGER FAILED for in_fasta = $in_fasta\n";
				return $numContigs_after_merger;
				exit;
			}

			# The merged fasta become fasta1 fopr the next round of mergers
			$bash_cmd = "cp $tmp_merged_fasta $tmp_fasta1";
			execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "die");
		}
 	}
	
	my $numContigs_after_merger = getNumContigs_in_fasta($tmp_merged_fasta);
	if ($numContigs_after_merger > 1) {
		print LOGFILE "MERGER FAILED for in_fasta = $in_fasta\n";
# 		return $numContigs_after_merger;
		exit;
	}
	else {
		$bash_cmd = "mv $tmp_merged_fasta $out_fasta";
		execute_and_log_bash_cmd(bash_cmd => $bash_cmd, print_to_log => "YES", do_on_fail => "die");
# 		return $numContigs_after_merger;
	} 	
	## Cleanup
	if (-e $tmp_fasta1) {	`rm -f $tmp_fasta1`;	}
	if (-e $tmp_fasta2) {	`rm -f $tmp_fasta2`;	}
	if (-e $tmp_merged_fasta) {	`rm -f $tmp_merged_fasta`;	}
	if (-e $tmp_merged_aln) {	`rm -f $tmp_merged_fasta`;	}
	return $numContigs_after_merger;
	
}

##################### END #####################################

