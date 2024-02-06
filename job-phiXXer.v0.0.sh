#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 4 
#$ -l mem_free=128G 
# export all environment variables to SGE
#$ -V
# Parallel jobs to fire : specified by "-t" below
#$ -t 1-96

# set -ue; ## Commented out as it creates problem with mkdir in parallel mode


echo "Running script $0 on `hostname`";
echo "Running in folder `pwd`";
#echo "Job is:"
################################################
#cat $0;
################################################

################################################
### User-specified parameters: Compulsory ######
################################################
PHIXXER_SRC_FOLDER="."; # Folder where the phiXXer*.pl script and pCC1 fasta reeference files are
PHIXXER_SCRIPT=$PHIXXER_SRC_FOLDER"/phiXXer.v0.0.pl";

input_raw_reads_folder="nanopore_out_raw"; # Demultiplexed raw reads from nanopore run
input_metadata="barcode-to-cloneID-map.tsv"; # Table with nanopore barcode_number and the corresponding clone_ID

pCC1_trailer_leader_fasta=$PHIXXER_SRC_FOLDER"/pCC1.trailer_leader_catted.fasta";
my_DEBUG_mode="OFF"; # OFF by defualt; Change to "ON" for debugging purposes

################################################
################################################
my_output_dir="phixxer_${input_raw_reads_folder}"; # Folder where the outputs will be stored

### User-specified parameters: Optional to change values from default
CANU_MIN_COVERAGE=7; # Default CANU_MIN_COVERAGE=7
echo;  echo "CANU_MIN_COVERAGE = $CANU_MIN_COVERAGE (Used as values for canu parameters stopOnLowCoverage and minInputCoverage)."; echo;

#### Actual script : Begin
my_job_name="job-phiXXer.v0.0.sh";
# my_job_name=$0;
echo;echo "BEGIN $my_job_name : `date`";echo;

NUM_CPU=4;

current_BARCODE=$((SGE_TASK_ID));
echo "current_BARCODE = ${current_BARCODE}";

echo "my_output_dir = $my_output_dir";
[ -d "$my_output_dir" ] || mkdir $my_output_dir;

CMD="md5sum $PHIXXER_SCRIPT";
echo;echo "Running: $CMD [`date`]";eval ${CMD};echo;

source activate phiXXer_env;

perl $PHIXXER_SCRIPT --barcodeNum=${current_BARCODE} \
        --raw_reads_folder=${input_raw_reads_folder} \
        --output_folder=${my_output_dir} \
        --in_pCC1_trailer_leader=${pCC1_trailer_leader_fasta} \
        --CANU_MIN_COVERAGE=${CANU_MIN_COVERAGE} \
        --my_RUN_metadata_tsv=${input_metadata} \
        --my_DEBUG_mode=${my_DEBUG_mode};

echo;echo "$my_job_name DONE: `date`";echo;
############### END OF SCRIPT #################################

