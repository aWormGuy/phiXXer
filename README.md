# phiXXer
phiXXer: Assemble fosmid inserts from nanopore sequencing of phiX amplified products

## Installation
One-time setup: Create a conda environment `phiXXer_env` by running `sh.install_phiXXer_env.sh`

## Running the pipeline
For each new run, the qsub script `job-phiXXer.v0.0.sh` needs to be edited as follows:
1. Prepare and provide a tab-separated file which should have two columns: (1) the nanopore_barcode_numnber, and (2) user-specified fosmid_clone_ID. This filepath should be provided to as the value for the variable `input_metadata` in the qsub script.
2. Update all other `User-specified parameters` (lines 25 to 32) in the qsub script.
3. As of now, 96 samples can be multiplexed using nanopore barcodes. The corresponding data can be analyzed in paralled on SGE by specifying `-t 1-96` (line 10) in this qsub script.
4. Submit the job via `qsub`

## Citing phiXXer
The complete experimental and computational workflows are described in the following preprint: 
- Chuzel L, Sinha A, Cunningham CV, Taron CH. High-throughput nanopore DNA sequencing of large insert fosmid clones directly from bacterial colonies. bioRxiv; 2024. [doi:10.1101/2024.02.05.578990](https://doi.org/10.1101/2024.02.05.578990)


