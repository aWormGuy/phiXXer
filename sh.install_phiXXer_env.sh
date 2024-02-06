#!/bin/bash

## This script needs to be run only once during installation
conda create --yes --name phiXXer_env canu minimap2 samtools;
	conda activate phiXXer_env;
	conda install --yes -c bioconda perl-bioperl;
	conda install --yes mummer4 cd-hit
	conda install --yes seqkit seqtk;
	conda install --yes emboss;
	conda install --yes -c bioconda gepard;
# 	conda install --yes last;
# 	conda install --yes nextflow;
# 	conda install --yes circlator; 
# 	conda install --yes flye;
conda deactivate;
