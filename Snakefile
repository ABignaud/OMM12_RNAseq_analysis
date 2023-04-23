#!/bin/env snakemake -s

# snakemake --rulegraph | dot -Tsvg > images/rulegraph.svg
# snakemake --dag | dot -Tsvg > images/dag.svg

# This file can be run using snakemake. It was tested on snakemake 5.10.0.
# It orchestrates the analysis of the gene trasncription impact on the HiC map
# from different species.

import numpy as np
import pandas as pd
from os.path import join
from snakemake.utils import validate

# Set parameters.
shell.prefix("set -euo pipefail;")

# LOAD CONFIG FILES
configfile: 'config/config.yaml'

samples = pd.read_csv(
    config['samples'], 
    sep=';', 
    dtype=str,
    comment='#',
).set_index(['library'], drop=False)

# Set directory
OUT_DIR = join(config['base_dir'], config['out_dir'])
TMP = join(config['base_dir'], config['tmp_dir'])
REF_DIR = join(config['base_dir'], config['ref_dir'])
FASTQ_DIR = join(config['base_dir'], config['fastq_dir'])
FIG_DIR = join(config['base_dir'], config['fig_dir'])

library = np.unique(samples.library) 

wildcard_constraints:
    library = '|'.join(library),
    ref = '|'.join(config['ref']),
    reseq = "|".join(config['reseq']),

# Pipeline sub-workflows
include: 'rules/01_rna_processing.smk'
include: 'rules/02_feature_counts.smk'

rule all:
    input:
        # 01 - RNA processing
        expand(
            join(OUT_DIR, 'multiqc', '{ref}_multiqc_report.html'),
            ref = config['ref'],
        ),
        # 02 - FeatreCounts
        expand(
           join(OUT_DIR, 'featureCounts', '{bacteria}_featureCounts_corrected.txt'),
            bacteria = config['bacteria'],
        ),
