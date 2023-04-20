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

wildcard_constraints:
    library = '|'.join(np.unique(samples.library))
    reseq = "|".join(['isq']),

# Pipeline sub-workflows
include: 'rules/01_rna_processing.smk'

rule all:
    input:
        # 01 - RNA processing
        expand(
            join(TMP, 'salmon', '{ref}', '{rna_library}{reseq}.done'),
            rna_library=rna_library,
            reseq = '_nxq',
            ref=config['ref'],      
        ),
        join(OUT_DIR, 'multiqc', f'{config["title"]}_multiqc_report.html'),
