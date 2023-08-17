#!/bin/env snakemake -s

# Rules to generates the gene expression analysis.

# Transform ID key values for featureCounts
rule gff_to_gtf:
    input:
        gff = join(REF_DIR, '{bacteria}.gff'),
    output:
        gtf = join(REF_DIR, '{bacteria}.gtf'),
    threads: 1
    conda: '../envs/DGE.yaml'
    shell: "sed 's/ID/gene_id/g' {input.gff} > {output.gtf}"  

# Make the antisens annotation to have two different tables to correct the DNA 
# contamination and create one cleaned table. 
rule gtf_antisens:
    input:
        gtf = join(REF_DIR, '{bacteria}.gtf'),
    output:
        anti = join(REF_DIR, '{bacteria}.antisens.gtf'),
        sens = join(REF_DIR, '{bacteria}.sens.gtf'),
    threads: 1
    conda: '../envs/DGE.yaml'
    shell: 
        """
        sed 's/+/____/' {input.gtf} | \
            sed 's/-/+/' | \
            sed 's/____/-/' > {output.anti}
        ln -s {input.gtf} {output.sens}
        """  

# Launch featureCounts
rule featureCounts_pe:
    input:
        fasta = join(REF_DIR, '{bacteria}.fa'),
        gff = join(REF_DIR, '{bacteria}.{sens}.gtf'),
        bam = expand(
            join(TMP, 'bam', f'{config["ref"][0]}', '{library}{pe_reseq}.bam'),
            library=library,
            pe_reseq=config['pe_reseq'],
        ),
        bam2 = join(TMP, 'bam', f'{config["ref"][0]}', 'DCXXX_nxq2.bam'),
    output:
        join(OUT_DIR, 'featureCounts', '{bacteria}_pe_featureCounts.{sens}.txt'),
    params:
        outdir = join(OUT_DIR, 'featureCounts'),
        tmp = join(OUT_DIR, 'featureCounts', '{bacteria}_featureCounts.{sens}.tmp.txt'),
        file_path = join(TMP, 'bam', f'{config["ref"][0]}/'),
    threads: config['threads']
    conda: '../envs/DGE.yaml'
    shell:
        """
        mkdir -p {params.outdir}
        featureCounts \
            -p \
            -O \
            -s 2 \
            -t CDS \
            --ignoreDup \
            -T {threads} \
            -G {input.fasta} \
            -a {input.gff} \
            -o {params.tmp} \
            {input.bam} {input.bam2}
        FILE_PATH={params.file_path}
        sed -i "s|$FILE_PATH||g" {params.tmp}
        sed 's/.bam//g' {params.tmp} > {output}
        rm {params.tmp}
        """

rule featureCounts_se:
    input:
        fasta = join(REF_DIR, '{bacteria}.fa'),
        gff = join(REF_DIR, '{bacteria}.{sens}.gtf'),
        bam = expand(
            join(TMP, 'bam', f'{config["ref"][0]}', '{library}{se_reseq}.bam'),
            library=library,
            se_reseq=config['se_reseq'],
        ),
    output:
        join(OUT_DIR, 'featureCounts', '{bacteria}_se_featureCounts.{sens}.txt'),
    params:
        outdir = join(OUT_DIR, 'featureCounts'),
        tmp = join(OUT_DIR, 'featureCounts', '{bacteria}_featureCounts.{sens}.tmp.txt'),
        file_path = join(TMP, 'bam', f'{config["ref"][0]}/'),
    threads: config['threads']
    conda: '../envs/DGE.yaml'
    shell:
        """
        mkdir -p {params.outdir}
        featureCounts \
            -O \
            -s 2 \
            -t CDS \
            --ignoreDup \
            -T {threads} \
            -G {input.fasta} \
            -a {input.gff} \
            -o {params.tmp} \
            {input.bam}
        FILE_PATH={params.file_path}
        sed -i "s|$FILE_PATH||g" {params.tmp}
        sed 's/.bam//g' {params.tmp} > {output}
        rm {params.tmp}
        """

# Rule to correct DNA contamination based on antisens transcripts which should 
# be mainly DNA.
rule correct_conta:
    input:
        sens_file = join(OUT_DIR, 'featureCounts', '{bacteria}_{end}_featureCounts.sens.txt'),
        anti_file = join(OUT_DIR, 'featureCounts', '{bacteria}_{end}_featureCounts.antisens.txt'),
    output:
        join(OUT_DIR, 'featureCounts', '{bacteria}_{end}_featureCounts_corrected.txt'),
    threads: 1
    conda: '../envs/DGE.yaml'
    script: '../scripts/diff_antisens.py'

# Join PE and SE tables
# rule join_featurecounts:
#     input:
#         se = join(OUT_DIR, 'featureCounts', '{bacteria}_se_featureCounts_corrected.txt'),
#         pe = join(OUT_DIR, 'featureCounts', '{bacteria}_pe_featureCounts_corrected.txt'),
#     output:
#         join(OUT_DIR, 'featureCounts', '{bacteria}_featureCounts_corrected.txt'),
#     threads: 1
#     conda: '../envs/DGE.yaml'
#     shell: '../scripts/diff_antisens.py'
# awk -v OFS=" "  '{$2=$3=$4=$5=$6=""; print $0}' b_coccoides_se_featureCounts_corrected.txt | sed 's/      / /g' | sed 's/ /\t/g' > tmp 
# join b_coccoides_pe_featureCounts_corrected.txt tmp  | sed 's/ /\t/g' > b_coccoides_featureCounts_corrected.txt 
# rm tmp

# Run DESeq2
rule deseq2:
    input:
        table = join(OUT_DIR, 'featureCounts', '{bacteria}_featureCounts_corrected.txt'),
        metadata = join(config['base_dir'], 'OMM12_RNAseq_analysis', 'config', 'metadata.tsv'),
    params:
        pval = config['pval'],
        fc = config['fold_change'],
        script = join(config['base_dir'], 'OMM12_RNAseq_analysis', 'scripts', 'deseq.r'),
        outdir = join(OUT_DIR, 'DESeq2', '{bacteria}'),
    output:
        touch(join(OUT_DIR, 'DESeq2', '{bacteria}', 'DESeq2.done')),
    threads: 1
    conda: '../envs/DGE.yaml'
    shell: 
        """
        module add R/4.3.0 
        mkdir -p {params.outdir}
        cd {params.outdir}
        Rscript {params.script} \
            {input.metadata} \
            {input.table} \
            {params.fc} \
            {params.pval}
        """

