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
rule featureCounts:
    input:
        fasta = join(REF_DIR, '{bacteria}.fa'),
        gff = join(REF_DIR, '{bacteria}.{sens}.gtf'),
        bam = expand(
            join(TMP, 'bam', f'{config["ref"][0]}', '{library}.bam'),
            library=library,
        ),
    output:
        join(OUT_DIR, 'featureCounts', '{bacteria}_featureCounts.{sens}.txt'),
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
        sens_file = join(OUT_DIR, 'featureCounts', '{bacteria}_featureCounts.sens.txt'),
        anti_file = join(OUT_DIR, 'featureCounts', '{bacteria}_featureCounts.antisens.txt'),
    output:
        join(OUT_DIR, 'featureCounts', '{bacteria}_featureCounts_corrected.txt'),
    threads: 1
    conda: '../envs/DGE.yaml'
    script: '../scripts/diff_antisens.py'

