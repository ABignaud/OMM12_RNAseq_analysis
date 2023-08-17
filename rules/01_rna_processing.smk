#!/bin/env snakemake -s

# Rule to generates the RNA tracks

# FastQC
rule rna_fastqc:
    input:
        join(FASTQ_DIR, '{library}{reseq}_R{end}.fq.gz'),
    output:
        join(OUT_DIR, 'fastqc', '{library}{reseq}_R{end}_fastqc.html'),
    params:
        outdir = join(OUT_DIR, 'fastqc')
    threads: config['threads']
    conda: "../envs/qc.yaml"
    shell: 'fastqc --quiet --threads {threads} --outdir {params.outdir} {input}'        

# Trim adapters
rule rna_trimgalore_pe:
    input:
        R1 = join(FASTQ_DIR, '{library}{pe_reseq}_R1.fq.gz'),
        R2 = join(FASTQ_DIR, '{library}{pe_reseq}_R2.fq.gz'),
    output:
        R1 = join(TMP, 'trimgalore', '{library}{pe_reseq}_R1.fq.gz'),
        R2 = join(TMP, 'trimgalore', '{library}{pe_reseq}_R2.fq.gz'),
    params:
        basename = '{library}{pe_reseq}',
        outdir = join(TMP, 'trimgalore'),
        out_R1 = join(TMP, 'trimgalore', '{library}{pe_reseq}_val_1.fq.gz'),
        out_R2 = join(TMP, 'trimgalore', '{library}{pe_reseq}_val_2.fq.gz'),
    threads: config['threads_large']
    conda: "../envs/qc.yaml"
    shell: 
        """
        trim_galore \
            --basename {params.basename} \
            --output_dir {params.outdir} \
            --illumina \
            --quality 20 \
            --stringency 3 \
            --fastqc \
            --paired \
            -j {threads} \
            --gzip \
            {input.R1} \
            {input.R2}
        mv {params.out_R1} {output.R1}
        mv {params.out_R2} {output.R2}
        """

rule rna_trimgalore_se:
    input:
        R1 = join(FASTQ_DIR, '{library}{se_reseq}_R1.fq.gz'),
    output:
        R1 = join(TMP, 'trimgalore', '{library}{se_reseq}_R1.fq.gz'),
    params:
        basename = '{library}{se_reseq}',
        outdir = join(TMP, 'trimgalore'),
        out_R1 = join(TMP, 'trimgalore', '{library}{se_reseq}_trimmed.fq.gz'),
    threads: config['threads']
    conda: "../envs/qc.yaml"
    shell: 
        """
        trim_galore \
            --basename {params.basename} \
            --output_dir {params.outdir} \
            --illumina \
            --quality 20 \
            --stringency 3 \
            -j {threads} \
            --gzip \
            {input.R1} 
        mv {params.out_R1} {output.R1}
        """

# fastq screen
rule rna_fastq_screen:
    input:
        fastq = join(TMP, 'trimgalore', '{library}{reseq}_R{end}.fq.gz'),
    output:
        join(OUT_DIR, 'fastq_screen', '{library}{reseq}_R{end}_screen.html'),
    params:
        outdir =  join(OUT_DIR, 'fastq_screen'),
        conf = config['fastq_screen_conf']
    threads: config['threads']
    conda: "../envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        fastq_screen \
            --outdir {params.outdir} \
            --conf {params.conf} \
            --threads {threads} \
            --force \
            {input.fastq} 
        """

# Build bowtie2 index of the reference genome.
rule bt2_index:
    input: join(REF_DIR, '{ref}.fa'),
    output: touch(join(REF_DIR, 'index', '{ref}.bt2_index.done')),
    params:
        idx = join(REF_DIR, 'index', '{ref}'),
    threads: config['threads']
    conda: "../envs/gen_tracks.yaml"
    shell: "bowtie2-build --threads {threads} {input} {params.idx}"

# Map the reads and sort them.
rule align_rna_pe:
    input:
        index_flag = join(REF_DIR, 'index', '{ref}.bt2_index.done'),
        R1 = join(TMP, 'trimgalore', '{library}{pe_reseq}_R1.fq.gz'),
        R2 = join(TMP, 'trimgalore', '{library}{pe_reseq}_R2.fq.gz'),
    output: 
        join(TMP, 'bam', '{ref}', '{library}{pe_reseq}.bam')
    params:
        index = join(REF_DIR, 'index', '{ref}'),
        bt2_presets = config['bowtie2_args'],
        outdir = join(TMP, 'align', '{ref}'),
        sam = join(TMP, 'align', '{ref}', '{library}{pe_reseq}.sam'),
        tmp1 = join(TMP, 'align', '{ref}', '{library}{pe_reseq}.bam.tmp1'),
        tmp2 = join(TMP, 'align', '{ref}', '{library}{pe_reseq}.bam.tmp2'),
    threads: config['threads_large']
    conda: "../envs/gen_tracks.yaml"
    log: join(OUT_DIR, 'logs', 'rna',  '{ref}', '{library}{pe_reseq}.log')
    shell:
        """
        mkdir -p {params.outdir}
        bowtie2 {params.bt2_presets} \
                -p {threads} \
                -x {params.index} \
                --maxins 1000 \
                -1 {input.R1} \
                -2 {input.R2} 2> {log} > {params.sam}
        samtools sort -@ {threads} -n -O BAM {params.sam} -o {params.tmp1}
        samtools fixmate -@ {threads} --output-fmt bam -m {params.tmp1} {params.tmp2}
        samtools view -@ {threads} --output-fmt bam -f 2 -q 10 -1 -b {params.tmp2} -o {params.tmp1}
        samtools sort -@ {threads} --output-fmt bam -l 9 {params.tmp1} -o {output}
        rm {params.sam} {params.tmp1} {params.tmp2}
        """

# Map the reads and sort them.
rule align_rna_se:
    input:
        index_flag = join(REF_DIR, 'index', '{ref}.bt2_index.done'),
        R1 = join(TMP, 'trimgalore', '{library}{se_reseq}_R1.fq.gz'),
    output: 
        join(TMP, 'bam', '{ref}', '{library}{se_reseq}.bam')
    params:
        index = join(REF_DIR, 'index', '{ref}'),
        bt2_presets = config['bowtie2_args'],
        outdir = join(TMP, 'align', '{ref}'),
        sam = join(TMP, 'align', '{ref}', '{library}{se_reseq}.sam'),
        tmp1 = join(TMP, 'align', '{ref}', '{library}{se_reseq}.bam.tmp1'),
        tmp2 = join(TMP, 'align', '{ref}', '{library}{se_reseq}.bam.tmp2'),
    threads: config['threads']
    conda: "../envs/gen_tracks.yaml"
    log: join(OUT_DIR, 'logs', 'rna',  '{ref}', '{library}{se_reseq}.log')
    shell:
        """
        mkdir -p {params.outdir}
        bowtie2 {params.bt2_presets} \
                -p {threads} \
                -x {params.index} \
                -U {input.R1} 2> {log} > {params.sam}
        samtools view -@ {threads} --output-fmt bam -q 10 -1 -b {params.sam} -o {params.tmp1}
        samtools sort -@ {threads} --output-fmt bam -l 9 {params.tmp1} -o {output}
        rm {params.sam} {params.tmp1}
        """

# Preseq QC:
rule preseq:
    input:
        bam = join(TMP, 'bam', '{ref}', '{library}{reseq}.bam')
    output:
        preseq = join(OUT_DIR, 'preseq', '{ref}_{library}{reseq}_preseq.txt'),
    threads: 1
    conda: "../envs/qc.yaml"
    shell: 'preseq lc_extrap -v -B {input.bam} -o {output.preseq}'

# Build RNA tracks
rule rna_coverage:
    input:
        bam = join(TMP, 'bam', '{ref}', '{library}{reseq}.bam'),
    output:
        unstranded = join(OUT_DIR, 'RNA_tracks', '{library}{reseq}_{ref}_unstranded.bw'),
        fw = join(OUT_DIR, 'RNA_tracks', '{library}{reseq}_{ref}_forward.bw'),
        rv = join(OUT_DIR, 'RNA_tracks', '{library}{reseq}_{ref}_reverse.bw'),
    threads: config['threads']
    conda: "../envs/gen_tracks.yaml"
    shell:
        """
        samtools index {input.bam} -@ {threads}
        bamCoverage --bam {input.bam} \
            --outFileName {output.unstranded} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM \
            --extendReads 100
        bamCoverage --bam {input.bam} \
            --outFileName {output.fw} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM \
            --extendReads 100 \
            --filterRNAstrand forward
        bamCoverage --bam {input.bam} \
            --outFileName {output.rv} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM \
            --extendReads 100 \
            --filterRNAstrand reverse
        """

# Generate gff of the the mixed bacteria.
rule join_annotation:
    input:
        gff = expand(
            join(REF_DIR, '{bacteria}.gff'),
            bacteria=config['bacteria'],
        ),
    output:
        gff = join(REF_DIR, '{ref}.gff'),
        bed = join(REF_DIR, '{ref}.bed'),
    threads: 1
    conda: "../envs/qc.yaml"
    shell: 
        """
        echo "##gff-version 3" > {output.gff}
        for i in {input.gff}
            do grep "##sequence-region" $i >> {output.gff}
        done
        for i in {input.gff}
            do grep -v "^#" $i >> {output.gff}
        done
        grep "CDS" {output.gff} | \
            awk -v OFS='\t' '{{print $1,$4,$5,$3,"0",$7,$4,$5,"255,0,0","1",$5-$4,"0"}}' \
            > {output.bed}
        """

# RSeQC 
rule rseqc:
    input:
        bam = join(TMP, 'bam', '{ref}', '{library}{reseq}.bam'),
        bw = join(OUT_DIR, 'RNA_tracks', '{library}{reseq}_{ref}_unstranded.bw'),
        bed = join(REF_DIR, '{ref}.bed'), 
    output:
        duprate = join(OUT_DIR, 'RSeQC', '{ref}_{library}{reseq}.read_duplication.DupRate_plot.pdf'),
        genebody = join(OUT_DIR, 'RSeQC', '{ref}_{library}{reseq}.geneBodyCoverage.pdf'),
    params:
        basename =  join(OUT_DIR, 'RSeQC', '{ref}_{library}{reseq}'),
    threads: 1
    conda: "../envs/qc.yaml"
    shell:
        """
        samtools index {input.bam} -@ {threads}
        infer_experiment.py -i {input.bam} -r {input.bed} \
            > {params.basename}.infer_experiment.txt
        bam_stat.py -i {input.bam} 2> {params.basename}.bam_stat.txt
        inner_distance.py \
            -i {input.bam} \
            -o {params.basename}.rseqc \
            -r {input.bed}
        read_distribution.py -i {input.bam} -r {input.bed} \
            > {params.basename}.read_distribution.txt
        read_duplication.py \
            -i {input.bam} \
            -o {params.basename}.read_duplication
        geneBody_coverage2.py \
            --input-file {input.bw} \
            --refgene {input.bed} \
            --out-prefix {params.basename}
        """

# MultiQC report
rule multiQC_rna_report:
    input:
        join(OUT_DIR, 'RSeQC', '{ref}_DCXXX_nxq2.geneBodyCoverage.pdf'),
        join(OUT_DIR, 'preseq', '{ref}_DCXXX_nxq2_preseq.txt'),
        expand(
            join(OUT_DIR, 'fastqc', 'DCXXX_nxq2_R{end}_fastqc.html'),
            end=[1, 2],
        ),
        expand(
            join(OUT_DIR, 'fastq_screen', 'DCXXX_nxq2_R{end}_screen.html'),
            end=[1, 2],
        ),
        expand(
            join(OUT_DIR, 'RSeQC', '{ref}_{library}{reseq}.geneBodyCoverage.pdf'),
            library=library,
            reseq=config['reseq'],
            ref=config['ref'],
        ),
        expand(
            join(OUT_DIR, 'preseq', '{ref}_{library}{reseq}_preseq.txt'),
            library=library,
            reseq=config['reseq'],
            ref=config['ref'],
        ),
        expand(
            join(OUT_DIR, 'fastqc', '{library}{pe_reseq}_R{end}_fastqc.html'),
            library=library,
            pe_reseq = config['pe_reseq'],
            end = [1, 2],
        ),
        expand(
            join(OUT_DIR, 'fastqc', '{library}{se_reseq}_R1_fastqc.html'),
            library=library,
            se_reseq = config['se_reseq'],
        ),
        expand(
            join(OUT_DIR, 'fastq_screen', '{library}{pe_reseq}_R{end}_screen.html'),
            library=library,
            pe_reseq = config['pe_reseq'],
            end = [1, 2],
        ),
        expand(
            join(OUT_DIR, 'fastq_screen', '{library}{se_reseq}_R1_screen.html'),
            library=library,
            se_reseq = config['se_reseq'],
        ),
    output:
        join(OUT_DIR, 'multiqc', '{ref}_multiqc_report.html'),
    params:
        title = '{ref}',
        outdir = join(OUT_DIR, 'multiqc'),
        fastq_dir = FASTQ_DIR,
        fastqc_dir = join(OUT_DIR, 'fastqc'),
        fastq_screen_dir = join(OUT_DIR, 'fastq_screen'),
        rseqc_dir = join(OUT_DIR, 'RSeQC'),
        trim_galore_dir = join(TMP, 'trimgalore'),
        preseq_dir = join(OUT_DIR, 'preseq'),
    threads: 1
    conda: "../envs/qc.yaml"
    shell:
        """
        multiqc \
            --title {params.title} \
            --outdir {params.outdir} \
            --force \
            --module fastqc \
            --module fastq_screen \
            --module rseqc \
            --module cutadapt \
            --module preseq \
            ../fastq \
            ../results/fastqc \
            ../results/fastq_screen/ \
            ../results/RSeQC/ \
            ../tmp/trimgalore \
            ../results/preseq
        """
