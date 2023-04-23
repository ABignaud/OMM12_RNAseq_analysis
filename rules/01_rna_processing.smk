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
rule rna_trimgalore:
    input:
        R1 = join(FASTQ_DIR, '{library}{reseq}_R1.fq.gz'),
        R2 = join(FASTQ_DIR, '{library}{reseq}_R2.fq.gz'),
    output:
        R1 = join(TMP, 'trimgalore', '{library}{reseq}_R1.fq.gz'),
        R2 = join(TMP, 'trimgalore', '{library}{reseq}_R2.fq.gz'),
    params:
        basename = '{library}{reseq}',
        outdir = join(TMP, 'trimgalore'),
        out_R1 = join(TMP, 'trimgalore', '{library}{reseq}_val_1.fq.gz'),
        out_R2 = join(TMP, 'trimgalore', '{library}{reseq}_val_2.fq.gz'),
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
            --fastqc \
            --paired \
            -j {threads} \
            --gzip \
            {input.R1} \
            {input.R2}
        mv {params.out_R1} {output.R1}
        mv {params.out_R2} {output.R2}
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
rule align_rna:
    input:
        index_flag = join(REF_DIR, 'index', '{ref}.bt2_index.done'),
        R1 = join(TMP, 'trimgalore', '{library}{reseq}_R1.fq.gz'),
        R2 = join(TMP, 'trimgalore', '{library}{reseq}_R2.fq.gz'),
    output: 
        join(TMP, 'align', '{ref}', '{library}{reseq}.bam'),
    params:
        index = join(REF_DIR, 'index', '{ref}'),
        bt2_presets = config['bowtie2_args'],
        outdir = join(TMP, 'align', '{ref}'),
        sam = join(TMP, 'align', '{ref}', '{library}{reseq}.sam'),
        tmp1 = join(TMP, 'align', '{ref}', '{library}{reseq}.bam.tmp1'),
        tmp2 = join(TMP, 'align', '{ref}', '{library}{reseq}.bam.tmp2'),
    threads: config['threads']
    conda: "../envs/gen_tracks.yaml"
    log: join(OUT_DIR, 'logs', 'rna',  '{ref}', '{library}{reseq}.log')
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
        samtools sort -@ {threads} -O BAM {params.tmp2} -o {params.tmp1}
        samtools markdup -@ {threads} --output-fmt bam -r {params.tmp1} {params.tmp2}
        samtools view -@ {threads} --output-fmt bam -f 2 -q 10 -1 -b {params.tmp2} -o {params.tmp1}
        samtools sort -@ {threads} --output-fmt bam -l 9 {params.tmp1} -o {output}
        rm {params.sam} {params.tmp1} {params.tmp2}
        """

# Merge RNA-seq alignements
rule merge_reseq_rna_alignments:
    input: 
        lambda w: [join(
            TMP,
            'align',
            w.ref,
            f'{i.split("_R")[0]}.bam',
        ) for i in list(filter(
            lambda x:w.library in x,  list(filter(
                lambda y: 'R1' in y, os.listdir(FASTQ_DIR)
            ))
        ))]
    output: join(TMP, 'bam', '{ref}', '{library}.bam')
    threads: config['threads']
    conda: "../envs/gen_tracks.yaml"
    shell:
        """
        samtools merge -n -O BAM -@ {threads} - {input} | 
            samtools sort -@ {threads} --output-fmt bam -l 9 -o {output}
        """

# Build index of transcripts using salmon
rule salmon_index:
    input:
        ref = join(REF_DIR, '{ref}.ffn'),
    output:
        touch(join(REF_DIR, 'index', '{ref}.salmon_index.done')),
    params:
        index = join(REF_DIR, '{ref}_salmon_index'),
    threads: config['threads']
    conda: "../envs/gen_tracks.yaml"
    shell: 'salmon index -t {input.ref} -i {params.index}'

# Map librairies using salmon
rule run_salmon:
    input:
        index_flag = join(REF_DIR, '{ref}' + '.salmon_index.done'),
        R1 = join(TMP, 'trimgalore', '{library}{reseq}_R1.fq.gz'),
        R2 = join(TMP, 'trimgalore', '{library}{reseq}_R2.fq.gz'),
    output:
        touch(join(TMP, 'salmon', '{ref}', '{library}{reseq}.done'))
    params:
        index = join(REF_DIR, '{ref}_salmon_index'),
        outdir = join(TMP, 'salmon', '{ref}', '{library}{reseq}')
    threads: config['threads']
    conda: "../envs/gen_tracks.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        salmon quant \
            --libType ISF \
            --index {params.index} \
            --mates1 {input.R1} \
            --mates2 {input.R2} \
            --threads {threads} \
            --output {params.outdir}
        """

# Preseq QC:
rule preseq:
    input:
        bam = join(TMP, 'bam', '{ref}', '{library}.bam'),
    output:
        preseq = join(OUT_DIR, 'preseq', '{ref}_{library}_preseq.txt'),
    threads: 1
    conda: "../envs/qc.yaml"
    shell: 'preseq lc_extrap -v -B {input.bam} -o {output.preseq}'

# Build RNA tracks
rule rna_coverage:
    input:
        bam = join(TMP, 'bam', '{ref}', '{library}.bam'),
    output:
        unstranded = join(OUT_DIR, 'RNA_tracks', '{library}_{ref}_unstranded.bw'),
        fw = join(OUT_DIR, 'RNA_tracks', '{library}_{ref}_forward.bw'),
        rv = join(OUT_DIR, 'RNA_tracks', '{library}_{ref}_reverse.bw'),
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
            --extendReads
        bamCoverage --bam {input.bam} \
            --outFileName {output.fw} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM \
            --extendReads \
            --filterRNAstrand forward
        bamCoverage --bam {input.bam} \
            --outFileName {output.rv} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM \
            --extendReads \
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
        bam = join(TMP, 'bam', '{ref}', '{library}.bam'),
        bw = join(OUT_DIR, 'RNA_tracks', '{library}_{ref}_unstranded.bw'),
        bed = join(REF_DIR, '{ref}.bed'), 
    output:
        duprate = join(OUT_DIR, 'RSeQC', '{ref}_{library}.read_duplication.DupRate_plot.pdf'),
        genebody = join(OUT_DIR, 'RSeQC', '{ref}_{library}.geneBodyCoverage.pdf'),
    params:
        basename =  join(OUT_DIR, 'RSeQC', '{ref}_{library}'),
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
        expand(
            join(OUT_DIR, 'RSeQC', '{ref}_{library}.geneBodyCoverage.pdf'),
            library=library,
            ref=config['ref'],
        ),
        expand(
            join(OUT_DIR, 'preseq', '{ref}_{library}_preseq.txt'),
            library=library,
            ref=config['ref'],
        ),
        expand(
            join(OUT_DIR, 'fastqc', '{library}{reseq}_R{end}_fastqc.html'),
            library=library,
            reseq = config['reseq'],
            end = [1, 2],
        ),
        expand(
            join(OUT_DIR, 'fastq_screen', '{library}{reseq}_R{end}_screen.html'),
            library=library,
            reseq = config['reseq'],
            end = [1, 2],
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
