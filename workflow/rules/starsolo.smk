# Generate a test output
rule star_solo:
    input: 
        read_1 = lambda wildcards: umi_bc_read[wildcards.id],
        read_2 = lambda wildcards: tx_read[wildcards.id]
    output: "results/star_solo/{id}/{id}.SJ.out.tab"
    conda: "../envs/star_env.yaml"
    benchmark: "results/benchmarks/{id}_star_solo.benchmark"
    log: "results/logs/{id}_star_solo.log"
    threads: 64
    params:
        outdir = "results/star_solo/{id}/{id}.",
        star_index = config['star_index'],
        ref_gtf = config['ref_gtf'],
        bc_whitelist = config['bc_whitelist']
    shell:
        '''
        STAR --runMode alignReads \
        --runThreadN {threads} \
        --genomeDir {params.star_index} \
        --sjdbGTFfile {params.ref_gtf} \
        --readFilesIn {input.read_2} \
        {input.read_1} \
        --readFilesCommand gunzip -c \
        --outFileNamePrefix {params.outdir} \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes CB UB sS GN NH HI AS nM NM MD jM jI MC ch \
        --twopassMode Basic \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist {params.bc_whitelist} \
        --soloUMIlen 12 \
        --soloUMIdedup 1MM_CR \
        --soloUMIfiltering MultiGeneUMI_CR \
        --soloCellFilter None \
        --soloCellReadStats Standard 2> {log}
        '''