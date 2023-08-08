# Rules related to alignment to mirbase database
rule mirdeep2_mapper:
    """
    Data-processing step to align reads to the reference genome. 
    For more information, check out its github repository:
    https://github.com/rajewsky-lab/mirdeep2
    @Input:
        Cleaned FASTA file of Reads (scatter)
    @Output:
        Reads mapped to the reference genome (ARF file),
        FASTA file of processed reads,
        mirdeep2 mapper logfile,
        mirdeep2 per-sample alignment stats file,
        bowtie/1.X log file
    """
    input:
        reads_fa  = join(workpath, "trim", "{sample}_trimmed_cleaned.fa"),
    output:
        arf       = join(workpath, "mirdeep2", "mapper", "{sample}_mapped.arf"),
        collapsed = join(workpath, "mirdeep2", "mapper", "{sample}_collapsed.fa"),
        map_log   = join(workpath, "mirdeep2", "mapper", "{sample}", "{sample}.mapper.log"),
        map_tsv   = join(workpath, "mirdeep2", "mapper", "{sample}", "{sample}.mapper.tsv"),
        new_log   = join(workpath, "mirdeep2", "mapper", "{sample}", "{sample}.bowtie.log"),
    params:
        rname = "mapper",
        # Minimum read length also
        # used for trimming reads
        min_len = min_read_length,
        bw_index = config['references'][genome]['bowtie1_index'],
        tmpdir   = join(workpath, "mirdeep2", "mapper", "{sample}"),
        sample   = "{sample}",
    envmodules: config['tools']['bowtie'],
    container: config['images']['mir-seek'],
    threads: int(allocated("threads", "mirdeep2_mapper", cluster)),
    shell: """
    # Setups temporary directory for
    # intermediate files, miRDeep2
    # output directories rely on
    # timestamps, this helps avoid
    # collision due to multiple runs
    # of the same sample, needed for
    # the bowtie/1.X log files
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    cd "${{tmp}}"

    # Aligns reads to the reference 
    # genome with Bowtie/1.X, allows
    # one mismatch in the alignment
    mapper.pl \\
        {input.reads_fa} \\
        -c \\
        -j \\
        -l {params.min_len} \\
        -m \\
        -n \\
        -q \\
        -p {params.bw_index} \\
        -s {output.collapsed} \\
        -t {output.arf} \\
        -v \\
        -o {threads} \\
    > {output.map_log} 2>&1

    # Extract mapper statistics
    paste \\
        <(echo -e "sample\\n{params.sample}") \\
        <(grep -A1 --color=never '^#desc' \\
            {output.map_log} \\
            | tr ' ' '\\t' \\
            | cut -f2-
        ) \\
    > {output.map_tsv}

    # Rename bowtie/1.X log file
    mv "${{tmp}}/bowtie.log" "{output.new_log}"
    """


rule mirdeep2_gather_mapper_statistics:
    """
    Data-processing step to gather mirDeep2 mapper statistics across all 
    samples. These mapping stats will be added to the MultiQC report, as 
    the bowtie/1.X mapping statistics that are included in the report are
    from deduplicated/collapsed reads. These alignment statistics represent
    uniquely aligned reads, meaning multi-mapping reads are only reported 
    or counted once.  
    @Input:
        mirdeep2 per-sample alignment stats file (gather)
    @Output:
        mirdeep2 alignment stats table
    """
    input:
        map_tsvs = expand(
            join(workpath, "mirdeep2", "mapper", "{sample}", "{sample}.mapper.tsv"),
            sample=samples
        ),
    output:
        map_tsv  = join(workpath, "mirdeep2", "mapper", "mirdeep2.mapper.tsv"),
    params:
        rname = "maptable",
    envmodules: config['tools']['bowtie'],
    container: config['images']['mir-seek'],
    threads: int(allocated("threads", "mirdeep2_gather_mapper_statistics", cluster)),
    shell: """
    # Create table from per-sample
    # mideep2 mapper files
    i=0
    for f in {input.map_tsvs}; do
        if [ "$i" -eq 0  ]; then
            # Add header to output file 
            head -1 "${{f}}" \\
            > {output.map_tsv}
        fi
        awk 'NR=="2" {{print}}' "${{f}}" \\
        >> {output.map_tsv}
        i=$((i + 1))
    done
    """
