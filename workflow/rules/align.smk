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
        FASTA file of processed reads 
    """
    input:
        reads_fa  = join(workpath, "trim", "{sample}_trimmed_cleaned.fa"),
    output:
        arf       = join(workpath, "mirdeep2", "mapper", "{sample}_mapped.arf"),
        collapsed = join(workpath, "mirdeep2", "mapper", "{sample}_collapsed.fa"),
    params:
        rname = "mapper",
        # Minimum read length also
        # used for trimming reads
        min_len = min_read_length,
        bw_index = config['references'][genome]['bowtie1_index'],
    envmodules: config['tools']['bowtie'],
    threads: int(allocated("threads", "mirdeep2_mapper", cluster)),
    shell: """
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
        -o {threads}
    """