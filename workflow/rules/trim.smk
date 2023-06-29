# Rules related to trimming and pre-procesing
rule fastp:
    """
    Data-processing step to adapter trimming and read filtering.
    For more information, check out its github repository:
    https://github.com/OpenGene/fastp
    @Input:
        Raw FastQ file (scatter)
    @Output:
        Trimmed FastQ file and reports
    """
    input:
        raw_fq  = join(workpath,"{sample}.R1.fastq.gz"),
    output:
        trim_fq     = join(workpath,"trim", "{sample}_trimmed.fastq.gz"),
        html_report = join(workpath,"trim", "{sample}_fastp_report.html"),
        json_report = join(workpath,"trim", "{sample}_fastp_report.json"),
    params:
        rname = "trim",
        adapters = config['references'][genome]['adapters'],
        min_len = min_read_length,
        max_len = max_read_length,
    envmodules: config['tools']['fastp'],
    threads: int(allocated("threads", "fastp", cluster)),
    shell: """
    fastp \\
        --thread {threads} \\
        --in1 {input.raw_fq} \\
        --out1 {output.trim_fq} \\
        --json {output.json_report} \\
        --html {output.html_report} \\
        --adapter_fasta {params.adapters} \\
        -l {params.min_len} \\
        --max_len1 {params.max_len}
    """


rule seqkit_fq2fa:
    """Data processing step to convert trimmed FastQ file to FASTA format.
    miRDeep2 is very particular about the format of the input FASTA files.
    As so, they need to be cleaned of any white-space, tab, an asterick
    characters. For more information, check out its github repository:
    https://github.com/shenwei356/seqkit
    @Input:
        Trimmed FastQ file (scatter)
    @Output:
        Cleaned FASTA file of Reads 
    """
    input:
        trim_fq  = join(workpath,"trim", "{sample}_trimmed.fastq.gz")
    output:
        raw_fa   = join(workpath,"trim", "{sample}_trimmed.fa"),
        clean_fa = join(workpath,"trim", "{sample}_trimmed_cleaned.fa"),
    params:
        rname = "fq2fa",
    envmodules: config['tools']['seqkit'],
    threads: int(allocated("threads", "seqkit_fq2fa", cluster)),
    shell: """
    # Covert FastQ to FASTA format
    seqkit fq2fa --threads {threads} \\
        {input.trim_fq} \\
        -o {output.raw_fa} 
    
    # Clean sequence identifiers 
    # to replace spaces, tabs, and 
    # asterisks with underscores
    sed '/^>/ s/\\s/_/g' {output.raw_fa} \\
        | sed '/^>/ s/\\t/_/g' \\
        | sed '/^>/ s/*/_/g' \\
    > {output.clean_fa}
    """
