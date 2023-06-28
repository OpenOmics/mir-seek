# Rules related to trimming and pre-procesing
rule fastp:
    """
    Data-processing step to adapter trimming and read filtering. 
    @Input:
        Raw FastQ files (scatter)
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
        adapters = config['references']['adapters']
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
