# Quality-control rules
# Pre-alignment quality control
rule fastqc_raw:
    """
    Quality-control step to assess sequencing quality of raw data.
    For more information, please visit their documentation:
    https://github.com/s-andrews/FastQC
    @Input:
        Raw FastQ file (scatter)
    @Output:
        FastQC html report and zip file on raw data
    """
    input:
        fq  = join(workpath,"{sample}.R1.fastq.gz"),
    output:
        html     = join(workpath, "fastqc", "{sample}.R1_fastqc.html"),
        zip_file = join(workpath, "fastqc", "{sample}.R1_fastqc.zip")
    params:
        rname  = "rawfqc", 
        outdir = join(workpath, "fastqc"),
        tmpdir = tmpdir,
        adapters = config['references'][genome]['adapters'],
    container: config['images']['mir-seek'],
    threads: int(allocated("threads", "fastqc_raw", cluster))
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT

    # Create fastqc adapters file,
    # file needs to be in TSV format
    paste - - < "{params.adapters}" \\
        | sed 's/^>//g' \\
    > "${{tmp}}/adapters.tsv"

    # Running fastqc with local
    # disk or a tmpdir, fastqc
    # has been observed to lock
    # up gpfs filesystems, adding
    # this on request by HPC staff
    fastqc \\
        -a "${{tmp}}/adapters.tsv" \\
        -t {threads} \\
        -o "${{tmp}}" \\
        {input.fq}
    
    # Copy output files from tmpdir
    # to output directory
    find "${{tmp}}" \\
        -type f \\
        \\( -name '*.html' -o -name '*.zip' \\) \\
        -exec cp {{}} {params.outdir} \\;    
    """


rule fastqc_trim:
    """
    Quality-control step to assess sequencing quality of trimmed data.
    For more information, please visit their documentation:
    https://github.com/s-andrews/FastQC
    @Input:
        Trimmed FastQ file (scatter)
    @Output:
        FastQC html report and zip file on trimmed data
    """
    input:
        fq  = join(workpath, "trim", "{sample}_trimmed.fastq.gz"),
    output:
        html     = join(workpath, "fastqc", "{sample}_trimmed_fastqc.html"),
        zip_file = join(workpath, "fastqc", "{sample}_trimmed_fastqc.zip")
    params:
        rname  = "filtfqc", 
        outdir = join(workpath, "fastqc"),
        tmpdir = tmpdir,
        adapters = config['references'][genome]['adapters'],
    container: config['images']['mir-seek'],
    threads: int(allocated("threads", "fastqc_trim", cluster))
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT

    # Create fastqc adapters file,
    # file needs to be in TSV format
    paste - - < "{params.adapters}" \\
        | sed 's/^>//g' \\
    > "${{tmp}}/adapters.tsv"

    # Running fastqc with local
    # disk or a tmpdir, fastqc
    # has been observed to lock
    # up gpfs filesystems, adding
    # this on request by HPC staff
    fastqc \\
        -a "${{tmp}}/adapters.tsv" \\
        -t {threads} \\
        -o "${{tmp}}" \\
        {input.fq}
    
    # Copy output files from tmpdir
    # to output directory
    find "${{tmp}}" \\
        -type f \\
        \\( -name '*.html' -o -name '*.zip' \\) \\
        -exec cp {{}} {params.outdir} \\; 
    """


# Post-alignment quality control
rule multiqc:
    """
    Reporting step to aggregate QC information across all tools.
    For more information, please visit their documentation: 
    https://multiqc.info/
    @Input:
        Fastp trimmed fastq files (gather)
        FastQC html reports (gather),
        miRDeep Mapper/Bowtie/1.X Alignment files (gather)
    @Output:
        MulitQC html report
    """
    input:
        expand(join(workpath, "trim", "{sample}_fastp_report.json"), sample=samples),
        expand(join(workpath, "fastqc", "{sample}.R1_fastqc.zip"), sample=samples),
        expand(join(workpath, "fastqc", "{sample}_trimmed_fastqc.zip"), sample=samples),
        expand(join(workpath, "mirdeep2", "mapper", "{sample}_mapped.arf"), sample=samples),
        expand(join(workpath, "mirdeep2", "mapper", "{sample}", "{sample}.bowtie.log"), sample=samples),
        join(workpath, "mirdeep2", "mapper", "mirdeep2.mapper.tsv"),
    output:
        html = join(workpath, "reports", "multiqc_report.html"),
    params:
        rname  = "multiqc", 
        wdir   = join(workpath),
        outdir = join(workpath, "reports"),
        conf   = join(workpath, "resources", "multiqc_config.yaml"),
    container: config['images']['mir-seek'],
    threads: int(allocated("threads", "mulitqc", cluster))
    shell: """
    multiqc \\
        --ignore '*/.singularity/*' \\
        -f \\
        -c {params.conf} \\
        --interactive \\
        --outdir {params.outdir} \\
        {params.wdir}
    """