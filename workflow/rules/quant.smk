# Rules for quantification of known and novel miRNAs
rule mirdeep2_run:
    """
    Data-processing step to detect known and novel miRNA expression.
    For more information, check out its github repository:
    https://github.com/rajewsky-lab/mirdeep2
    @Input:
        Reads mapped to the reference genome (scatter),
        FASTA file of processed reads 
    @Ouput:
        Known and novel miRNA expression
    """
    input:
        arf       = join(workpath, "mirdeep2", "mapper", "{sample}_mapped.arf"),
        collapsed = join(workpath, "mirdeep2", "mapper", "{sample}_collapsed.fa"),
    output:
        mirna     = join(workpath, "mirdeep2", "{sample}_miRNA_expressed.tsv"),
    log: 
        report    = join(workpath, "mirdeep2", "{sample}_mirdeep2.log")
    params:
        rname = "mirdeep2",
        fasta   = config['references'][genome]['fasta'],
        mature  = config['references'][genome]['mature'],
        hairpin = config['references'][genome]['hairpin'],
        species = config['references'][genome]['species'],
        tmpdir = join(worthpath, "mirdeep2", "run", "{sample}"),
    envmodules: config['tools']['bowtie'],
    shell: """
    # Setups temporary directory for
    # intermediate files, miRDeep2
    # output directories rely on
    # timestamps, this helps avoid
    # collision due to multiple runs
    # of the same sample
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    cd "${{tmp}}"

    # Run miRDeep2 to detect known
    # and novel miRNA expression
    miRDeep2.pl \\
        {input.collapsed} \\
        {params.fasta} \\
        {input.arf} \\
        {params.mature} \\
        none \\
        {params.hairpin} \\
        -t {params.species} \\
        -P \\
        -v \\
    2> {log.report}

    # Link expression results from 
    # miRDeep2 timestamp directory
    exp=$(
        find "${{tmp}}"/expression_analyses/ \\
            -type f \\
            -name "miRNA_expressed.csv" \\
            -print \\
            -quit
    )
    ln -sf "${{exp}}" {output}
    """