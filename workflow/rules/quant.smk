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
        mirna     = join(workpath, "mirdeep2", "counts", "{sample}_miRNA_expressed.tsv"),
    log: 
        report    = join(workpath, "mirdeep2", "{sample}_mirdeep2.log")
    params:
        rname   = "mirdeep2",
        fasta   = config['references'][genome]['genome'],
        mature  = config['references'][genome]['mature'],
        hairpin = config['references'][genome]['hairpin'],
        species = config['references'][genome]['species'],
        tmpdir  = join(workpath, "mirdeep2", "run", "{sample}"),
    envmodules: config['tools']['bowtie'],
    threads: int(allocated("threads", "mirdeep2_run", cluster)),
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
        find "${{tmp}}/expression_analyses/" \\
            -type f \\
            -iname "miRNA_expressed.csv" \\
            -print \\
            -quit
    )
    ln -sf "${{exp}}" {output}
    """


rule mature_expression:
    """
    Data-processing step to calculate avergae expression across the same mature miRNA.
    The relationship between mature to precursor miRNA is 1:many. This step averages the
    expression of multiple precursor miRNA pointing to the same mature miRNA to find
    average mature miRNA expression.
    @Input:
        Known and novel miRNA expression (scatter)
    @Ouput:
        Average mature miRNA expression
    """
    input:
        mirna   = join(workpath, "mirdeep2", "counts", "{sample}_miRNA_expressed.tsv"),
    output:
        avg_exp = join(workpath, "mirdeep2", "counts", "{sample}_mature_miRNA_expression.tsv"),
    params:
        rname = "matrexp",
    threads: int(allocated("threads", "mature_expression", cluster)),
    shell: """
    # Removes comment character from 
    # header and calculates average 
    # mature miRNA expression
    head -1 {input.mirna} \\
        | sed '1 s/^#//g' \\
        | cut -f1,2 \\    
    > {output.avg_exp}
    # Cut on prefix of miRBase identifer
    # to get mature miRNA identifers for 
    # aggregation/averaging. These identifers
    # more compatible with downstream tools.
    tail -n+2 {input.mirna} \\
        | cut -f1,2 \\
        | awk -F '\\t' -v OFS='\\t' '{{split($1,a,"MIMA"); print a[1], $NF}}' \\
        | awk -F '\\t' -v OFS='\\t' '{{seen[$1]+=$2; count[$1]++}} END {{for (x in seen) print x, seen[x]/count[x]}}' \\
    >> {output.avg_exp}
    """
