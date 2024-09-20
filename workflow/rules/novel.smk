# Rules related to novel miR prediction and
# quantification. miRDeep2 does not produce 
# a joint counts matrix with novel miR counts.
# It does produce a strangely formatted file
# with several counts matrices; however, the 
# identifers produced do not match across
# different samples. To produce a joint set 
# of novel miRs (with matching identifers),
# we run mirdeep2 in two passes.
# Here is a general overview of the approach:
#   0. Concatenate all the trimmed reads 
#      across all samples with `cat`, 
#      (gather-per-cohort)
#   1. Run mirdeep2 `mapper.pl` on the 
#      concatenated trimmed reads fasta file
#      to collapse reads and align them against 
#      the reference genome, (singleton)
#   2. Run mirdeep2 `miRDeep2.pl` on the mapper
#      output to get a fasta file of novel miRs.
#   3. Run mirdeep2 `quantifier.pl` for each
#      sample with the cohort's novel miR
#      fasta file to get novel miR counts, 
#      (scatter-per-sample)

# First-pass rules for novel miR quantification
rule mirdeep2_novel_p1_concatenate:
    """Data processing step to concatenate the cleaned
    reads across all samples. These cleaned reads are 
    used to generate a joint callset of novel miRs. 
    Please see this issue for more information:
    https://github.com/rajewsky-lab/mirdeep2/issues/120
    @Input:
        Cleaned FASTA file of Reads (gather-per-cohort)
    @Output:
        FASTA file of all reads across all samples 
    """
    input:
        fastas = expand(
            join(workpath,"trim", "{sample}_trimmed_cleaned.fa"),
            sample=samples
        ),
    output:
        concat_fa = join(workpath, "novel", "cohort_trimmed_cleaned.fa"),
    params:
        rname = "novel_concat",
        trimdir =  join(workpath, "trim"),
    envmodules: config['tools']['seqkit'],
    container: config['images']['mir-seek'],
    threads: int(allocated("threads", "mirdeep_novel_p1_concatenate", cluster)),
    shell: """
    # Concatenate cleaned, trimmed 
    # reads across all samples, using
    # find with cat to avoids ARG_MAX
    # issue which limits max length
    # of a given command. 
    find "{params.trimdir}" \\
        -type f \\
        -iname '*_trimmed_cleaned.fa' \\
        -exec cat {{}} \\; \\
    > {output.concat_fa}  
    """


rule mirdeep2_novel_p1_mapper:
    """
    Data-processing step to align concatenated reads to the 
    reference genome. This step run mirdeep2 `mapper.pl` on 
    the concatenated trimmed reads fasta file to generate a 
    joint callset of novel miRs.
    Please see this issue for more information:
    https://github.com/rajewsky-lab/mirdeep2/issues/120
    @Input:
        FASTA file of all reads across all samples (gather-per-cohort)
    @Output:
        Reads mapped to the reference genome (ARF file),
        FASTA file of processed reads,
        mirdeep2 mapper logfile,
        mirdeep2 per-sample alignment stats file,
        bowtie/1.X log file
    """
    input:
        reads_fa  = join(workpath, "novel", "cohort_trimmed_cleaned.fa"),
    output:
        arf       = join(workpath, "novel", "mapper", "cohort_mapped.arf"),
        collapsed = join(workpath, "novel", "mapper", "cohort_collapsed.fa"),
    params:
        rname = "novel_mapper",
        # Minimum read length also
        # used for trimming reads
        min_len = min_read_length,
        bw_index = config['references'][genome]['bowtie1_index'],
        tmpdir   = join(workpath, "novel", "mapper", "cohort"),
        map_log = join(workpath, "novel", "mapper", "cohort", "cohort.mapper.log"),
    envmodules: config['tools']['bowtie'],
    container: config['images']['mir-seek'],
    threads: int(allocated("threads", "mirdeep2_novel_p1_mapper", cluster)),
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
    > {params.map_log} 2>&1

    # Delete bowtie2 log so it does
    # not get picked up by multiQC
    rm -f "${{tmp}}/bowtie.log"
    """


rule mirdeep2_novel_p1_run:
    """
    Data-processing step to detect novel miRNA expression using
    the joint alignments of the cleaned reads across all samples.
    This step runs mirdeep2 `miRDeep2.pl` on the pass-1 mapper 
    output to generate a joint callset of novel miRs. This step
    produces a FASTA file of novel miRs which can be used to
    quantify novel miRs in each sample using `quantifier.pl`.
    Please see this issue for more information:
    https://github.com/rajewsky-lab/mirdeep2/issues/120
    @Input:
        All collapsed reads mapped to the reference genome (singleton),
        FASTA file of processed reads 
    @Ouput:
        A joint callset of novel mature miRNAs in FASTA format,
        A joint callset of novel hairpin miRNAs in FASTA format
    """
    input:
        arf       = join(workpath, "novel", "mapper", "cohort_mapped.arf"),
        collapsed = join(workpath, "novel", "mapper", "cohort_collapsed.fa"),
    output:
        mature    = join(workpath, "novel", "pass1", "cohort_novel_mature_miRNA.tsv"),
        hairpin   = join(workpath, "novel", "pass1", "cohort_novel_hairpin_miRNA.tsv"),
    log: 
        report    = join(workpath, "novel", "pass1", "mirdeep2.log")
    params:
        rname   = "novel_run",
        fasta   = config['references'][genome]['genome'],
        mature  = config['references'][genome]['mature'],
        hairpin = config['references'][genome]['hairpin'],
        # Building miRDeep2 -t species option,
        # To get a list of supported species names,
        # please run the following command: 
        #   $ miRDeep2.pl -u
        # This is not a required option to miRDeep2
        # so if your organism is not on the list, then
        # please set the "species" key in the following
        # file, config/genome.json, to an empty string.
        # This will ensure mirDeep2 is run without this
        # option to avoid any errors associated with 
        # providing an invalid species name.
        species_option = lambda _: "-t {0}".format(
            config['references'][genome]['species']
        ) if config['references'][genome]['species'] else "",
        tmpdir  = join(workpath, "novel", "pass1", "cohort"),
    envmodules: config['tools']['bowtie'],
    container: config['images']['mir-seek'],
    threads: int(allocated("threads", "mirdeep2_novel_p1_run", cluster)),
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

    # Run miRDeep2 get a joint callset
    # of novel miRNAs across all samples
    miRDeep2.pl \\
        {input.collapsed} \\
        {params.fasta} \\
        {input.arf} \\
        {params.mature} \\
        none \\
        {params.hairpin} \\
        -P {params.species_option} \\
        -v \\
    2> {log.report}

    # Get the novel miRNA mature and hairpin
    # sequences from the miRDeep2 output
    mature=$(
        find "${{tmp}}/" \\
            -type f \\
            -iname "novel_mature_*.fa" \\
            -print \\
            -quit
    )
    star=$(
        find "${{tmp}}/" \\
            -type f \\
            -iname "novel_star_*.fa" \\
            -print \\
            -quit
    )
    # Rename the identifier to contain the
    # suffix _star to indicate if its a star
    # sequence before merging the two files
    sed -i '/^>.*/s/$/_star/' "${{star}}"
    cat "${{mature}}" "${{star}}" > {output.mature}
    
    # Create a symbolic link to the novel 
    # hairpin sequences for quantification
    hairpin=$(
        find "${{tmp}}/" \\
            -type f \\
            -iname "novel_pres_*.fa" \\
            -print \\
            -quit
    )
    ln -sf "${{hairpin}}" {output.hairpin}
    """


# Second-pass rules for novel miR quantification
rule mirdeep2_novel_p2_quantifier:
    """
    Data-processing step to quantify novel miRNAs in each sample
    using the novel miR fasta file that was generated in the first
    pass. This step runs mirdeep2 `quantifier.pl` on the novel miR
    fasta file to get novel miR counts.
    Please see this issue for more information:
    https://github.com/rajewsky-lab/mirdeep2/issues/120
    @Input:
        FASTA file of collapsed reads mapped (scatter-per-sample),
        Reads mapped to the reference genome (scatter-per-sample),
        FASTA file of novel mature miRNAs,
        FASTA file of novel hairpin miRNAs
    @Output:
        Novel miRNA expression
    """
    input:
        arf       = join(workpath, "mirdeep2", "mapper", "{sample}_mapped.arf"),
        collapsed = join(workpath, "mirdeep2", "mapper", "{sample}_collapsed.fa"),
        mature    = join(workpath, "novel", "pass1", "cohort_novel_mature_miRNA.tsv"),
        hairpin   = join(workpath, "novel", "pass1", "cohort_novel_hairpin_miRNA.tsv"),
    output:
        mirna     = join(workpath, "novel", "counts", "{sample}_novel_miRNA_expressed.tsv"),
    params:
        rname   = "novel_quantifier",
        tmpdir  = join(workpath, "novel", "counts", "{sample}"),
    envmodules: config['tools']['bowtie'],
    container: config['images']['mir-seek'],
    threads: int(allocated("threads", "mirdeep2_novel_p2_quantifier", cluster)),
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

    # Run miRDeep2 to quantify novel
    # miRNAs in each sample
    quantifier.pl \\
        -p {input.hairpin} \\
        -m {input.mature} \\
        -r {input.collapsed} \\
        -T {threads}
    
    # Link novel expression results from 
    # miRDeep2 timestamp directory
    exp=$(
        find "${{tmp}}/" \\
            -type f \\
            -iname "miRNA_expressed.csv" \\
            -print \\
            -quit
    )
    ln -sf "${{exp}}" {output.mirna}
    """


rule mirdeep2_novel_p2_mature_expression:
    """
    Data-processing step to calculate avergae expression across the same mature novel miRNA.
    The relationship between mature to precursor novel miRNA is 1:many. This step averages 
    the expression of multiple precursor miRNA pointing to the same mature miRNA to find
    average novel mature miRNA expression.
    @Input:
        Known novel miRNA expression (scatter)
    @Ouput:
        Average novel mature miRNA expression
    """
    input:
        mirna   = join(workpath, "novel", "counts", "{sample}_novel_miRNA_expressed.tsv"),
    output:
        avg_exp = join(workpath, "novel", "counts", "{sample}_novel_mature_miRNA_expression.tsv"),
    params:
        rname = "novel_matrexp",
    container: config['images']['mir-seek'],
    threads: int(allocated("threads", "mirdeep2_novel_p2_mature_expression", cluster)),
    shell: """
    # Removes comment character from 
    # header and calculates average 
    # mature miRNA expression
    head -1 {input.mirna} \\
        | sed '1 s/^#//g' \\
        | cut -f1,2 \\
    > {output.avg_exp}
    # Find average expression due to
    # 1:many relationship between 
    # mature and precursor miRNA
    tail -n+2 {input.mirna} \\
        | cut -f1,2 \\
        | awk -F '\\t' -v OFS='\\t' '{{seen[$1]+=$2; count[$1]++}} END {{for (x in seen) print x, seen[x]/count[x]}}' \\
    >> {output.avg_exp}
    """


rule mirdeep2_novel_p2_merge_results:
    """
    Data-processing step to aggreagte per-sample mature miRNA 
    counts into a counts matrix.
    @Input:
        Average novel mature miRNA expression files (gather)
    @Output:
        Average novel mature miRNA counts matrix
    """
    input:
        counts = expand(join(workpath, "novel", "counts", "{sample}_novel_mature_miRNA_expression.tsv"), sample=samples),
    output:
        matrix = join(workpath, "novel", "counts", "miRDeep2_novel_mature_miRNA_counts.tsv"),
    params:
        rname   = "novel_mirmatrix",
        script  = join("workflow", "scripts", "create_matrix.py"),
    container: config['images']['mir-seek'],
    threads: int(allocated("threads", "mirdeep2_novel_p2_merge_results", cluster)),
    shell: """
    # Create counts matrix of mature miRNAs
    {params.script} \\
        --input {input.counts} \\
        --output {output.matrix} \\
        --join-on miRNA \\
        --extract read_count \\
        --clean-suffix '_novel_mature_miRNA_expression.tsv' \\
        --nan-values 0.0
    """