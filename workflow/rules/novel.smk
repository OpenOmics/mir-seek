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
        map_log   = join(workpath, "novel", "mapper", "cohort", "cohort.mapper.log"),
        map_tsv   = join(workpath, "novel", "mapper", "cohort", "cohort.mapper.tsv"),
        new_log   = join(workpath, "novel", "mapper", "cohort", "cohort.bowtie.log"),
    params:
        rname = "novel_mapper",
        # Minimum read length also
        # used for trimming reads
        min_len = min_read_length,
        bw_index = config['references'][genome]['bowtie1_index'],
        tmpdir   = join(workpath, "novel", "mapper", "cohort"),
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
    > {output.map_log} 2>&1
    """
