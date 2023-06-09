# mir-seek pipeline
# Author: Yue (Gary) Zhang, Ph.D
from snakemake.utils import R
import glob
import os
from os import listdir
from os.path import isfile, isdir, join
import os.path as path
import sys

# Set the paths
path_root = config['path_root']
path_raw_reads = config['path_raw_reads']
path_analysis = config['path_analysis']


# Creates a list of raw fastq files 
sample_list = [os.path.basename(file).split('_R1')[0] for file in glob.glob(path_raw_reads+'/*R1_001.fastq.gz')]                  # Basenames of the input samples
print(sample_list)


# Final targets of pipeline
rule all:
    input:
        # run step 1 trim
        expand("{path}/Trimmed/{sample}_trimmed.fastq.gz", path=path_analysis, sample=sample_list),
        expand("{path}/Trimmed/{sample}_trimmed_curated.fa", path=path_analysis, sample=sample_list),
        expand("{path}/mirdeep2/temporary1/{sample}_collapsed.fa", path=path_analysis, sample=sample_list),
        expand("{path}/mirdeep2/temporary1/{sample}_mapperARF.arf", path=path_analysis, sample=sample_list),
        expand("{path}/mirdeep2/{sample}_miRNA_expressed.tsv", path=path_analysis, sample=sample_list),
        expand("{path}/mirdeep2/{sample}_miRNA_expressed_hashtagRmd_mature.tsv", path=path_analysis, sample=sample_list),


# 1. trim raw data
rule trim:
    input:
        path.join(path_raw_reads,"{sample}_R1_001.fastq.gz"),
    output:
        path.join(path_analysis,"Trimmed/{sample}_trimmed.fastq.gz"),
    envmodules:
        config['tools']['trimmomatic'],
        config['tools']['java']
    log: 
        "logs/{sample}_trim.log",
    params:
        rname = "trim",
        tag = "{sample}",
        t = "4",
        adapters = config['references']['adapters']
    shell: """
        /usr/bin/java -jar /usr/local/apps/trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar \
        SE -threads {params.t} -phred33 {input} {output} \
        ILLUMINACLIP:{params.adapters}:2:36:10 \
        LEADING:10 TRAILING:10 MAXINFO:50:0.97 MINLEN:18 CROP:27 > {log} 2>&1
    """


# 2. convert fq to fa AND remove special characters
rule fq2fa:
    input:
        path.join(path_analysis,"Trimmed/{sample}_trimmed.fastq.gz"),
    output:
        path.join(path_analysis,"Trimmed/{sample}_trimmed_curated.fa"),
    envmodules:
        config['tools']['seqkit'],
    params:
        rname = "fq2fa",
        tag = "{sample}",
        t = "4",
        fa = path.join(path_analysis,"Trimmed/{sample}_trimmed.fa")
    shell: """
        seqkit fq2fa {input} -o {params.fa} --threads {params.t}
        # Rename sequence identifiers to replace whitespace characters and asterisks with underscore
        sed '/^>/ s/\s/_/g' {params.fa} | sed '/^>/ s/\t/_/g' | sed '/^>/ s/*/_/g' > {output}
    """


# 3. miRDeep2 step1: mapper
rule mapper:
    input:
        path.join(path_analysis,"Trimmed/{sample}_trimmed_curated.fa"),
    output:
        collapsed = path.join(path_analysis,"mirdeep2/temporary1/{sample}_collapsed.fa"),
        arf = path.join(path_analysis,"mirdeep2/temporary1/{sample}_mapperARF.arf"),
    log: 
        log = "logs/{sample}_mapper.log",
        error = "logs/{sample}_mapper.err",
    envmodules:
        config['tools']['bowtie'],
    params:
        rname = "mapper",
        tag = "{sample}",
        t = "4",
        hg38_index = config['references']['hg38_index']
    shell: """
        mapper.pl {input} -c -j -l 18 -m -q -p {params.hg38_index} \
        -s {output.collapsed} -t {output.arf} \
        -v -o {params.t} 2>{log.log} 1>{log.error}
    """


# 4. miRDeep2 step2: main script
rule miRDeep2:
    input:
        collapsed = path.join(path_analysis,"mirdeep2/temporary1/{sample}_collapsed.fa"),
        arf = path.join(path_analysis,"mirdeep2/temporary1/{sample}_mapperARF.arf"),
    output:
        path.join(path_analysis,"mirdeep2/{sample}_miRNA_expressed.tsv"),
    log: 
        log = path.join(path_root,"workflow/logs/{sample}_mirdeep2.log"),
        error = path.join(path_root,"workflow/logs/{sample}_mirdeep2.err"),
    envmodules:
        config['tools']['bowtie'],
    params:
        rname = "mirdeep2",
        tag = "{sample}",
        hg38_fasta = config['references']['hg38_fasta'],
        mature = config['references']['mature'],
        hairpin = config['references']['hairpin'],
        mir_dir = path.join(path_analysis,"mirdeep2/temporary2/{sample}")
    shell: """
        mkdir -p {params.mir_dir}
        cd {params.mir_dir}
        miRDeep2.pl {input.collapsed} {params.hg38_fasta} {input.arf} \
        {params.mature} none {params.hairpin} -t Human -P -v 2>{log.log} 1>{log.error}
        cp {params.mir_dir}/expression_analyses/expression_analyses_*/miRNA_expressed.csv {output}
    """


# 5. parse mideep2 output
rule parseTable:
    input:
        path.join(path_analysis,"mirdeep2/{sample}_miRNA_expressed.tsv"),
    output:
        path.join(path_analysis,"mirdeep2/{sample}_miRNA_expressed_hashtagRmd_mature.tsv")
    envmodules:
        config['tools']['R']
    params:
        rname = "parseTable",
        tag = "{sample}",
        mir_dir_main = path.join(path_analysis,"mirdeep2"),
        hashtagRmd = path.join(path_analysis,"mirdeep2/{sample}_miRNA_expressed_hashtagRmd.tsv"),
        merged = path.join(path_analysis,"mirdeep2/mirdeep2_counts_merged.csv"),
        collapse_mature_count = path.join(path_root,"workflow/scripts/create_mature_miRNA_count.R"),
        bind_column = path.join(path_root,"workflow/scripts/bindcolumns.R"),
    shell: """
        # remove hashtag at the front of the first line
        # remove suffix in miRNA column and aggregate rows to calculate mean for each miRNA
        tail -c +2 {input} > {params.hashtagRmd}
        Rscript {params.collapse_mature_count} {params.hashtagRmd} {output}
        # merge data
        Rscript {params.bind_column} '{params.mir_dir_main}' 'mature' 'tsv' 'miRNA' 'read_count' {params.merged}
    """
