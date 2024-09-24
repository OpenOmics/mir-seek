# <code>mir-seek <b>run</b></code>

## 1. About 
The `mir-seek` executable is composed of several inter-related sub commands. Please see `mir-seek -h` for all available options.

This part of the documentation describes options and concepts for <code>mir-seek <b>run</b></code> sub command in more detail. With minimal configuration, the **`run`** sub command enables you to start running mir-seek pipeline. 

Setting up the mir-seek pipeline is fast and easy! In its most basic form, <code>mir-seek <b>run</b></code> only has *three required inputs*.

## 2. Synopsis
```text
$ mir-seek run [--help] \
      [--dry-run] [--job-name JOB_NAME] [--mode {slurm,local}] \
      [--sif-cache SIF_CACHE] [--singularity-cache SINGULARITY_CACHE] \
      [--silent] [--threads THREADS] [--tmp-dir TMP_DIR] \
      [--min-read-length MIN_READ_LENGTH] \
      [--max-read-length MAX_READ_LENGTH] \
      [--novel-mir-identification] \
      --input INPUT [INPUT ...] \
      --output OUTPUT \
      --genome {hg38,mm10}
```

The synopsis for each command shows its arguments and their usage. Optional arguments are shown in square brackets.

A user **must** provide a list of single-end FastQ files (globbing is supported) to analyze via `--input` argument, an output directory to store results via `--output` argument, and a reference genome for alignment and annotation via `--genome` argument.

Use you can always use the `-h` option for information on a specific sub command. 

### 2.1 Required arguments

Each of the following arguments are required. Failure to provide a required argument will result in a non-zero exit-code.

  `--input INPUT [INPUT ...]`  
> **Input FastQ file(s).**  
> *type: file(s)*  
> 
> One or more single-end FastQ files can be provided. The pipeline does NOT support paired-end data. From the command-line, each input file should seperated by a space. Globbing is supported! This makes selecting FastQ files easy. Please note the input FastQ files should always be gzipp-ed.
> 
> ***Example:*** `--input .tests/*.R1.fastq.gz`

---  
  `--output OUTPUT`
> **Path to an output directory.**   
> *type: path*
>   
> This location is where the pipeline will create all of its output files, also known as the pipeline's working directory. If the provided output directory does not exist, it will be created automatically.
> 
> ***Example:*** `--output /data/$USER/mir-seek_out`

---  
  `--genome {hg38,mm10}`
> **Reference genome.**   
> *type: string*
>   
> This option defines the reference genome of the samples. mir-seek does comes bundled with pre-built reference files from GENCODE and miRBase (v22) for human and mouse samples. Please select from one of the following options: `hg38`, `mm10`. Please note that `hg38` is a human reference genome, while `mm10` is a reference genomes available for mouse.
> ***Example:*** `Example: --genome hg38`

### 2.2 Analysis options

Each of the following arguments are optional, and do not need to be provided. 


  `--min-read-length MIN_READ_LENGTH`   
> **Minimum read length.**  
> *type: int*  
> *default: 17*
> 
> After trimming adapters, reads shorter than this length will be discarded. ENCODE discards reads shorter than 16 bp; however, miRDeep2 enforces a minimum read-length of 17 bp. If you feel 17 bp is too permissive, setting this option's value to `18` is also a good alternative. 
> 
> ***Example:*** `---min-read-length 17`

---  
  `--max-read-length MAX_READ_LENGTH`   
> **Maximum read length.**  
> *type: int*  
> *default: 27*
> 
> After trimming adapters, reads that exceed this length will be trimmed again from their 3'-end to this exact length. This option may result in the recovery of additional aligned reads, as miRDeep2's aligner only allows for  1 mismatch. This means if a trimmed read's length is 26 bps and it perfectly aligns to a 24 bp micro-RNA, then it will remain unaligned (2 mismatches); however, if this option were set to 25, then it would not be discarded. If you do not want to cropped any of the  reads, you can set this options value to a higher  number, like 999. If you want to maximize the number of aligned reads at the cost of accuracy/precision, then you can set this option to a lower number, like 25. Please take care before overriding the default  value to this option. The maximum sequence length of mature miRNAs in miRBase (v22) for hsa/human is 28, while the maximum sequence length of hairpin miRNAs miRBase (v22) for hsa/human is 180. The maximum sequence length of mature miRNA in miRBase (v22) in mmu/mouse is 27, while the maximum sequence length of hairpin miRNAs miRBase (v22) for mmu/mouse is 147.
> 
> ***Example:*** `---max-read-length 27`

---  
  `--novel-mir-identification`  
> **Quantify novel microRNAs.**  
> *type: boolean flag*  
>  
> This option will enable the pipeline to identify novel microRNAs.  If this option is provided, the pipeline will run mirdeep2 using a two-pass approach to create a counts matrix of novel miR expression. In the first-pass,  novel miRs are identified using information across all samples. In the second-pass, the expression of each novel is quantified using *quantifier.pl* and a counts matrix is produced.  If you are interested in identifying and quantifying novel microRNAs, then please provide this option. Please note that this option will increase the runtime of the pipeline.  
>  
> ***Example:*** `--novel-mir-identification`  

### 2.3 Orchestration options

Each of the following arguments are optional, and do not need to be provided. 

  `--dry-run`            
> **Dry run the pipeline.**  
> *type: boolean flag*
> 
> Displays what steps in the pipeline remain or will be run. Does not execute anything!
>
> ***Example:*** `--dry-run`

---  
  `--silent`            
> **Silence standard output.**  
> *type: boolean flag*
> 
> Reduces the amount of information directed to standard output when submitting master job to the job scheduler. Only the job id of the master job is returned.
>
> ***Example:*** `--silent`

---  
  `--mode {slurm,local}`  
> **Execution Method.**  
> *type: string*  
> *default: slurm*
> 
> Execution Method. Defines the mode or method of execution. Vaild mode options include: slurm or local. 
> 
> ***slurm***    
> The slurm execution method will submit jobs to the [SLURM workload manager](https://slurm.schedmd.com/). It is recommended running mir-seek in this mode as execution will be significantly faster in a distributed environment. This is the default mode of execution.
>
> ***local***  
> Local executions will run serially on compute instance. This is useful for testing, debugging, or when a users does not have access to a high performance computing environment. If this option is not provided, it will default to a local execution mode. 
> 
> ***Example:*** `--mode slurm`

---  
  `--job-name JOB_NAME`  
> **Set the name of the pipeline's master job.**  
> *type: string*
> *default: pl:mir-seek*
> 
> When submitting the pipeline to a job scheduler, like SLURM, this option always you to set the name of the pipeline's master job. By default, the name of the pipeline's master job is set to "pl:mir-seek".
> 
> ***Example:*** `--job-name pl_id-42`

---  
  `--singularity-cache SINGULARITY_CACHE`  
> **Overrides the $SINGULARITY_CACHEDIR environment variable.**  
> *type: path*  
> *default: `--output OUTPUT/.singularity`*
>
> Singularity will cache image layers pulled from remote registries. This ultimately speeds up the process of pull an image from DockerHub if an image layer already exists in the singularity cache directory. By default, the cache is set to the value provided to the `--output` argument. Please note that this cache cannot be shared across users. Singularity strictly enforces you own the cache directory and will return a non-zero exit code if you do not own the cache directory! See the `--sif-cache` option to create a shareable resource. 
> 
> ***Example:*** `--singularity-cache /data/$USER/.singularity`

---  
  `--sif-cache SIF_CACHE`
> **Path where a local cache of SIFs are stored.**  
> *type: path*  
>
> Uses a local cache of SIFs on the filesystem. This SIF cache can be shared across users if permissions are set correctly. If a SIF does not exist in the SIF cache, the image will be pulled from Dockerhub and a warning message will be displayed. The `mir-seek cache` subcommand can be used to create a local SIF cache. Please see `mir-seek cache` for more information. This command is extremely useful for avoiding DockerHub pull rate limits. It also remove any potential errors that could occur due to network issues or DockerHub being temporarily unavailable. We recommend running mir-seek with this option when ever possible.
> 
> ***Example:*** `--singularity-cache /data/$USER/SIFs`

---  
  `--threads THREADS`   
> **Max number of threads for each process.**  
> *type: int*  
> *default: 2*
> 
> Max number of threads for each process. This option is more applicable when running the pipeline with `--mode local`.  It is recommended setting this vaule to the maximum number of CPUs available on the host machine.
> 
> ***Example:*** `--threads 12`

---  
  `--tmp-dir TMP_DIR`   
> **Max number of threads for each process.**  
> *type: path*  
> *default: `/lscratch/$SLURM_JOBID`*
> 
> Path on the file system for writing temporary output files. By default, the temporary directory is set to '/lscratch/$SLURM_JOBID' for backwards compatibility with the NIH's Biowulf cluster; however, if you are running the pipeline on another cluster, this option will need to be specified. Ideally, this path should point to a dedicated location on the filesystem for writing tmp files. On many systems, this location is set to somewhere in /scratch. If you need to inject a variable into this string that should NOT be expanded, please quote this options value in single quotes.
> 
> ***Example:*** `--tmp-dir /scratch/$USER/`

### 2.4 Miscellaneous options  
Each of the following arguments are optional, and do not need to be provided. 

  `-h, --help`            
> **Display Help.**  
> *type: boolean flag*
> 
> Shows command's synopsis, help message, and an example command
> 
> ***Example:*** `--help`

## 3. Example
```bash 
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./mir-seek run --input .tests/*.fastq.gz \
                  --output /data/$USER/output \
                  --genome hg38 \
                  --novel-mir-identification \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the mir-seek pipeline
# The slurm mode will submit jobs to 
# the cluster. It is recommended running 
# the pipeline in this mode.
./mir-seek run --input .tests/*.fastq.gz \
                  --output /data/$USER/output \
                  --genome hg38 \
                  --novel-mir-identification \
                  --mode slurm
```