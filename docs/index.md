<div align="center">

  <h1 style="font-size: 250%">mir-seek 🔬</h1>

  <b><i>An awesome microRNA-sequencing pipeline</i></b><br>
  <a href="https://doi.org/10.5281/zenodo.13834237">
    <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.13834237.svg" alt="DOI">
  </a>
  <a href="https://github.com/OpenOmics/mir-seek/releases">
    <img alt="GitHub release" src="https://img.shields.io/github/v/release/OpenOmics/mir-seek?color=blue&include_prereleases">
  </a>
  <a href="https://hub.docker.com/repository/docker/skchronicles/mir-seek">
    <img alt="Docker Pulls" src="https://img.shields.io/docker/pulls/skchronicles/mir-seek">
  </a><br> 
  <a href="https://github.com/OpenOmics/mir-seek/actions/workflows/main.yaml">
    <img alt="tests" src="https://github.com/OpenOmics/mir-seek/workflows/tests/badge.svg">
  </a>
  <a href="https://github.com/OpenOmics/mir-seek/actions/workflows/docs.yml">
    <img alt="docs" src="https://github.com/OpenOmics/mir-seek/workflows/docs/badge.svg">
  </a>
  <a href="https://github.com/OpenOmics/mir-seek/issues">
    <img alt="GitHub issues" src="https://img.shields.io/github/issues/OpenOmics/mir-seek?color=brightgreen">
  </a>
  <a href="https://github.com/OpenOmics/mir-seek/blob/main/LICENSE">
    <img alt="GitHub license" src="https://img.shields.io/github/license/OpenOmics/mir-seek">
  </a>

  <p>
    This is the home of the pipeline, mir-seek. Its long-term goals: to accurately estimate microRNA expression, to perform quality-control, and to boldly discover known and novel miRNAs like no pipeline before!
  </p>

</div>  


## Overview

Welcome to mir-seek's documentation! This guide is the main source of documentation for users that are getting started with the [miR-sequencing pipeline](https://github.com/OpenOmics/mir-seek/). 

The **`./mir-seek`** pipeline is composed several inter-related sub commands to setup and run the pipeline across different systems. Each of the available sub commands perform different functions: 

<section align="center" markdown="1" style="display: flex; flex-wrap: row wrap; justify-content: space-around;">

!!! inline custom-grid-button ""

    [<code style="font-size: 1em;">mir-seek <b>run</b></code>](usage/run.md)   
    Run the mir-seek pipeline with your input files.

!!! inline custom-grid-button ""

    [<code style="font-size: 1em;">mir-seek <b>unlock</b></code>](usage/unlock.md)  
    Unlocks a previous runs output directory.

</section>

<section align="center" markdown="1" style="display: flex; flex-wrap: row wrap; justify-content: space-around;">


!!! inline custom-grid-button ""

    [<code style="font-size: 1em;">mir-seek <b>install</b></code>](usage/install.md)  
    Download remote reference files locally.


!!! inline custom-grid-button ""

    [<code style="font-size: 1em;">mir-seek <b>cache</b></code>](usage/cache.md)  
    Cache remote software containers locally.  

</section>

**mir-seek** is a comprehensive microRNA-sequencing pipeline. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

The pipeline is compatible with data generated from Illumina short-read sequencing technologies. As input, it accepts a set of single-end FastQ files and can be run locally on a compute instance or on-premise using a cluster. A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM (more coming soon!). A hybrid approach ensures the pipeline is accessible to all users.

Before getting started, we highly recommend reading through the [usage](usage/run.md) section of each available sub command.

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](faq/questions.md) prior to [opening an issue on Github](https://github.com/OpenOmics/mir-seek/issues).

## Contribute

This site is a living document, created for and by members like you. mir-seek is maintained by the members of NCBR and is improved by continous feedback! We encourage you to contribute new content and make improvements to existing content via pull request to our [GitHub repository :octicons-heart-fill-24:{ .heart }](https://github.com/OpenOmics/mir-seek).

## Citation

Please note that more citation formats and releases can be found on the [mir-seek zenodo page](https://doi.org/10.5281/zenodo.13834237).

If you use this software, please cite it as below:  

=== "BibTex"

    ```
    @software{Kuhn_OpenOmics_mir-seek_2024,
      author       = {Skyler, Kuhn and
                      Gary, Zhang},
      title        = {OpenOmics/mir-seek: v0.3.0},
      month        = sep,
      year         = 2024,
      publisher    = {Zenodo},
      version      = {v0.3.0},
      doi          = {10.5281/zenodo.13834238},
      url          = {https://doi.org/10.5281/zenodo.13834238}
    }
    ```

=== "APA"

    ```
    Skyler, K., & Gary, Z. (2024). OpenOmics/mir-seek: v0.3.0 (v0.3.0). Zenodo. https://doi.org/10.5281/zenodo.13834238
    ```

## References

<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**2.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>  
