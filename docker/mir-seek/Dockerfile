# Base image for mir-seek,
# uses Ubuntu Jammy (LTS)
FROM ubuntu:22.04

# Depedencies of mir-seek:
#  - fastp/0.23.2         # binary from opengene.org, installed: 0.23.4
#  - seqkit/2.2.0         # from github, installed: 2.4.0
#  - bowtie/1.3.1         # from apt-get, installed: 1.3.1-1
#  - fastqc/0.11.9        # apt-get, installed: 0.11.9
#  - miRDeep2             # from github, installed c6440e2
#      @requires: bowtie/1.X
#      @requires: Vienna and RNAfold
#      @requires: SQUID
#      @requires: randfold
#      @requires: PDF::API2 perl package 
#  - multiqc/1.14         # from pypi, installed: 1.14
#      @requires: python3
#      @requires: python3-pip
LABEL maintainer="Skyler Kuhn <kuhnsa AT nih.gov>" \
    base_image="ubuntu:22.04" \
    version="v0.1.0"   \
    software="mir-seek" \
    about.summary="An awesome microRNA-sequencing pipeline" \
    about.home="https://github.com/OpenOmics/mir-seek" \
    about.documentation="https://openomics.github.io/mir-seek/" \
    about.tags="Transcriptomics"

############### INIT ################
# Create Container filesystem specific 
# working directory and opt directories
# to avoid collisions with the host's
# filesystem, i.e. /opt and /data
RUN mkdir -p /opt2 && mkdir -p /data2
WORKDIR /opt2 

# Set time zone to US east coast 
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
        && echo $TZ > /etc/timezone

############### SETUP ################
# This section installs system packages 
# required for your project. If you need 
# extra system packages add them here.
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        apt-transport-https \
        autotools-dev \
        automake \
        autoconf \
        build-essential \
        ca-certificates \
        cmake \
        cpanminus \
        curl \
        gawk \
        gfortran \
        git \
        grep \
        gzip \
        locales \
        libtool \
        libpdf-api2-perl \
        libfont-ttf-perl \
        libcompress-raw-zlib-perl \
        liblapack-dev \
        liblapacke-dev \
        make \
        unzip \
        wget \
        perl \
        pkg-config \
        python3 \
        python3-pip \
        software-properties-common \
        unzip \
        zlib1g-dev \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Make python3 the default interpreter
# and install Python Packages 
# and install multiqc/1.14
RUN ln -sf /usr/bin/python3 /usr/bin/python
RUN pip3 install --upgrade pip \
	&& pip3 install argparse \
    && pip3 install matplotlib \
    && pip3 install multiqc==1.14 \
    && pip3 install pandas \
    && pip3 install scipy

# Perl fix issue
RUN cpanm FindBin Term::ReadLine LWP::Simple

############### INSTALL ################
# Install any bioinformatics tools
# available with apt-get on Ubuntu/22.04
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        # bowtie/1.3.1-1
        bowtie \
        # fastqc/0.11.9
        fastqc \
        # samtools/1.13
        samtools \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
WORKDIR /opt2

############### MANUAL ################
# Install tools from src manually
# Installs fastp/0.23.4, for more info see:
# https://github.com/OpenGene/fastp/
RUN mkdir -p /opt2/fastp/0.23.4/ \
    && wget http://opengene.org/fastp/fastp.0.23.4 -P /opt2/fastp/0.23.4/ \
    && mv /opt2/fastp/0.23.4/fastp.0.23.4 /opt2/fastp/0.23.4/fastp \
    && chmod -R a+rx /opt2/fastp/0.23.4/fastp
ENV PATH="${PATH}:/opt2/fastp/0.23.4"

# Installs seqkit/2.4.0, for more info see:
# https://github.com/shenwei356/seqkit/
RUN mkdir -p /opt2/seqkit/2.4.0/ \
    && wget https://github.com/shenwei356/seqkit/releases/download/v2.4.0/seqkit_linux_amd64.tar.gz -P /opt2/seqkit/2.4.0/ \
    && tar -xvzf /opt2/seqkit/2.4.0/seqkit_linux_amd64.tar.gz -C /opt2/seqkit/2.4.0/ \
    && rm /opt2/seqkit/2.4.0/seqkit_linux_amd64.tar.gz \
    && chmod a+rx /opt2/seqkit/2.4.0/seqkit
ENV PATH="${PATH}:/opt2/seqkit/2.4.0"

# Install randfold/2.0.1, this is a 
# dependency of miRDeep2, that fails 
# to install correctly with miRDeep2's
# install.pl script, for more info see:
# http://bioinformatics.psb.ugent.be/supplementary_data/erbon/nov2003/
RUN mkdir -p /opt2/randfold/src /opt2/randfold/bin \
    && wget http://bioinformatics.psb.ugent.be/supplementary_data/erbon/nov2003/downloads/randfold-2.0.1.tar.gz -P /opt2/randfold/src \
    && tar -xvzf /opt2/randfold/src/randfold-2.0.1.tar.gz -C /opt2/randfold/src \
    && rm /opt2/randfold/src/randfold-2.0.1.tar.gz \
    && mv /opt2/randfold/src/randfold-2.0.1-a7feeeaeba2afe567dbd061b9f4965646386bc98 /opt2/randfold/src/randfold-2.0.1 \
    && chmod -x /opt2/randfold/src/randfold-2.0.1/randfold.c /opt2/randfold/src/randfold-2.0.1/Readme.txt \
    && cd /opt2/randfold/src/randfold-2.0.1/ \
    && tar -xvzf /opt2/randfold/src/randfold-2.0.1/squid.tar.gz \
    && cd /opt2/randfold/src/randfold-2.0.1/squid-1.9g/ \
    && ./configure CFLAGS="-O3" --prefix=/opt2/randfold/bin \
    && make \
    && make install \
    && cd .. \
    && sed -i "s@INCLUDE=-I\.@INCLUDE=-I. -I/opt2/randfold/bin/include@" Makefile \
    && sed -i "s@LIBS=-lm -lsquid@LIBS=-lm -lsquid -L/opt2/randfold/bin/lib@" Makefile \
    && make \
    && mv randfold /opt2/randfold/bin
ENV PATH="${PATH}:/opt2/randfold/bin"
WORKDIR /opt2

# Install ViennaRNA/2.6.2, this is a 
# dependency of miRDeep2, that fails 
# to install correctly with miRDeep2's
# install.pl script, for more info see:
RUN mkdir -p /opt2/vienna/2.6.2 \
    && wget https://github.com/ViennaRNA/ViennaRNA/releases/download/v2.6.2/ViennaRNA-2.6.2.tar.gz -P /opt2/vienna/2.6.2 \
    && tar -xvzf /opt2/vienna/2.6.2/ViennaRNA-2.6.2.tar.gz -C /opt2/vienna/2.6.2 \
    && rm /opt2/vienna/2.6.2/ViennaRNA-2.6.2.tar.gz \
    && cd /opt2/vienna/2.6.2/ViennaRNA-2.6.2 \
    && ./configure  \
    && make \
    && make install
WORKDIR /opt2

# Installs miRDeep2/0.1.3-c6440e2, for more info:
# https://github.com/rajewsky-lab/mirdeep2/tree/c6440e298795579ad62351bba7aff9cc43c50c68
RUN mkdir -p /opt2/mirdeep2/0.1.3-c6440e2/ \
    && git clone https://github.com/rajewsky-lab/mirdeep2.git /opt2/mirdeep2/0.1.3-c6440e2/ \
    && cd /opt2/mirdeep2/0.1.3-c6440e2/ \
    && perl install.pl \
    && touch /opt2/mirdeep2/0.1.3-c6440e2/install_successful


################ POST #################
# Add Dockerfile and export environment 
# variables and update permissions
ADD Dockerfile /opt2/mir-seek.dockerfile
RUN chmod -R a+rX /opt2
ENV PATH="/opt2:$PATH"
WORKDIR /data2
