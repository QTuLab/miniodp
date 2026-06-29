FROM ubuntu:22.04

ARG VCS_REF=unknown
ARG BUILD_DATE=unknown

LABEL org.opencontainers.image.revision=$VCS_REF \
      org.opencontainers.image.created=$BUILD_DATE

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    curl \
    default-jre \
    gawk \
    fastqc \
    gffread \
    git \
    graphviz \
    hisat2 \
    bowtie2 \
    stringtie \
    samtools \
    bedtools \
    trim-galore \
    picard-tools \
    sra-toolkit \
    python3 \
    python3-pip \
    python-is-python3 \
    python3-dev \
    python3-venv \
    python3-deeptools \
    pkg-config \
    unzip \
    wget \
    zlib1g-dev \
    libbz2-dev \
    libcurl4-openssl-dev \
    libgraphviz-dev \
    liblzma-dev \
    libssl-dev \
    tabix \
    r-base \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir \
      chardet \
      "numpy<2" \
      "distributed==2023.9.1" \
      "MACS2==2.2.9.1" \
      "rseqc==5.0.4" \
      "git+https://github.com/vanheeringen-lab/ANANSE@v0.5.1"

RUN wget -q --output-document /tmp/Bismark.zip https://github.com/FelixKrueger/Bismark/archive/refs/tags/0.24.1.zip && \
    unzip -q /tmp/Bismark.zip -d /usr/local/share/ && \
    find /usr/local/share/Bismark-0.24.1/ -type f -executable -exec ln -s {} /usr/local/bin/ \; && \
    rm -f /tmp/Bismark.zip

RUN wget -q --output-document /usr/local/bin/bedGraphToBigWig https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig && \
    wget -q --output-document /usr/local/bin/bigWigAverageOverBed https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed && \
    wget -q --output-document /usr/local/bin/bedClip https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedClip && \
    chmod +x /usr/local/bin/bedGraphToBigWig /usr/local/bin/bigWigAverageOverBed /usr/local/bin/bedClip

RUN wget -q --output-document /tmp/ROSE-1.3.2.zip https://github.com/stjude/ROSE/archive/refs/tags/1.3.2.zip && \
    unzip -q /tmp/ROSE-1.3.2.zip -d /usr/local/share/ && \
    ln -s /usr/local/share/ROSE-1.3.2/lib/* /usr/local/bin/ && \
    ln -s /usr/local/share/ROSE-1.3.2/bin/* /usr/local/bin/ && \
    rm -f /tmp/ROSE-1.3.2.zip

RUN mkdir -p /work
WORKDIR /work

CMD ["/bin/bash"]
