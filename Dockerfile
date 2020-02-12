FROM ubuntu:19.04 as shorah-builder
ENV LANG=C.UTF-8
RUN apt-get update -qq && \
    apt-get install -qqy python3 python3-pip python3-dev build-essential curl
RUN apt-get install -qqy perl pkg-config zlib1g zlib1g-dev libgsl-dev libhts-dev libboost-dev
RUN pip3 install virtualenv
RUN virtualenv /opt/shorah
RUN /opt/shorah/bin/pip install Biopython numpy
RUN curl -sSLO https://github.com/cbg-ethz/shorah/releases/download/v1.9.95/shorah-1.9.95.tar.bz2 && \
    tar xf shorah-1.9.95.tar.bz2 && \
    cd shorah-1.9.95 && \
    ./configure --prefix=/opt/shorah PYTHON=/opt/shorah/bin/python3.7 && \
    make -j4 && \
    make install

FROM ubuntu:19.04
ENV LANG=C.UTF-8
RUN apt-get update -qq && \
    apt-get install -qqy perl bowtie2 bwa curl python3 openjdk-11-jre \
    zlib1g libgsl23 libhts2 libboost-dev samtools bcftools pigz ncbi-blast+
RUN cd /tmp && \
    curl -sSLO https://bioinformatics.cvr.ac.uk/wp-content/uploads/2019/11/Tanoti-1.2-Linux.tar.gz && \
    tar xf Tanoti-1.2-Linux.tar.gz && \
    mv /tmp/Tanoti-1.2-Linux/* /usr/bin && \
    rm -rf /tmp/Tanoti-1.2*
COPY --from=shorah-builder /opt/shorah /opt/shorah
COPY entrypoints/sam2bam entrypoints/tanoti-wrapper entrypoints/fastq2fasta entrypoints/blastfilter entrypoints/blastfilter_fasta /usr/bin/
RUN ln -s /opt/shorah/bin/shorah /usr/bin/shorah && \
    ln -s ../../../../bin /opt/shorah/lib/python3.7/site-packages/shorah && \
    chmod +x /usr/bin/sam2bam /usr/bin/tanoti-wrapper /usr/bin/fastq2fasta /usr/bin/blastfilter /usr/bin/blastfilter_fasta
RUN mkdir -p /opt/indelfixer && cd /opt/indelfixer && \
    curl -sSLO https://github.com/cbg-ethz/InDelFixer/releases/download/v1.1/InDelFixer.jar
RUN cd /tmp && \
    curl -sSLO http://snap.cs.berkeley.edu/downloads/snap-beta.18-linux.tar.gz && \
    tar xf snap-beta.18-linux.tar.gz && \
    mv /tmp/snap-aligner /usr/bin/snap-aligner && \
    chmod +x /usr/bin/snap-aligner && \
    rm -rf /tmp/snap-beta.18-linux.tar.gz
WORKDIR /workspace
VOLUME /workspace
