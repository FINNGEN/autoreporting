FROM ubuntu:22.04

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

ARG HTSLIB_VER=1.21
ARG BCFTOOLS_VER=1.21
ENV HTSLIB_CONFIGURE_OPTIONS="--enable-gcs"
ENV HTSLIB_LIBRARY_DIR=/usr/local/lib
ENV HTSLIB_INCLUDE_DIR=/usr/local/include

RUN apt-get update && apt-get install -qqy wget unzip bzip2 curl libbz2-dev liblzma-dev python3 clang python3-pip zlib1g-dev gnupg2 lsb-release --yes \
&& apt-get install curl tabix --yes && \
    apt-get clean

RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_latest.zip && unzip plink_linux_x86_64_latest.zip  && cp plink /bin  
COPY requirements.txt /usr/local/requirements.txt

#install google cloud sdk
ENV CLOUDSDK_INSTALL_DIR=/usr/local/gcloud/
RUN curl -sSL https://sdk.cloud.google.com | bash
ENV PATH=$PATH:/usr/local/gcloud/google-cloud-sdk/bin


# install HTSLIB
RUN curl -LO https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VER}/bcftools-${BCFTOOLS_VER}.tar.bz2 && \
    tar -xvjf bcftools-${BCFTOOLS_VER}.tar.bz2 && cd bcftools-${BCFTOOLS_VER} && \
    ./configure && make && make install && cd .. && rm -rf bcftools-*


RUN pip3 install  -r usr/local/requirements.txt
ADD Scripts /usr/local/bin
RUN chmod +x /usr/local/bin/main.py
COPY Scripts/wdl_processing_scripts/process_serial.py /usr/local/bin
RUN chmod +x /usr/local/bin/process_serial.py
RUN chmod +x /usr/local/bin/meta_filter_top.py /usr/local/bin/post_process_hits.py 
