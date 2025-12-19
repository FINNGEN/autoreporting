FROM ubuntu:22.04

ENV HTSLIB_CONFIGURE_OPTIONS="--enable-gcs"

# Install gcsfuse

RUN apt-get update && apt-get install -qqy wget unzip bzip2 curl libbz2-dev liblzma-dev python3 clang python3-pip zlib1g-dev gnupg2 lsb-release --yes \
&& apt-get install curl tabix --yes && \
    apt-get clean
#RUN export GCSFUSE_REPO=gcsfuse-`lsb_release -c -s` && \
#    echo "deb http://packages.cloud.google.com/apt $GCSFUSE_REPO main" | tee /etc/apt/sources.list.d/gcsfuse.list && \
#    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
#    echo "deb https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
#    apt-get update --fix-missing && \
#    apt-get install -y gcsfuse google-cloud-cli && \
#    apt-get clean && \
#    rm -rf /var/lib/apt/lists/*

RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_latest.zip && unzip plink_linux_x86_64_latest.zip  && cp plink /bin  
COPY requirements.txt /usr/local/requirements.txt
RUN pip3 install --no-binary -r usr/local/requirements.txt
ADD Scripts /usr/local/bin
RUN chmod +x /usr/local/bin/main.py
COPY Scripts/wdl_processing_scripts/process_serial.py /usr/local/bin
RUN chmod +x /usr/local/bin/process_serial.py
RUN chmod +x /usr/local/bin/meta_filter_top.py /usr/local/bin/post_process_hits.py 
