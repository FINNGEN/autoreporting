FROM ubuntu:22.10

RUN apt-get update && apt-get install -qqy wget unzip bzip2 curl libbz2-dev liblzma-dev python3 python3-pip zlib1g-dev --yes \
&& apt-get install curl tabix --yes && \
    apt-get clean
RUN wget http://www.christianbenner.com/ldstore_v1.1_x86_64.tgz && tar xvf ldstore_v1.1_x86_64.tgz && \
    wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_latest.zip && unzip plink_linux_x86_64_latest.zip  && cp plink /bin  
RUN cp ldstore_v1.1_x86_64/ldstore /usr/local/bin/
RUN chmod +x /usr/local/bin/ldstore
COPY requirements.txt /usr/local/requirements.txt
RUN pip3 install -r usr/local/requirements.txt
RUN cp ldstore_v1.1_x86_64/ldstore /usr/local/bin/
RUN chmod +x /usr/local/bin/ldstore
ADD Scripts /usr/local/bin
RUN chmod +x /usr/local/bin/main.py
COPY Scripts/wdl_processing_scripts/process_serial.py /usr/local/bin
RUN chmod +x /usr/local/bin/process_serial.py
RUN chmod +x /usr/local/bin/meta_filter_top.py /usr/local/bin/post_process_hits.py 
