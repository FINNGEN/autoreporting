FROM ubuntu:18.04
RUN apt-get update && apt-get install -qqy wget unzip bzip2 curl python3.6 python3-venv python3-pip\
&& apt-get install curl make tabix python3 python3-pip zlib1g-dev --yes && \
    apt-get clean && \
    pip3 install pysam requests numpy pandas && \
wget http://www.christianbenner.com/ldstore_v1.1_x86_64.tgz && tar xvf ldstore_v1.1_x86_64.tgz && \
wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_latest.zip && unzip plink_linux_x86_64_latest.zip  && cp plink /bin  
RUN cp ldstore_v1.1_x86_64/ldstore /usr/local/bin/
RUN chmod +x /usr/local/bin/ldstore
ADD Scripts /usr/local/bin
RUN chmod +x /usr/local/bin/main.py
COPY Scripts/wdl_processing_scripts/process_serial.py /usr/local/bin
RUN chmod +x /usr/local/bin/process_serial.py
