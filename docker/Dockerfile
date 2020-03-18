FROM ubuntu:18.04
RUN apt-get update && apt-get install -qqy wget unzip bzip2 curl python3.6 python3-venv python3-pip git nano \
&& apt-get install curl make gcc g++ cmake tabix python3 python3-pip zlib1g-dev --yes && \
    apt-get clean && \
    pip3 install pytabix requests numpy pandas && \
wget http://www.christianbenner.com/ldstore_v1.1_x86_64.tgz && tar xvf ldstore_v1.1_x86_64.tgz && \
wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_latest.zip && unzip plink_linux_x86_64_latest.zip  && cp plink /bin  
RUN cp ldstore_v1.1_x86_64/ldstore /usr/local/bin/
RUN chmod +x /usr/local/bin/ldstore
COPY Scripts/main.py /usr/local/bin/
RUN chmod +x /usr/local/bin/main.py
COPY Scripts/gws_fetch.py /usr/local/bin/
COPY Scripts/annotate.py /usr/local/bin/
COPY Scripts/compare.py /usr/local/bin/
COPY Scripts/autoreporting_utils.py /usr/local/bin/
COPY Scripts/gwcatalog_api.py /usr/local/bin/
COPY Scripts/linkage.py /usr/local/bin/
COPY Scripts/wdl_processing_scripts/process_serial.py /usr/local/bin
RUN chmod +x /usr/local/bin/process_serial.py
COPY Scripts/wdl_processing_scripts/process_wdl_results.py /usr/local/bin
RUN chmod +x /usr/local/bin/process_wdl_results.py
ADD data ../data