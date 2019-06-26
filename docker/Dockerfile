FROM ubuntu:18.04
RUN apt-get update && apt-get install -qqy wget unzip bzip2 curl python3.6 python3-venv python3-pip git nano \
&& apt-get install curl make gcc g++ cmake tabix python3 python3-pip plink1.9 zlib1g-dev --yes && \
    apt-get clean && \
    pip3 install pytabix requests numpy pandas && \
wget http://www.christianbenner.com/ldstore_v1.1_x86_64.tgz && tar xvf ldstore_v1.1_x86_64.tgz  && cp ldstore_v1.1_x86_64/ldstore /usr/local/bin/ldstore
COPY Scripts/main.py /usr/local/bin/
COPY Scripts/gws_fetch.py /usr/local/bin/
COPY Scripts/annotate.py /usr/local/bin/
COPY Scripts/compare.py /usr/local/bin/
COPY Scripts/autoreporting_utils.py /usr/local/bin/
COPY Scripts/gwcatalog_api.py /usr/local/bin/
