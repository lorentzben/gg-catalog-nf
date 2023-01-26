FROM ubuntu:20.04
WORKDIR /filtlong_build

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update && apt-get -y upgrade && \
    apt install -y procps && \
    apt install -y python3 && \
    apt-get install -y python3-pandas && \
	apt install -y wget unzip make g++ cmake git && \
    apt install -y libz-dev
    
    
RUN cd /opt && git clone https://github.com/lh3/minimap2 && \
    cd minimap2 && make &&\
    cp minimap2 /usr/local/bin && \
    chmod +x /usr/local/bin/minimap2 && chmod 755 /usr/local/bin/minimap2 

ENV PATH="/usr/local/bin/minimap2:${PATH}"