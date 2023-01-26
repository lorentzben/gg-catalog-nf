FROM staphb/filtlong:0.2.1
WORKDIR /filtlong_build

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update && apt-get -y upgrade && \
    apt install -y procps && \
    apt install -y python3 && \
    apt-get install -y python3-pandas && \
	apt install -y wget unzip make g++ cmake git libc-bin && \
    apt install -y libz-dev
    
    
