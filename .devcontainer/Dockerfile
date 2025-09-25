FROM mcr.microsoft.com/devcontainers/python:1-3.10

ARG RMLMAPPER_JAR_URL="https://github.com/RMLio/rmlmapper-java/releases/download/v7.3.3/rmlmapper-7.3.3-r374-all.jar"

USER root 
RUN apt-get update && apt-get -y install --no-install-recommends \
        ca-certificates \
        libcurl4-openssl-dev \
        libglpk40 \
        libglpk-dev \
        poppler-utils \
        swi-prolog \
        minizinc && \
update-ca-certificates
RUN mkdir -p /opt/tools && cd /opt/tools && curl -LO $RMLMAPPER_JAR_URL

RUN curl -L https://rig.r-pkg.org/deb/rig.gpg -o /etc/apt/trusted.gpg.d/rig.gpg && \
        sh -c 'echo "deb http://rig.r-pkg.org/deb rig main" > /etc/apt/sources.list.d/rig.list' && \
        apt update && \
        apt install r-rig && \
        rig add 4.2.3 && \
        rig default 4.2.3

ARG MINIZINC_TAR_URL="https://github.com/MiniZinc/MiniZincIDE/releases/download/2.8.7/MiniZincIDE-2.8.7-bundle-linux-x86_64.tgz"
RUN cd /opt/tools && curl -LO $MINIZINC_TAR_URL && \
    tar xf MiniZincIDE-2.8.7-bundle-linux-x86_64.tgz && \
    ln -s /opt/tools/MiniZincIDE-2.8.7-bundle-linux-x86_64/bin/minizinc /usr/local/bin/minizinc

