ARG UBUNTU_VERSION=20.04

FROM ubuntu:${UBUNTU_VERSION} AS base

ARG DEBIAN_FRONTEND=noninteractive

ENV TZ=America/Vancouver \
    RUBY_VERSION=2.2.2 \
    BUNDLER_VERSION=1.17.3

RUN apt-get -y update &&\
    apt-get -y upgrade &&\
    apt-get install -y curl gpg

RUN echo "deb http://security.ubuntu.com/ubuntu/ bionic-security main" >> /etc/apt/sources.list &&\
    echo "deb http://security.ubuntu.com/ubuntu/ bionic main" >> /etc/apt/sources.list &&\
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 3B4FE6ACC0B21F32 &&\
    apt update &&\
    apt install -y \
    gawk=1:4.1.4+dfsg-1build1

RUN curl -sSL https://rvm.io/mpapis.asc | gpg --import - && \
    curl -sSL https://rvm.io/pkuczynski.asc | gpg --import - && \
    curl -sSL https://get.rvm.io | bash -s stable && \
    usermod -a -G rvm root

RUN /bin/bash -l -c "source /etc/profile.d/rvm.sh && rvm pkg install openssl"
RUN /bin/bash -l -c "source /etc/profile.d/rvm.sh && rvm requirements && rvm install ${RUBY_VERSION} --with-openssl-dir=/usr/local/rvm/usr && rvm --default use ${RUBY_VERSION}"
