FROM ubuntu:14.04

RUN apt update && apt install -y autoconf libtool python libgsl0-dev libpython-dev libboost-dev g++ python-numpy python-scipy make libboost-python-dev libboost-regex-dev python-h5py

RUN mkdir -p /srv/
ADD . /srv/epdp
WORKDIR /srv/epdp/
RUN ./autogen.sh && ./configure && make && make install
