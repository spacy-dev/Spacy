FROM ubuntu:16.04
MAINTAINER Lars Lubkoll "lars.lubkoll@posteo.de"

RUN apt-get update
RUN apt-get install -y software-properties-common
RUN add-apt-repository ppa:ubuntu-toolchain-r/test
RUN apt-get update
RUN apt-get install -y cmake lcov doxygen graphviz git rbenv gcc-6 g++-6 gcc-5 g++-5 libboost-dev
RUN gem install coveralls-lcov
RUN cd /home && git clone https://github.com/eigenteam/eigen-git-mirror.git && cd cd eigen-git-mirror && mkdir build && cd build && cmake .. && sudo make install

WORKDIR /home
