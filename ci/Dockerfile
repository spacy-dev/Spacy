FROM ubuntu:16.04
MAINTAINER Lars Lubkoll "lars.lubkoll@posteo.de"

RUN apt-get update
RUN apt-get install -y software-properties-common
RUN add-apt-repository ppa:ubuntu-toolchain-r/test
RUN apt-get update
RUN apt-get install -y cmake lcov doxygen graphviz git rbenv gcc-6 g++-6 gcc-5 g++-5
RUN gem install coveralls-lcov

WORKDIR /home
