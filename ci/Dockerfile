FROM ubuntu:20.04
MAINTAINER Lars Lubkoll "lars.lubkoll@posteo.de"

RUN apt-get update
RUN apt-get install -y cmake lcov doxygen graphviz git rbenv gcc-6 g++-6 gcc-5 g++-5 libboost-dev g++-10 gcc-10 clang-10 libeigen3-dev liblapack-dev
RUN gem install coveralls-lcov

WORKDIR /home
