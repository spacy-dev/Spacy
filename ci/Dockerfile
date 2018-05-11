FROM ubuntu:16.04
MAINTAINER Lars Lubkoll "lars.lubkoll@posteo.de"

RUN apt-get update
RUN apt-get install -y cmake lcov doxygen graphviz git rbenv
RUN gem install coveralls-lcov

WORKDIR /home
