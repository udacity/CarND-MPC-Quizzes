#! /bin/bash

# update
sudo apt-get update

# gfortran dependency

# get unzip
sudo apt-get install unzip

# Ipopt: get, install, unzip
sudo apt-get install gfortran
wget https://www.coin-or.org/download/source/Ipopt/Ipopt-3.12.7.zip && unzip Ipopt-3.12.7.zip && rm Ipopt-3.12.7.zip
./install_ipopt.sh

# CppAD
sudo apt-get install cppad

# python and matplotlib
sudo apt-get install python-matplotlib
sudo apt-get install python2.7-dev
