#! /bin/bash

# This script is related to install-ubuntu_MPC.sh
# workspaces log in as root and calls to sudo
# are not accepted
# matplotlib and python are not required but are 
# left in for future development

# update
apt-get update

# gfortran dependency

# get unzip
apt-get install unzip

# Ipopt: get, install, unzip
apt-get install gfortran
wget https://www.coin-or.org/download/source/Ipopt/Ipopt-3.12.7.zip && unzip Ipopt-3.12.7.zip && rm Ipopt-3.12.7.zip
./install_ipopt.sh

# CppAD
apt-get install cppad

# python and matplotlib
apt-get install python-matplotlib
apt-get install python2.7-dev
