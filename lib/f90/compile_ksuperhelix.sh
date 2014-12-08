#!/bin/bash

#Set the location of the lawson library file
gfortran -c liblawson.f90
ar rcv liblawson.a liblawson.o
ranlib liblawson.a
rm liblawson.o

#Compile once the libraries have been made, where liblawson.a is in the pr-wd
gfortran $1 -L. -llawson -O3 -o eli55a.o
