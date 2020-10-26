# HYCOM_TOOLS Auxiliary code for HYCOM
#
# Created: July 8, 2020 by M. Solano

This repository contains code to manage HYCOM input/output, 
and provides code to read binary files (*.a and *.BinF)

The repository is divided into programming languages: 

FORTRAN 
extract_vars.f - Main fortran file to extract HYCOM variables

BASH 
make_lis.com - bash script to create list of files to read
compile_extract_vars.com - makefile for extract_vars.f 
run_extract_vars.com - run submit script for DSRC HPCs. 

MATLAB

