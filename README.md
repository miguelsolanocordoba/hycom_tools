# HYCOM_TOOLS Auxiliary code for HYCOM
#
# Created: July 8, 2020 by M. Solano

This repository contains code to manage HYCOM input/output, 
and provides code to read binary files (*.a and *.BinF)

The repository is divided by programming language): 

FORTRAN 
extract_vars.f    - Main fortran file to extract HYCOM variables
extract_vars2D.f  - Extract ATLc0.02 2D variables (8 total) 
extract_grid.f    - Extract ATLc0.02 grid (griddata)
makefile          - custom makefile for all (.f)

BASH 
make_lis.com - bash script to create list of files to read
compile_extract_vars.com - makefile for extract_vars.f 
run_extract_vars.com - run submit script for DSRC HPCs. 

MATLAB
read_hgrida.m - read Atlc0.02 grid
read_hvarsa.m - read ATLc0.02 variables (2D only)


#------------------ EXTRACT_VARIABLES --------------------# 

To extract HYCOM output (*.a) into tiles (*.BinF), you will 
need to use code from BASH and FORTRAN. Below is an example
of how to do this for the ATLc0.02, expt_04.3:

1. Modify make_lis.com to create the list (.lis) file. 

Make sure the variable root (FIN) is correct. 
Modify run number (RNM). 
Modify day/hour start/end (is,ie,js,je). 
Make sure the paths for the HYCOM output is correct, and that 
the appropriate permissions are granted to read/execute. 

2. Run the script file
>>./make_lis.com
Make sure the list created is correct.

2. Modify extract_vars2D.f to choose variable to extract. 
Use the flags at the headers to toggle on/ff (0=off, 1=on)
Make sure the parameters idm/jdm/kdm are correct. 


3. Compile the code.
To compile the code on a Cray computer (e.g. onyx), first start
by swaping from the Cray to the Intel compiler: 

>> module swap PrgEnv-cray PrgEnv-intel

To compile extract_vars2D.f, type:
>> make vars

 
4. Modify the run script run_extract_vars.com to submit job.


5. Run the code.
>> ./run_extractvars.com 
