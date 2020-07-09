## Makefile: extract_vars.f 
# NOTE: This makefile works on DSRC HPC which uses 'ftn'. 
# 'ftn' is a wrapper for all fortran compilers, make sure 
# module 'PrgEnv-intel' is loaded (>> module list). 
# Type 'ftn -V' for more info. 
# 
# Created: M. Solano, July 9, 2020

# Options
FC = ftn                   # (F)ortran (C)ompiler 
PROGRAM = extract_vars.f   # source code 
EXECNAME = extract_vars.x  # executable

# Compile flags (Intel)  
FCFLAGS = -traceback -g -O3 -fp-model source -warn nogeneral -convert big_endian -assume byterecl -mcmodel=small -fPIC

all: 
	$(FC) -o $(EXECNAME) $(PROGRAM) $(FCFLAGS)
clean: 
	rm $(EXECNAME) 


#ftn -o extract_vars_TS_c004_v12.x extract_vars_TS_c004_v12.f -traceback -g -O3 -fp-model source -warn nogeneral -convert big_endian -assume byterecl -mcmodel=small -fPIC
