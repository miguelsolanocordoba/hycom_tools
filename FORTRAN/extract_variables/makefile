## Makefile: extract_vars2D.f extract_grid.f 
# NOTE: This makefile works on DSRC HPC which uses 'ftn'. 
# 'ftn' is a wrapper for all fortran compilers, make sure 
# module 'PrgEnv-intel' is loaded (>> module list). 
# On Cray computers 'PrgEnv-cray' is the default. If so, you
# can change modules typing:
#
# >> module swap PrgEnv-cray PrgEnv-intel
#
# Type 'ftn -V' for more info. 
# 
# Created: M. Solano, July 9, 2020

# Options
FC = ftn                   # (F)ortran (C)ompiler 
PROGRAM = extract_vars2D.f   # source code 
PROGRAM2 = extract_grid.f   # source code 
EXECNAME = extract_vars2D.x  # executable
EXECNAME2 = extract_grid.x  # executable

# Compile flags (Intel)  
FCFLAGS = -traceback -g -O3 -fp-model source -warn nogeneral -convert big_endian -assume byterecl -mcmodel=small -fPIC

grid: 
	$(FC) -o $(EXECNAME2) $(PROGRAM2) $(FCFLAGS)
vars: 
	$(FC) -o $(EXECNAME) $(PROGRAM) $(FCFLAGS)
all: 
	$(FC) -o $(EXECNAME) $(PROGRAM) $(FCFLAGS)
	$(FC) -o $(EXECNAME2) $(PROGRAM2) $(FCFLAGS)
clean: 
	rm $(EXECNAME) $(EXECNAME2) 


#ftn -o extract_vars_TS_c004_v12.x extract_vars_TS_c004_v12.f -traceback -g -O3 -fp-model source -warn nogeneral -convert big_endian -assume byterecl -mcmodel=small -fPIC
