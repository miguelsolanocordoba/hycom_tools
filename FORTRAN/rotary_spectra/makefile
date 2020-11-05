## Makefile: rotary_spectra_v01.f rotary_spectr_v02.f
#
# NOTE: This makefile works on NAVY DSRC HPC: 
#     - (KOEHR) koehr.navydsrc.hpc.mil
#     - (GAFFNEY) gaffney.navydsrc.hpc.mil
#
# Before compiling, make sure the following modules are loaded:
#     - costinit 
#     - fftw3
#
# On koehr/gaffney, the default include and library paths provided 
# by module costinit cannot be linked statically. Compilation does
# not work unless the fftw3 libraries are linked dynamically. 
#
# >> ifort $(CPPFLAGS) $(LDFLAGS) -lfftw3 -o $(EXECNAME) $(FROMGRAM) 
#
# Make sure the include directory (CPPFLAGS) and the lib directory
# $(LDFLAGS) are correct, or override it by uncommenting these 
# variables below. 
#
# ----------------------------------------------------------- # 
# EXAMPLE: 
# Compiling rotary_spectra_v02.f 
#
# >> module add costinit 
# >> module add fftw3
# >> make rotv2
#
# Clean previous compilations
#
# >> make clean
#
# ----------------------------------------------------------- # 
# 
## Created: M. Solano, July 9, 2020

# Options
FC = ifort                        # (F)ortran (C)ompiler 
PROGRAM1 = rotary_spectra_v01.f   # source code v01
PROGRAM2 = rotary_spectra_v02.f   # source code v02
EXECNAME1 = rotary_spectra_v01.x  # executable v01
EXECNAME2 = rotary_spectra_v02.x  # executable v02

# lib and include directories for fftw3 (dynamic linking)
#CPPFLAGS = -I/p/app/COST/fftw3/3.3.5/intel/include
#LDFLAGS = -L/p/app/COST/fftw3/3.3.5/intel/lib 

# Compile flags (Intel) !DO NOT CHANGE! 
FCFLAGS = -traceback -xHost -O3 -fp-model precise -ftz -align array64byte -convert big_endian -assume byterecl -warn nogeneral -diag-disable 10212

# Directives
rotv1: $(PROGRAM1)
	$(FC) $(CPPFLAGS) $(LDFLAGS) -lfftw3 -o $(EXECNAME1) $(PROGRAM1) $(FCFLAGS)
rotv2: $(PROGRAM2)
	$(FC) $(CPPFLAGS) $(LDFLAGS) -lfftw3 -o $(EXECNAME2) $(PROGRAM2) $(FCFLAGS)
all: rotv1 rotv2
clean: 
	rm $(EXECNAME1)
	rm $(EXECNAME2)
