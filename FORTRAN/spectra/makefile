#!/bin/bash 
## Makefile: fft_uvpModes.f rotary_spectra_uvModes.f KE_freqband_uvModes.f
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
## Created: M. Solano, February 24, 2021

# Options
FC = ifort                           # (F)ortran (C)ompiler 

# Programs
PROGRAM1 = fft_uvpModes.f             # fft of u/v/p modal amplitudes 
PROGRAM2 = rotary_spectra_uvModes.f   # rotary spectra u/v modal amplitudes
PROGRAM3 = KE_freqband_uvModes.f      # kinetic energy over frequency band

# Executables (*need to be soft linked in running dir*)
EXECNAME1 = fft_uvpModes.x            # fft of u/v/p modal amplitudes
EXECNAME2 = rotary_spectra_uvModes.x  # rotary spectra u/v modal amplitudes
EXECNAME3 = KE_freqband_uvModes.x     # kinetic energy over frequency band

# lib and include directories for fftw3 (dynamic linking - optional)
#CPPFLAGS = -I/p/app/COST/fftw3/3.3.5/intel/include
#LDFLAGS = -L/p/app/COST/fftw3/3.3.5/intel/lib 

# Compile flags (Intel) !DO NOT CHANGE! 
FCFLAGS = -traceback -xHost -O3 -fp-model precise -ftz -align array64byte -convert big_endian -assume byterecl -warn nogeneral -diag-disable 10212

# Directives
fft: $(PROGRAM1)
	$(FC) -o $(EXECNAME1) $(PROGRAM1) $(FCFLAGS)
#	$(FC) $(CPPFLAGS) $(LDFLAGS) -lfftw3 -o $(EXECNAME1) $(PROGRAM1) $(FCFLAGS)
rotspec: $(PROGRAM2)
	$(FC) -o $(EXECNAME2) $(PROGRAM2) $(FCFLAGS)
#	$(FC) $(CPPFLAGS) $(LDFLAGS) -lfftw3 -o $(EXECNAME2) $(PROGRAM2) $(FCFLAGS)
ke: $(PROGRAM3)
	$(FC) -o $(EXECNAME3) $(PROGRAM3) $(FCFLAGS)
#	$(FC) $(CPPFLAGS) $(LDFLAGS) -lfftw3 -o $(EXECNAME3) $(PROGRAM3) $(FCFLAGS)
all: fft rotspec ke
clean: 
	rm $(EXECNAME1)
	rm $(EXECNAME2)
	rm $(EXECNAME3)
