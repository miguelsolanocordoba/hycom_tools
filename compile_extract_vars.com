ftn -o extract_vars_TS_c004_v12.x extract_vars_TS_c004_v12.f -traceback -g -O3 -fp-model source -warn nogeneral -convert big_endian -assume byterecl -mcmodel=small -fPIC
