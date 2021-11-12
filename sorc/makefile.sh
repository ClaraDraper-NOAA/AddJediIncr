#!/bin/ksh
set -x
#-----------------------------------------------------
#-use standard module.
#-----------------------------------------------------

##export DEBUG='-ftrapuv -check all -check nooutput_conversion -fp-stack-check -fstack-protector -traceback -g'
export INCS="-I$IP_INCd ${NETCDF_INCLUDE}"
export FFLAGS="$INCS -O3 -fp-model precise -r8 -convert big_endian -traceback -g"
export OMPFLAG=-qopenmp
export LDFLG=-qopenmp

export LIBSM="${NETCDF_LDFLAGS_F}"

make -f Makefile clean
make -f Makefile
make -f Makefile install
make -f Makefile clean
