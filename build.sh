#! /usr/bin/env bash
set -eux

#source ./machine-setup.sh > /dev/null 2>&1

if ( ! eval module help > /dev/null 2>&1 ) ; then
    echo load the module command 1>&2
    source /apps/lmod/lmod/init/$__ms_shell
fi
target=hera
module purge
MOD_PATH=/scratch2/NCEPDEV/nwprod/NCEPLIBS/modulefiles

cwd=`pwd`

export MOD_PATH
source ./global_cycle.$target             > /dev/null 2>&1

# Check final exec folder exists
if [ ! -d "./exec" ]; then
  mkdir ./exec
fi

cd ${cwd}/sorc/
./makefile.sh
