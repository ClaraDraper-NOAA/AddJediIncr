#! /usr/bin/env bash
set -eux

source hera_modules

export FCMP=mpiifort

# Check final exec folder exists
if [ ! -d "./exec" ]; then
  mkdir ./exec
fi

cd ./sorc/
./makefile.sh
