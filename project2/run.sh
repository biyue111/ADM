#!/bin/sh
# run findprime

mpirun -n $1 ./findprime $2 $3
