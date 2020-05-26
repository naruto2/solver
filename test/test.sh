#!/bin/sh

cd ../test

../src/solver *.mtx

exit 

for i in *.mtx; do ../src/solver $i; done
