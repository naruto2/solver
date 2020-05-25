#!/bin/sh

cd ../test
for i in *.mtx; do ../src/solver $i; done
