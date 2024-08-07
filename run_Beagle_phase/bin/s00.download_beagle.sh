#!/bin/bash
if [ ! -f beagle.01Mar24.d36.jar ]; then
  echo
  echo "Downloading beagle.01Mar24.d36.jar"
  wget http://faculty.washington.edu/browning/beagle/beagle.01Mar24.d36.jar
fi

if [ ! -f bref3.01Mar24.d36.jar ]; then
  echo
  echo "Downloading bref3.01Mar24.d36.jar"
  wget http://faculty.washington.edu/browning/beagle/bref3.01Mar24.d36.jar
fi

echo

if [ ! -f test.01Mar24.d36.vcf.gz ]; then
    echo
    echo "*** Downloading some 1000 Genomes Project data to file: test.01Mar24.d36.vcf.gz ***"
    wget http://faculty.washington.edu/browning/beagle/test.01Mar24.d36.vcf.gz
fi
