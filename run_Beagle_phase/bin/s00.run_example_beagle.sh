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

echo
echo "*** Creating test files: ref.01Mar24.d36.vcf.gz target.01Mar24.d36.vcf.gz ***"
echo
zcat test.01Mar24.d36.vcf.gz | cut -f1-190 | tr '/' '|' | gzip > ref.01Mar24.d36.vcf.gz
zcat test.01Mar24.d36.vcf.gz | cut -f1-9,191-200 | gzip > target.01Mar24.d36.vcf.gz

echo
echo "*** Running test analysis with \"gt=\" argument ***"
echo
java -jar beagle.01Mar24.d36.jar gt=test.01Mar24.d36.vcf.gz out=out.gt

echo
echo "*** Running test analysis with \"ref=\" and \"gt=\" arguments ***"
echo
java -jar beagle.01Mar24.d36.jar ref=ref.01Mar24.d36.vcf.gz gt=target.01Mar24.d36.vcf.gz out=out.ref

echo
echo "*** Making \"bref3\" file ***"
echo
java -jar bref3.01Mar24.d36.jar ref.01Mar24.d36.vcf.gz > ref.01Mar24.d36.bref3

echo
echo "*** Running test analysis with \"bref3\" file ***"
echo
java -jar beagle.01Mar24.d36.jar ref=ref.01Mar24.d36.bref3 gt=target.01Mar24.d36.vcf.gz out=out.bref3
