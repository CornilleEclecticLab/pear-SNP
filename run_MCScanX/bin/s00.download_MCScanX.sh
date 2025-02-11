#!/usr/bin/bash


cd ~/software

wget https://github.com/wyp1125/MCScanX/archive/refs/tags/v1.0.0.zip

mv v1.0.0.zip MCScanX_v1.0.0.zip

unzip MCScanX_v1.0.0.zip

cd MCScanX-1.0.0

module load java-jdk

make
