#!/bin/bash

for i in tf_2d tf_3d
do
   cd $i
   make cpu
   cd ..
done
