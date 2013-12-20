#!/bin/bash

for ((i=1;i<=31;i++)) ; do
  for ((j=1;j<=i;j++)) ; do
    echo "------------------------------"
    echo "downloading atom/ion:"$i","$j
    echo "------------------------------"
    ./download-one-atom-ion.sh $i $j
  done
done