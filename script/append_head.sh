#!/bin/sh

src=$1
echo $src
shift
for arg in "$@" 
do
  cp $src $src.temp
  cat $arg >>$src.temp
  cat $src.temp > $arg
  rm $src.temp
  echo $arg
done
