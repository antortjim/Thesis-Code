#! /bin/bash

DIR=$1
FILES=$DIR/*mgf
OUT_DIR=$2
echo ${FILES[*]}


for MGF_FILE in $FILES
do
  ~/opt/tpp/tpp/bin/msconvert $MGF_FILE --mzML -o $OUT_DIR
done
