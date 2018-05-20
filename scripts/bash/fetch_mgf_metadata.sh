#! /bin/bash

MGF_PATH=$1
SCRIPT=$(readlink -f "$0")
# Absolute path this script is in, thus /home/user/bin
source load_flags.sh > /dev/null 2>&1
echo $MGF_PATH
SEARCH_OUTPUT_DIR=$ROOT_DIR/searchgui_out/

for MGF_FILE in $MGF_PATH/*.mgf; do
  FILENAME=$(basename $MGF_FILE)
  SAMPLE_NAME="${FILENAME%.*}"
  #ROOT_DIR=$(dirname $0)
  METADATA_DIR=$ROOT_DIR/data/mgf_metadata/
  METADATA_TABLE=$METADATA_DIR/${SAMPLE_NAME}_metadata_table.tsv
  SEARCH_OUTPUT=$SEARCH_OUTPUT_DIR/${SAMPLE_NAME}.tsv
  
  
  echo mgf_file:
  echo $MGF_FILE
  echo filename:
  echo $FILENAME
  echo sample_name:
  echo $SAMPLE_NAME
  echo path:
  echo $MGF_PATH
  echo $ROOT_DIR
  echo $SEARCH_OUTPUT
  
  grep [A-Z] $MGF_FILE > $METADATA_DIR/${SAMPLE_NAME}_metadata.txt
  grep SCANS  $METADATA_DIR/${SAMPLE_NAME}_metadata.txt | cut -f2 -d= > $METADATA_DIR/${SAMPLE_NAME}_scans.txt
  grep PEPMASS  $METADATA_DIR/${SAMPLE_NAME}_metadata.txt | cut -f2 -d=  | cut -f1 -d' ' > $METADATA_DIR/${SAMPLE_NAME}_pepmass_mz.txt
  grep PEPMASS  $METADATA_DIR/${SAMPLE_NAME}_metadata.txt | cut -f2 -d=  | cut -f2 -d' ' > $METADATA_DIR/${SAMPLE_NAME}_pepmass_intensity.txt
  grep RTINSECONDS  $METADATA_DIR/${SAMPLE_NAME}_metadata.txt | cut -f2 -d= > $METADATA_DIR/${SAMPLE_NAME}_rt.txt

    
  bash tomzid.sh $SEARCH_OUTPUT_DIR $SAMPLE_NAME
  echo $SEARCH_OUTPUT_DIR
  echo $SAMPLE_NAME
  cut -f 3 -d\t $SEARCH_OUTPUT | tail -n +2 | sort | uniq > $METADATA_DIR/${SAMPLE_NAME}_identified.txt
  
  echo -e "scan\tpepmass_mz\tpepmass_intensity\trt" > $METADATA_TABLE
  paste $METADATA_DIR/${SAMPLE_NAME}_scans.txt $METADATA_DIR/${SAMPLE_NAME}_pepmass_mz.txt $METADATA_DIR/${SAMPLE_NAME}_pepmass_intensity.txt  $METADATA_DIR/${SAMPLE_NAME}_rt.txt >> $METADATA_TABLE
  rm $METADATA_DIR/${SAMPLE_NAME}_scans.txt $METADATA_DIR/${SAMPLE_NAME}_pepmass_mz.txt $METADATA_DIR/${SAMPLE_NAME}_pepmass_intensity.txt  $METADATA_DIR/${SAMPLE_NAME}_rt.txt $MGF_FILE > $METADATA_DIR/${SAMPLE_NAME}_metadata.txt
done
