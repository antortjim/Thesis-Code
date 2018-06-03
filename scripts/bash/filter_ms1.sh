#! /bin/bash

# Absolute path this script is in, thus /home/user/bin
SOFT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SOFT_DIR/load_flags.sh > /dev/null 2>&1

MGF_FILE=$1
FILENAME=$(basename $MGF_FILE)
SAMPLE_NAME="${FILENAME%.*}"
#ROOT_DIR=$(dirname $0)
METADATA_DIR=$ROOT_DIR/$EXP_NAME/data/ms1/
METADATA_FILE=$METADATA_DIR/${SAMPLE_NAME}_ms1.txt
METADATA_TABLE=$METADATA_DIR/${SAMPLE_NAME}_ms1_table.tsv
  
  
echo mgf_file:
echo $MGF_FILE
echo filename:
echo $FILENAME
echo sample_name:
echo $SAMPLE_NAME
echo root_dir:
echo $ROOT_DIR

METADATA_FILE=$ROOT_DIR/$EXP_NAME/data/ms1/${SAMPLE_NAME}_metadata.txt
grep [A-Z] $MGF_FILE > $METADATA_FILE
grep SCANS   $METADATA_FILE  | cut -f2 -d= > $METADATA_DIR/${SAMPLE_NAME}_scans.txt
grep PEPMASS $METADATA_FILE  | cut -f2 -d=  | cut -f1 -d' ' > $METADATA_DIR/${SAMPLE_NAME}_pepmass_mz.txt
grep PEPMASS $METADATA_FILE  | cut -f2 -d=  | cut -f2 -d' ' > $METADATA_DIR/${SAMPLE_NAME}_pepmass_intensity.txt
grep RTINSECONDS  $METADATA_DIR/${SAMPLE_NAME}_metadata.txt | cut -f2 -d= > $METADATA_DIR/${SAMPLE_NAME}_rt.txt

echo -e "Spectrum.Scan.Number\tpepmass_mz\tpepmass_intensity\trt" > $METADATA_TABLE
paste $METADATA_DIR/${SAMPLE_NAME}_scans.txt $METADATA_DIR/${SAMPLE_NAME}_pepmass_mz.txt $METADATA_DIR/${SAMPLE_NAME}_pepmass_intensity.txt  $METADATA_DIR/${SAMPLE_NAME}_rt.txt >> $METADATA_TABLE
rm $METADATA_DIR/${SAMPLE_NAME}_scans.txt $METADATA_DIR/${SAMPLE_NAME}_pepmass_mz.txt $METADATA_DIR/${SAMPLE_NAME}_pepmass_intensity.txt  $METADATA_DIR/${SAMPLE_NAME}_rt.txt $METADATA_DIR/${SAMPLE_NAME}_metadata.txt
