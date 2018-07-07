#! /bin/bash

BASE_DIR=$1
SAMPLE_NAME=$2
INPUT_MZID=$BASE_DIR/$SAMPLE_NAME.mzid
OUTPUT_TSV=$BASE_DIR/$SAMPLE_NAME.tsv
MSGFPLUS_PATH="/z/home/aoj/opt/MSGFPlus/MSGFPlus.jar"

# Load MSGFPLUS_PATH 
# SOFT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# source $SOFT_DIR/load_flags.sh $1 $2 > /dev/null 2>&1

#if [ ! -f "$OUTPUT_TSV" ] && [ -f "$INPUT_MZID" ]
#then
  java -Xmx3500M -cp $MSGFPLUS_PATH edu.ucsd.msjava.ui.MzIDToTsv -i $INPUT_MZID -o $OUTPUT_TSV
#else
#  echo "MZ2ID not running for sample $SAMPLE_NAME"
#fi
