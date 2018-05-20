#! /bin/bash

BASE_DIR=$1
SAMPLE_NAME=$2
INPUT_MZID=$BASE_DIR/$SAMPLE_NAME.msgf.mzid
OUTPUT_TSV=$BASE_DIR/$SAMPLE_NAME.tsv
if [ ! -f "$OUTPUT_TSV" ] && [ -f "$INPUT_MZID" ]
then
  java -Xmx3500M -cp ~/opt/MSGFPlus/MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i $INPUT_MZID -o $OUTPUT_TSV
else
  echo "MZ2ID not running for sample $SAMPLE_NAME"
fi
