#! /bin/bash

SOFT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SOFT_DIR/load_flags.sh > /dev/null 2>&1


#########################################
## Data processing
#########################################

EXPERIMENTAL_DESIGN="$ROOT_DIR/$EXP_NAME/$SPECTRA/../experimental_design.tsv"

i=1
for FILEPATH in $ROOT_DIR/$EXP_NAME/$SPECTRA/*.mgf
do
  echo $FILEPATH
  SAMPLE_ID="sample_$i"
  FILENAME=$(basename $FILEPATH)
  echo $FILENAME
  SAMPLE_NAME="${FILENAME%.*}"
  echo $SAMPLE_NAME

  if [ -f $EXPERIMENTAL_DESIGN ]
  then
    CONDITION_NAME=$(grep $SAMPLE_NAME $EXPERIMENTAL_DESIGN | awk '{print $5}')
    REPLICATE=$(grep $SAMPLE_NAME $EXPERIMENTAL_DESIGN | awk '{print $1}')
  else
    CONDITION_NAME="condition_$i"
    REPLICATE="$i"
  fi
  
  echo $CONDITION_NAME
  echo $REPLICATE

  #IDENTIFICATION_FILES=$(find `pwd`/searchgui_out -maxdepth 1 | grep "$SAMPLE_NAME" | grep -v -e "\.mgf$" )
  #IDENTIFICATION_FILES_STRING="$(echo $IDENTIFICATION_FILES | sed 's/ /, /g')"
  MZID=$ROOT_DIR/$EXP_NAME/searchgui_out/$SAMPLE_NAME.msgf.mzid

  if [ ! -f "$PS_OUT/$SAMPLE_NAME.cpsx" ] && [ -f $MZID ]
  then
    echo "Calling peptideshaker"

    idconvert $MZID --pepXML -v -o $ROOT_DIR/$EXP_NAME/searchgui_out/
    sed 's/|\(\w*\)_REVERSED|/|REVERSED_\1|/' $ROOT_DIR/$EXP_NAME/searchgui_out/$SAMPLE_NAME.pepXML > $ROOT_DIR/$EXP_NAME/searchgui_out/${SAMPLE_NAME}_prefix.pepXML

    echo $ROOT_DIR/$EXP_NAME/searchgui_out/identification_files_$SAMPLE_NAME.txt

    find $ROOT_DIR/$EXP_NAME/searchgui_out -maxdepth 1 | grep "$SAMPLE_NAME" | grep -v -e "\.mgf$"  > $ROOT_DIR/$EXP_NAME/searchgui_out/identification_files_$SAMPLE_NAME.txt
    zip -j $ROOT_DIR/$EXP_NAME/searchgui_out/identification_files_$SAMPLE_NAME.zip -@ < $ROOT_DIR/$EXP_NAME/searchgui_out/identification_files_$SAMPLE_NAME.txt

    echo $FILEPATH

#    java -cp $PEPTIDESHAKER_PATH eu.isas.peptideshaker.cmd.PeptideShakerCLI \
#      -experiment $EXP_NAME \
#      -sample $CONDITION_NAME \
#      -replicate $REPLICATE  \
#      -identification_files $ROOT_DIR/$EXP_NAME/searchgui_out/identification_files_$SAMPLE_NAME.zip \
#      -spectrum_files \"$FILEPATH\" \
#      -id_params $ROOT_DIR/$EXP_NAME/$SETTINGS_DIR/$PARAMS_NAME.par \
#      -out $PS_OUT/$SAMPLE_NAME.cpsx

  fi
  if [ -f "$PS_OUT/$SAMPLE_NAME.cpsx" ]
  then
     # #########################################
     # ## Report 
     # #########################################
      
     java -cp $PEPTIDESHAKER_PATH eu.isas.peptideshaker.cmd.ReportCLI -in $PS_OUT/$SAMPLE_NAME.cpsx \
       -out_reports $ROOT_DIR/$EXP_NAME/$PS_OUT/reports -reports "0,1,2,3,4,5,6,7,8"
     
     # #########################################
     # ## Clean reports for easy R parsing
     # #########################################
      
     cd $ROOT_DIR/$EXP_NAME/$PS_OUT/reports
     for FILE in *.txt
     do
       cat $FILE | sed 's/#/Number_/g' | sed "s/'/prime/g" > ${FILE%.*}_clean.txt 
       #tail -n +2 $FILE  | sed "s/'/prime/g" >> ${FILE%.*}_clean.txt 
     done
     cd -

  #rm $FILEPATH
  fi

 ((i++))
done

