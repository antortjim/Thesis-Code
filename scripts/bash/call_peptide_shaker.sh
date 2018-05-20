#! /bin/bash

SOFT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SOFT_DIR/load_flags.sh > /dev/null 2>&1


#########################################
## Data processing
#########################################

i=1
for FILEPATH in $ROOT_DIR/$EXP_NAME/$SPECTRA/*.mgf
do
  SAMPLE_ID="sample_$i"
  FILENAME=$(basename $FILEPATH)
  SAMPLE_NAME=${FILENAME%.*}
  echo $FILEPATH
  echo $FILENAME
  echo $SAMPLE_NAME

  CONDITION_NAME=$(grep $SAMPLE_NAME  $ROOT_DIR/$EXP_NAME/$SPECTRA/../experimental_design.tsv | awk '{print $5}')
  REPLICATE=$(grep $SAMPLE_NAME  $ROOT_DIR/$EXP_NAME/$SPECTRA/../experimental_design.tsv | awk '{print $1}')
  
  echo $CONDITION_NAME
  echo $REPLICATE

  #IDENTIFICATION_FILES=$(find `pwd`/searchgui_out -maxdepth 1 | grep "$SAMPLE_NAME" | grep -v -e "\.mgf$" )
  #IDENTIFICATION_FILES_STRING="$(echo $IDENTIFICATION_FILES | sed 's/ /, /g')"

  if [ ! -f "$PS_OUT/$SAMPLE_NAME.cpsx" ] && [ -f $ROOT_DIR/$EXP_NAME/searchgui_out/$SAMPLE_NAME.msgf.mzid ]
  then
    echo "Calling peptideshaker"
    find $ROOT_DIR/$EXP_NAME/searchgui_out -maxdepth 1 | grep "$SAMPLE_NAME" | grep -v -e "\.mgf$"  > searchgui_out/identification_files_$SAMPLE_NAME.txt
    zip -j $ROOT_DIR/$EXP_NAME/searchgui_out/identification_files_$SAMPLE_NAME.zip -@ < $ROOT_DIR/$EXP_NAME/searchgui_out/identification_files_$SAMPLE_NAME.txt

    echo $FILEPATH

    java -cp $PEPTIDESHAKER_PATH eu.isas.peptideshaker.cmd.PeptideShakerCLI \
      -experiment $EXP_NAME \
      -sample $CONDITION_NAME \
      -replicate $REPLICATE  \
      -identification_files $ROOT_DIR/$EXP_NAME/searchgui_out/identification_files_$SAMPLE_NAME.zip \
      -spectrum_files \"$FILEPATH\" \
      -id_params $ROOT_DIR/$EXP_NAME/$SETTINGS_DIR/$PARAMS_NAME.par \
      -out $PS_OUT/$SAMPLE_NAME.cpsx

     # #########################################
     # ## Report 
     # #########################################
      
     java -cp $PEPTIDESHAKER_PATH eu.isas.peptideshaker.cmd.ReportCLI -in $PS_OUT/$SAMPLE_NAME.cpsx \
       -out_reports $PS_OUT/reports -reports "0,1,2,3,4,5,6,7,8"
     
     # #########################################
     # ## Clean reports for easy R parsing
     # #########################################
      
     cd $PS_OUT/reports
     
     for FILE in *.txt
     do
       cat $FILE | sed 's/#/Number_/g' | sed "s/'/prime/g" > ${FILE%.*}_clean.txt 
       #tail -n +2 $FILE  | sed "s/'/prime/g" >> ${FILE%.*}_clean.txt 
     done

  rm $FILEPATH

  fi
  ((i++))
done

