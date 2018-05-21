#! /bin/bash

SOFT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SOFT_DIR/load_flags.sh > /dev/null 2>&1


#########################################
## Data processing
#########################################

EXPERIMENTAL_DESIGN="$ROOT_DIR/$EXP_NAME/data/experimental_design.tsv"
FILEPATH=$1
i=$2
SEARCHGUI_OUT=$3

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

MZID=$ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/$SAMPLE_NAME.msgf.mzid
COMET_OUTPUT=$ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/$SAMPLE_NAME.comet.pep.xml

if [ ! -f "$PS_OUT/$SAMPLE_NAME.cpsx" ] && {[ -f $MZID ] || [ -f $COMET_OUTPUT ]}
then
  echo "Calling peptideshaker"

  #echo $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/identification_files_$SAMPLE_NAME.txt

  echo $FILEPATH

   java -cp $PEPTIDESHAKER_PATH eu.isas.peptideshaker.cmd.PeptideShakerCLI \
     -experiment $EXP_NAME \
     -sample $CONDITION_NAME \
     -replicate $REPLICATE  \
     -identification_files $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/identification_files_$SAMPLE_NAME.zip \
     -spectrum_files \"$FILEPATH\" \
     -id_params $ROOT_DIR/$EXP_NAME/$SETTINGS_DIR/$PARAMS_NAME.par \
     -out $PS_OUT/$SAMPLE_NAME.cpsx

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

