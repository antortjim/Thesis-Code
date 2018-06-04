#! /bin/bash

SOFT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SOFT_DIR/load_flags.sh > /dev/null 2>&1


#########################################
## Data processing
#########################################

EXPERIMENTAL_DESIGN="$ROOT_DIR/$EXP_NAME/data/experimental_design.tsv"
SPECTRUM_FILE=$1
i=$2
SEARCHGUI_OUT=$3

#echo $SPECTRUM_FILE
SAMPLE_ID="sample_$i"
FILENAME=$(basename $SPECTRUM_FILE)
#echo $FILENAME
SAMPLE_NAME="${FILENAME%.*}"
#echo $SAMPLE_NAME

if [ -f $EXPERIMENTAL_DESIGN ]
then
  CONDITION_NAME=$(grep $SAMPLE_NAME $EXPERIMENTAL_DESIGN | awk '{print $5}')
  REPLICATE=$(grep $SAMPLE_NAME $EXPERIMENTAL_DESIGN | awk '{print $1}')
else
  CONDITION_NAME="condition_1"
  REPLICATE=$i
fi

echo -e "${SAMPLE_NAME}\t${EXP_NAME}_${CONDITION_NAME}_${REPLICATE}"  >> $ROOT_DIR/$EXP_NAME/$PS_OUT/sample_mapping.txt

#echo $CONDITION_NAME
#echo $REPLICATE

MZID=$ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/$SAMPLE_NAME.msgf.mzid
echo $MZID
COMET_OUTPUT=$ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/$SAMPLE_NAME.comet.pep.xml
echo $COMET_OUTPUT

#echo "`date` $MZID $COMET_OUTPUT" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
#echo "`date` $PS_OUT" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log

#if [ ! -f "$ROOT_DIR/$EXP_NAME/$PS_OUT/$SAMPLE_NAME.cpsx" ] && ([ -f $MZID ] || [ -f $COMET_OUTPUT ])
if [ ! -f "$ROOT_DIR/$EXP_NAME/$PS_OUT/$SAMPLE_NAME.cpsx" ] || [ $FORCE_PS = "1" ]
then

  echo "`date` call_peptide_shaker.sh Calling peptide-shaker passing identification files in $SEARCHGUI_OUT/identification_files_$SAMPLE_NAME.zip" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
  echo "`date` call_peptide_shaker.sh Spectrum file: $SPECTRUM_FILE" >>  $ROOT_DIR/$EXP_NAME/log/pipeline.log
  echo "`date` call_peptide_shaker.sh Call:   java -cp $PEPTIDESHAKER_PATH eu.isas.peptideshaker.cmd.PeptideShakerCLI \
     -experiment $EXP_NAME \
     -sample $CONDITION_NAME \
     -replicate $REPLICATE \
     -identification_files $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/identification_files_$SAMPLE_NAME.zip \
     -spectrum_files $SPECTRUM_FILE \
     -id_params $ROOT_DIR/$EXP_NAME/$SETTINGS_DIR/$PARAMS_NAME.par \
     -out $ROOT_DIR/$EXP_NAME/$PS_OUT/$SAMPLE_NAME.cpsx"  >> $ROOT_DIR/$EXP_NAME/log/pipeline.log



   java -cp $PEPTIDESHAKER_PATH eu.isas.peptideshaker.cmd.PeptideShakerCLI \
     -experiment $EXP_NAME \
     -sample $CONDITION_NAME \
     -replicate $REPLICATE \
     -identification_files $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/identification_files_$SAMPLE_NAME.zip \
     -spectrum_files $SPECTRUM_FILE \
     -id_params $ROOT_DIR/$EXP_NAME/$SETTINGS_DIR/$PARAMS_NAME.par \
     -out $ROOT_DIR/$EXP_NAME/$PS_OUT/$SAMPLE_NAME.cpsx >> $ROOT_DIR/$EXP_NAME/log/pipeline.log


  echo "`date` call_peptide_shaker.sh Outputing to $PS_OUT/$SAMPLE_NAME.cpsx" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
else
  echo "`date` call_peptide_shaker.sh Skipping Peptide-Shaker" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
fi
if [ -f "$ROOT_DIR/$EXP_NAME/$PS_OUT/$SAMPLE_NAME.cpsx" ]
then
   # #########################################
   # ## Report 
   # #########################################
    
   echo "`date` call_peptide_shaker.sh Reporting for sample $SAMPLE_NAME" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
   java -cp $PEPTIDESHAKER_PATH eu.isas.peptideshaker.cmd.ReportCLI -in $ROOT_DIR/$EXP_NAME/$PS_OUT/$SAMPLE_NAME.cpsx \
     -out_reports $ROOT_DIR/$EXP_NAME/$PS_OUT/reports -reports "0,1,2,3,4,5,6,7,8"

   # Change prefix of reports to match $SAMPLE_NAME
   OLD_PREFIX=${EXP_NAME}_${CONDITION_NAME}_${REPLICATE}
   NEW_PREFIX=$SAMPLE_NAME
   cd $ROOT_DIR/$EXP_NAME/$PS_OUT/reports/
   for REPORT in ${OLD_PREFIX}*.txt
   do
      mv $REPORT $NEW_PREFIX${REPORT#$OLD_PREFIX}
   done

   # Generate custom report including MS1 intensity
   echo "`date` call_peptide_shaker.sh Creating custom_PSM_report for sample $SAMPLE_NAME" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
   ln ${SAMPLE_NAME}_Default_PSM_Report.txt ../PSM_reports/${SAMPLE_NAME}_Default_PSM_Report.txt
   
   # #########################################
   # ## Clean reports for easy R parsing
   # #########################################
    
   # for FILE in *.txt
   # do
   #   if [ $(grep "_clean.txt" <<< $FILE) ]
   #   then
   #     cat $FILE | sed 's/#/Number_/g' | sed "s/'/prime/g" > ${FILE%.*}_clean.txt 
   #     #tail -n +2 $FILE  | sed "s/'/prime/g" >> ${FILE%.*}_clean.txt 
   #   fi
   # done
   cd $ROOT_DIR
fi

