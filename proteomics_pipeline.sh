#! /bin/bash

DATABASE_NAMES=${1:-"Homo_sapiens NZ_Contaminants"}
IFS=', ' read -r -a DATABASE_NAMES <<< "$DATABASE_NAMES"
SPECTRA=${2:-data/mgf_input}
PARAMS_NAME=${3:-thp1}
SEARCHGUI_PATH=${4:-/z/home/aoj/opt/SearchGUI-3.3.1/SearchGUI-3.3.1.jar}
PEPTIDESHAKER_PATH=${5:-/z/home/aoj/opt/PeptideShaker-1.16.23/PeptideShaker-1.16.23.jar}
EXP_NAME=${6:-test_experiment}
PS_OUT=${7:-peptideShaker_out}
SETTINGS_DIR=${8:-settings}
ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR=${9:-$ROOT_DIR}

########################################
## Prepare database: create decoy
########################################

if [ ! -f databases/all.fasta ] 
then
  echo "Preparing decoy"
  rm -rf databases/all.fasta
  
  for DB in ${DATABASE_NAMES[*]}
  do
      cp /z/ms/mascot/mascot-2.5/sequence/$DB/current/*.fasta  databases/$DB.fasta
      if [ "$DB" == "NZ_Contaminants" ]
      then
        cat databases/NZ_Contaminants.fasta | sed -e 's/>con_\([0-9A-Z]\)/>generic\|\1/g' | sed 's/ /\|/' | sed 's/|/_CONTAMINANT|/2' > databases/NZ_Contaminants_uniprot.fasta
        cat databases/NZ_Contaminants_uniprot.fasta >> databases/all.fasta
      else
        cat databases/$DB.fasta >> databases/all.fasta
      fi
  done

java -cp $SEARCHGUI_PATH eu.isas.searchgui.cmd.FastaCLI -in $ROOT_DIR/databases/all.fasta -decoy
fi

#########################################
#Create search settings file
#########################################

# 0: specific, 1: semispecific, one of the ends could not be specific
# Allowed missed cleavages
# prec_tol in ppm
# frag_tol in Da
# prec_tol in ppm
# frag_tol in Da

if [ ! -f $ROOT_DIR/$SETTINGS_DIR/$PARAMS_NAME.par ]
then
  echo "Preparing settings file"

  java -cp $SEARCHGUI_PATH eu.isas.searchgui.cmd.IdentificationParametersCLI -out $SETTINGS_DIR/$PARAMS_NAME \
    -db $ROOT_DIR/databases/all_concatenated_target_decoy.fasta \
    -enzyme 'Trypsin' \
    -specificity '1' \
    -fixed_mods 'Carbamidomethylation of C' \
    -variable_mods 'Oxidation of M, Deamidation of N, Deamidation of Q' \
    -mc '2' \
    -prec_tol '10' \
    -frag_tol '0.5'\
    -fi 'b' \
    -ri 'y'
fi
  
########################################
## Search!
########################################

if [ ! -z "$(ls -A $ROOT_DIR/$SPECTRA)" ]
then
  echo "Searching"
  java -Djava.awt.headless=true -cp $SEARCHGUI_PATH eu.isas.searchgui.cmd.SearchCLI \
    -spectrum_files $ROOT_DIR/$SPECTRA \
    -output_folder $ROOT_DIR \
    -id_params $SETTINGS_DIR/$PARAMS_NAME.par \
    -output_data 1 \
    -msgf 1
  # -comet 1 -myrimatch 1 -xtandem 1
  unzip searchgui_out.zip -d searchgui_out
  mv searchgui_out.zip old_searches
fi


#########################################
## Data processing
#########################################

i=1
for FILEPATH in $ROOT_DIR/$SPECTRA/*.mgf
do
  SAMPLE_ID="sample_$i"
  FILENAME=$(basename $FILEPATH)
  SAMPLE_NAME=${FILENAME%.*}
  echo $FILEPATH
  echo $FILENAME
  echo $SAMPLE_NAME

  CONDITION_NAME=$(grep $SAMPLE_NAME  $ROOT_DIR/$SPECTRA/../experimental_design.tsv | awk '{print $5}')
  REPLICATE=$(grep $SAMPLE_NAME  $ROOT_DIR/$SPECTRA/../experimental_design.tsv | awk '{print $1}')
  
  echo $CONDITION_NAME
  echo $REPLICATE

  #IDENTIFICATION_FILES=$(find `pwd`/searchgui_out -maxdepth 1 | grep "$SAMPLE_NAME" | grep -v -e "\.mgf$" )
  #IDENTIFICATION_FILES_STRING="$(echo $IDENTIFICATION_FILES | sed 's/ /, /g')"

  if [ ! -f "$PS_OUT/$SAMPLE_NAME.cpsx" ] && [ -f $ROOT_DIR/searchgui_out/$SAMPLE_NAME.msgf.mzid ]
  then
    echo "Calling peptideshaker"
    find $ROOT_DIR/searchgui_out -maxdepth 1 | grep "$SAMPLE_NAME" | grep -v -e "\.mgf$"  > searchgui_out/identification_files_$SAMPLE_NAME.txt
    zip -j $ROOT_DIR/searchgui_out/identification_files_$SAMPLE_NAME.zip -@ < $ROOT_DIR/searchgui_out/identification_files_$SAMPLE_NAME.txt

    echo $FILEPATH

    java -cp $PEPTIDESHAKER_PATH eu.isas.peptideshaker.cmd.PeptideShakerCLI \
      -experiment $EXP_NAME \
      -sample $CONDITION_NAME \
      -replicate $REPLICATE  \
      -identification_files $ROOT_DIR/searchgui_out/identification_files_$SAMPLE_NAME.zip \
      -spectrum_files \"$FILEPATH\" \
      -id_params $ROOT_DIR/$SETTINGS_DIR/$PARAMS_NAME.par \
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

