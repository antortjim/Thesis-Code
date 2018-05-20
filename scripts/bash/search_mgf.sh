#! /bin/bash


SOFT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SOFT_DIR/load_flags.sh > /dev/null 2>&1

#######################################
## Check already searched
#######################################

for MZID in $ROOT_DIR/$EXP_NAME/searchgui_out/*mzid
do
  echo $MZID
  FILENAME=$(basename $MZID)
  # 2 times because msgf+ output has 2 extensions msgf.mzid
  SAMPLE_NAME="${FILENAME%.*}"
  SAMPLE_NAME="${SAMPLE_NAME%.*}"
  echo $SAMPLE_NAME
  echo $ROOT_DIR/$EXP_NAME/$SPECTRA/$SAMPLE_NAME.mgf
done

 
########################################
## Search!
########################################

if [ ! -z "$(ls -A $ROOT_DIR/$EXP_NAME/$SPECTRA)" ]
then
  echo "Searching"
  java -Djava.awt.headless=true -cp $SEARCHGUI_PATH eu.isas.searchgui.cmd.SearchCLI \
    -spectrum_files $ROOT_DIR/$EXP_NAME/$SPECTRA \
    -output_folder $ROOT_DIR/$EXP_NAME \
    -id_params $ROOT_DIR/$EXP_NAME/$SETTINGS_DIR/$PARAMS_NAME.par \
    -output_data 1 \
    -comet 1
#   # -comet 1 -myrimatch 1 -xtandem 1
  unzip $ROOT_DIR/$EXP_NAME/searchgui_out.zip -d $ROOT_DIR/$EXP_NAME/searchgui_out
  mv $ROOT_DIR/$EXP_NAME/searchgui_out.zip $ROOT_DIR/$EXP_NAME/old_searches
fi
