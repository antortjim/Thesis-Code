#! /bin/bash

source load_flags.sh > /dev/null 2>&1
 
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
    -msgf 1
  # -comet 1 -myrimatch 1 -xtandem 1
  unzip $ROOT_DIR/$EXP_NAME/searchgui_out.zip -d $ROOT_DIR/$EXP_NAME/searchgui_out
  mv $ROOT_DIR/$EXP_NAME/searchgui_out.zip $ROOT_DIR/$EXP_NAME/old_searches
fi
