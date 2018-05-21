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
  SPECTRA_ARRAY=$(ls $ROOT_DIR/$EXP_NAME/$SPECTRA/*mgf)
  i=1
  for SPECTRUM_FILE in ${SPECTRA_ARRAY[*]}
  do
    echo "Searching"
    SEARCHGUI_OUT=$SAMPLE_NAME_out
    java -Djava.awt.headless=true -cp $SEARCHGUI_PATH eu.isas.searchgui.cmd.SearchCLI \
      -spectrum_files $ROOT_DIR/$EXP_NAME/$SPECTRA/$SPECTRUM_FILE \
      -output_folder $ROOT_DIR/$EXP_NAME \
      -id_params $ROOT_DIR/$EXP_NAME/$SETTINGS_DIR/$PARAMS_NAME.par \
      -output_data 1 \
      $SEARCHGUI_ENGINES
    mv $ROOT_DIR/$EXP_NAME/searchgui_out.zip $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT.zip
    unzip $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT.zip -d $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT
    mv $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT $ROOT_DIR/$EXP_NAME/old_searches

    containsElement "msgf" ${SEARCHGUI_ENGINES[@]}
    if [ $? ]
    then
      idconvert $MZID --pepXML -v -o $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/$SAMPLE_NAME.msgf.mzid
      sed 's/|\(\w*\)_REVERSED|/|REVERSED_\1|/' $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/$SAMPLE_NAME.pepXML > $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/${SAMPLE_NAME}_prefix.pepXML
    fi

    containsElement "comet" ${SEARCHGUI_ENGINES[@]}
    if [ $? ]
    then
      sed 's/|\(\w*\)_REVERSED|/|REVERSED_\1|/' $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/$SAMPLE_NAME.comet.pep.xml > $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/${SAMPLE_NAME}_prefix.comet.pep.xml
    fi

    echo "`date` searchGUI finished for sample $SAMPLE_NAME with engines $SEARCH_ENGINES" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
    echo "Launching call PeptideShaker" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
    $SOFT_DIR/call_peptide_shaker.sh $ROOT_DIR/$EXP_NAME/$SPECTRA/$SPECTRUM_FILE $i $SEARCHGUI_OUT > $ROOT_DIR/$EXP_NAME/$EXP_NAME_peptide_shaker_$SAMPLE_NAME.out
    echo "`date` Peptide-Shaker finished for sample $SAMPLE_NAME" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
    ((i++))
  done
fi
