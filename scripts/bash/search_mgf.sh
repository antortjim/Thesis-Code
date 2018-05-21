#! /bin/bash


SOFT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SOFT_DIR/load_flags.sh > /dev/null 2>&1
 
########################################
## Search!
########################################

if [ ! -z "$(ls -A $ROOT_DIR/$EXP_NAME/$SPECTRA)" ]
then

  # Read spectra files in alphabetic order
  SPECTRA_ARRAY=$(ls $ROOT_DIR/$EXP_NAME/$SPECTRA/*mgf)
  i=1
  for SPECTRUM_FILE in ${SPECTRA_ARRAY[*]}
  do

    echo "Searching"
    SEARCHGUI_OUT=$SAMPLE_NAME_out
    UPDATED_SEARCH_ENGINES=$SEARCH_ENGINES

    #######################################
    ## Check already searched
    #######################################

    for SE in ${SEARCH_ENGINES[*]}
    do
      echo "Checking $SE is not already done"
      UPDATED_SEARCH_ENGINES=$(python $ROOT_DIR/scripts/Python/update_search_engines.py --updated_list ${UPDATED_SEARCH_ENGINES[*]} -se $SE --searchgui_out $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT --sample_name $SAMPLE_NAME)
    done
      UPDATED_SEARCH_ENGINES=$(python $ROOT_DIR/scripts/Python/definitive_list_search_engines.py --updated_list ${UPDATED_SEARCH_ENGINES[*]})
      echo $UPDATED_SEARCH_ENGINES
    
     echo "Using: ${UPDATED_SEARCH_ENGINES[*]}"
#    java -Djava.awt.headless=true -cp $SEARCHGUI_PATH eu.isas.searchgui.cmd.SearchCLI \
#      -spectrum_files $ROOT_DIR/$EXP_NAME/$SPECTRA/$SPECTRUM_FILE \
#      -output_folder $ROOT_DIR/$EXP_NAME \
#      -id_params $ROOT_DIR/$EXP_NAME/$SETTINGS_DIR/$PARAMS_NAME.par \
#      -output_data 1 \
#      $UPDATED_SEARCH_ENGINES
    exit
    mv $ROOT_DIR/$EXP_NAME/searchgui_out.zip $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT.zip
    unzip $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT.zip -d $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT
    mv $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT $ROOT_DIR/$EXP_NAME/old_searches

    ##Compatibility with TPP.. not working yet
    containsElement "msgf" ${SEARCHGUI_ENGINES[@]}
    if [ $? ]
    then
      #idconvert $MZID --pepXML -v -o $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/$SAMPLE_NAME.msgf.mzid
      #sed 's/|\(\w*\)_REVERSED|/|REVERSED_\1|/' $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/$SAMPLE_NAME.pepXML > $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/${SAMPLE_NAME}_prefix.pepXML
      echo $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/$SAMPLE_NAME.msgf.mzid >> $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/identification_files_$SAMPLE_NAME.txt
    fi

    containsElement "comet" ${SEARCHGUI_ENGINES[@]}
    if [ $? ]
    then
      #sed 's/|\(\w*\)_REVERSED|/|REVERSED_\1|/' $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/$SAMPLE_NAME.comet.pep.xml > $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/${SAMPLE_NAME}_prefix.comet.pep.xml
      echo $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/$SAMPLE_NAME.comet.pep.xml >> $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/identification_files_$SAMPLE_NAME.txt
    fi

    #find $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT -maxdepth 1 | grep "${SAMPLE_NAME}_prefix" | grep -v -e "\.mgf$"  > $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/identification_files_$SAMPLE_NAME.txt
    zip -j $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/identification_files_$SAMPLE_NAME.zip -@ < $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/identification_files_$SAMPLE_NAME.txt

    echo "`date` searchGUI finished for sample $SAMPLE_NAME with engines ${UPDATED_SEARCH_ENGINES[*]}" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
    echo "Launching call PeptideShaker" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
    $SOFT_DIR/call_peptide_shaker.sh $ROOT_DIR/$EXP_NAME/$SPECTRA/$SPECTRUM_FILE $i $SEARCHGUI_OUT > $ROOT_DIR/$EXP_NAME/$EXP_NAME_peptide_shaker_$SAMPLE_NAME.out
    echo "`date` Peptide-Shaker finished for sample $SAMPLE_NAME" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
    ((i++))
  done
fi
