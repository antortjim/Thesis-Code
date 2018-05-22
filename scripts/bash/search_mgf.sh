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
  for SPECTRUM_PATH in ${SPECTRA_ARRAY[*]}
  do
    SPECTRUM_FILE=$(basename $SPECTRUM_PATH)

    SAMPLE_NAME=$(basename ${SPECTRUM_FILE%.*})
    echo $SAMPLE_NAME
    SEARCHGUI_OUT=${SAMPLE_NAME}_out
    echo $SEARCHGUI_OUT
    UPDATED_SEARCH_ENGINES=$SEARCH_ENGINES
 
    #######################################
    ## Check already searched
    #######################################
 
    for SE in ${SEARCH_ENGINES[*]}
    do
      echo "`date` Checking $SE is not already done" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
      UPDATED_SEARCH_ENGINES=$(python $ROOT_DIR/scripts/Python/update_search_engines.py --updated_list "${UPDATED_SEARCH_ENGINES[*]}" --se $SE --searchgui_out $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT --sample_name $SAMPLE_NAME)
    done
    

   if [ $UPDATED_SEARCH_ENGINES != "None"  ]
   then
 
      UPDATED_SEARCH_ENGINES=$(python $ROOT_DIR/scripts/Python/definitive_list_search_engines.py --updated_list "${UPDATED_SEARCH_ENGINES[*]}")
      echo $UPDATED_SEARCH_ENGINES
      
      echo "Searching"
#      echo "`date` searchGUI started for sample $SAMPLE_NAME with engines ${UPDATED_SEARCH_ENGINES[*]}" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
#      java -Djava.awt.headless=true -cp $SEARCHGUI_PATH eu.isas.searchgui.cmd.SearchCLI \
#        -spectrum_files $ROOT_DIR/$EXP_NAME/$SPECTRA/$SPECTRUM_FILE \
#        -output_folder $ROOT_DIR/$EXP_NAME \
#        -id_params $ROOT_DIR/$EXP_NAME/$SETTINGS_DIR/$PARAMS_NAME.par \
#        -output_data 1 \
#        $UPDATED_SEARCH_ENGINES

      mv $ROOT_DIR/$EXP_NAME/searchgui_out.zip $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT.zip
      unzip -o $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT.zip -d $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT
    else
      echo "`date` Skipping search for $SAMPLE_NAME" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
    fi

    rm $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/identification_files_$SAMPLE_NAME.txt
    touch $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/identification_files_$SAMPLE_NAME.txt


    ## Add identification files for each search engine

    ##Compatibility with TPP.. not working yet
    containsElement "msgf" ${SEARCH_ENGINES[@]}
    if [ $? ]
    then
      #idconvert $MZID --pepXML -v -o $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/$SAMPLE_NAME.msgf.mzid
      #sed 's/|\(\w*\)_REVERSED|/|REVERSED_\1|/' $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/$SAMPLE_NAME.pepXML > $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/${SAMPLE_NAME}_prefix.pepXML
      echo $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/$SAMPLE_NAME.msgf.mzid >> $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/identification_files_$SAMPLE_NAME.txt
    fi

    containsElement "comet" ${SEARCH_ENGINES[@]}
    if [ $? ]
    then
      #sed 's/|\(\w*\)_REVERSED|/|REVERSED_\1|/' $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/$SAMPLE_NAME.comet.pep.xml > $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/${SAMPLE_NAME}_prefix.comet.pep.xml
      echo $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/$SAMPLE_NAME.comet.pep.xml >> $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/identification_files_$SAMPLE_NAME.txt
    fi

    #find $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT -maxdepth 1 | grep "${SAMPLE_NAME}_prefix" | grep -v -e "\.mgf$"  > $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/identification_files_$SAMPLE_NAME.txt
    zip -j $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/identification_files_$SAMPLE_NAME.zip -@ < $ROOT_DIR/$EXP_NAME/$SEARCHGUI_OUT/identification_files_$SAMPLE_NAME.txt

    echo "`date` searchGUI finished for sample $SAMPLE_NAME with engines ${SEARCH_ENGINES[*]}" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
    echo "`date` Launching call PeptideShaker" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
    $SOFT_DIR/call_peptide_shaker.sh $ROOT_DIR/$EXP_NAME/$SPECTRA/$SPECTRUM_FILE $i $SEARCHGUI_OUT > $ROOT_DIR/$EXP_NAME/$PS_OUT/$SAMPLE_NAME.out
    echo "`date` Peptide-Shaker finished for sample $SAMPLE_NAME" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
    ((i++))
  done
fi
