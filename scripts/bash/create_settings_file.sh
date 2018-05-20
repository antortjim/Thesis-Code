#! /bin/bash

source load_flags.sh > /dev/null 2>&1

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

  java -cp $SEARCHGUI_PATH eu.isas.searchgui.cmd.IdentificationParametersCLI -out $ROOT_DIR/$EXP_NAME/$SETTINGS_DIR/$PARAMS_NAME \
    -db $ROOT_DIR/$EXP_NAME/databases/all_concatenated_target_decoy.fasta \
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
