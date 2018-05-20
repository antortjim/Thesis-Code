#! /bin/bash

DATABASE_NAMES=${1:-"Homo_sapiens NZ_Contaminants"}
IFS=', ' read -r -a DATABASE_NAMES <<< "$DATABASE_NAMES"
SPECTRA=${2:-data/mgf_input}
PARAMS_NAME=${3:-thp1}
SEARCHGUI_PATH=${4:-/z/home/aoj/opt/SearchGUI-3.3.1/SearchGUI-3.3.1.jar}
PEPTIDESHAKER_PATH=${5:-/z/home/aoj/opt/PeptideShaker-1.16.23/PeptideShaker-1.16.23.jar}
EXP_NAME=${6:-thp1}
PS_OUT=${7:-peptideShaker_out}
SETTINGS_DIR=${8:-settings}
ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR=${9:-$ROOT_DIR/thp1}

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
