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
