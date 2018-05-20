#! /bin/bash

#DATABASE_NAMES=${1:-"Homo_sapiens NZ_Contaminants"}
#IFS=', ' read -r -a DATABASE_NAMES <<< "$DATABASE_NAMES"
#SPECTRA=${2:-data/mgf_input}
#PARAMS_NAME=${3:-thp1}
#SEARCHGUI_PATH=${4:-/z/home/aoj/opt/SearchGUI-3.3.1/SearchGUI-3.3.1.jar}
#PEPTIDESHAKER_PATH=${5:-/z/home/aoj/opt/PeptideShaker-1.16.23/PeptideShaker-1.16.23.jar}
#EXP_NAME=${6:-thp1}
#PS_OUT=${7:-peptideShaker_out}
#SETTINGS_DIR=${8:-settings}
#ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source load_flags.sh > /dev/null 2>&1

echo DATABASE_NAMES
echo ${DATABASE_NAMES[*]}
echo SPECTRA
echo $SPECTRA
echo PARAMS_NAME
echo $PARAMS_NAME
echo SEARCHGUI_PATH
echo $SEARCHGUI_PATH
echo PEPTIDESHAKER_PATH
echo $PEPTIDESHAKER_PATH
echo EXP_NAME
echo $EXP_NAME
echo PS_OUT
echo $PS_OUT
echo SETTINGS_DIR
echo $SETTINGS_DIR
echo ROOT_DIR
echo $ROOT_DIR

