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
## Prepare database: create decoy
########################################

if [ ! -f databases/all.fasta ] 
then
  echo "Preparing decoy"
  rm -rf databases/all.fasta
  
  for DB in ${DATABASE_NAMES[*]}
  do
      cp /z/ms/mascot/mascot-2.5/sequence/$DB/current/*.fasta  $ROOT_DIR/databases/$DB.fasta
      if [ "$DB" == "NZ_Contaminants" ]
      then
        cat $ROOT_DIR/databases/NZ_Contaminants.fasta | sed -e 's/>con_\([0-9A-Z]\)/>generic\|\1/g' | sed 's/ /\|/' | sed 's/|/_CONTAMINANT|/2' > $ROOT_DIR/databases/NZ_Contaminants_uniprot.fasta
        cat $ROOT_DIR/databases/NZ_Contaminants_uniprot.fasta >> $ROOT_DIR/databases/all.fasta
      else
        cat $ROOT_DIR/databases/$DB.fasta >> $ROOT_DIR/databases/all.fasta
      fi
  done

java -cp $SEARCHGUI_PATH eu.isas.searchgui.cmd.FastaCLI -in $ROOT_DIR/databases/all.fasta -decoy
fi
