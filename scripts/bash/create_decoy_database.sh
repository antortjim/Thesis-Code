#! /bin/bash

SOFT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SOFT_DIR/load_flags.sh $1 $2 > /dev/null 2>&1

#######################################
## Prepare database: create decoy
########################################

if [ ! -f databases/all.fasta ] 
then
  echo "Preparing decoy"
  rm -rf $ROOT_DIR/$EXP_NAME/databases/all*
  
  for DB in ${DATABASE_NAMES[*]}
  do
      echo $DB
      cp /z/ms/mascot/mascot-2.5/sequence/$DB/current/*.fasta  $ROOT_DIR/$EXP_NAME/databases/$DB.fasta
      if [ "$DB" == "NZ_Contaminants" ]
      then
        cat $ROOT_DIR/$EXP_NAME/databases/NZ_Contaminants.fasta | sed -e 's/>con_\([0-9A-Z]\)/>generic\|\1/g' | sed 's/ /\|/' | sed 's/|/_CONTAMINANT|/2' > $ROOT_DIR/$EXP_NAME/databases/NZ_Contaminants_uniprot.fasta
        cat $ROOT_DIR/$EXP_NAME/databases/NZ_Contaminants_uniprot.fasta >> $ROOT_DIR/$EXP_NAME/databases/all.fasta
      else
        cat $ROOT_DIR/$EXP_NAME/databases/$DB.fasta >> $ROOT_DIR/$EXP_NAME/databases/all.fasta
      fi
  done

  CMD="java -cp $SEARCHGUI_PATH eu.isas.searchgui.cmd.FastaCLI -in $ROOT_DIR/$EXP_NAME/databases/all.fasta -decoy"
  echo "`date` create_decoy_database.sh: Running $CMD" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
  #eval $CMD
fi
