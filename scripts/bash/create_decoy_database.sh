#! /bin/bash

SOFT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo $SOFT_DIR
source $SOFT_DIR/load_flags.sh $1 $2 > /dev/null 2>&1

#######################################
## Prepare database: create decoy
########################################

if [ ! -f databases/all.fasta ] 
then
  echo "Preparing decoy"
  rm -rf $ROOT_DIR/$EXP_NAME/databases/all*
  
  if [ $NZ_DB -eq 1 ]
  then
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
  else
    for DB in  $ROOT_DIR/$EXP_NAME/databases/*fasta
    do
      echo $DB
      FILE=$(basename -- "$DB")
      FILENAME="${FILE%.*}"
      if [ ${FILENAME:0:2} != "all" ]
      then
       echo "`date` create_decoy_database.sh: Adding database $FILENAME" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
       cat $DB >> $ROOT_DIR/$EXP_NAME/databases/all.fasta
      fi
    done
  fi

  CMD="java -cp $SEARCHGUI_PATH eu.isas.searchgui.cmd.FastaCLI -in $ROOT_DIR/$EXP_NAME/databases/all.fasta -decoy"
  echo "`date` create_decoy_database.sh: Running $CMD" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
  eval $CMD
fi
