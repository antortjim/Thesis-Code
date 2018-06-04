#! /bin/bash

source activate moFF

SOFT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SOFT_DIR/load_flags.sh > /dev/null 2>&1

######################################################################
## Perform match between runs and feature extraction using moFF
#######################################################################
EXP_DESIGN=$ROOT_DIR/$EXP_NAME/data/experimental_design.tsv
#echo $EXP_DESIGN
SAMPLE_INDEX=$(sed -n $'1s/\t/\\\n/gp' $EXP_DESIGN | grep -nx 'sample' | cut -f 1 -d:)
#FILE_NAME_INDEX=$(sed -n $'1s/\t/\\\n/gp' $EXP_DESIGN | grep -nx 'file' | cut -f 1 -d:)

#awk '{print $5}' $EXP_DESIGN | tail -n +2 

echo "`date` call_moFF.sh SAMPLE_NAME_INDEX: $FILE_NAME_INDEX" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
echo "`date` call_moFF.sh FILE_NAME_INDEX: $FILE_NAME_INDEX" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
SAMPLES=$(cut -f $SAMPLE_INDEX -d$'\t' $EXP_DESIGN | tail -n +2 | sort | uniq)
#SAMPLE_NAMES=$(cut -f $FILE_NAME_INDEX -d$'\t' $EXP_DESIGN | tail -n +2)


i=1
for S in ${SAMPLES[*]}
do
  echo $S
  echo "`date` call_moFF.sh Running grep $S $EXP_DESIGN | awk '{print $4}' | cut -f 1 -d. | xargs | sed -e 's/ /|/g'" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
  SAMPLE_NAMES=$(grep $S $EXP_DESIGN | awk '{print $4}' | cut -f 1 -d. | xargs | sed -e 's/ /|/g')
  SAMPLE_NAMES_ARRAY=($(grep $S $EXP_DESIGN | awk '{print $4}' | cut -f 1 -d. | xargs))
  OK=$(python $ROOT_DIR/scripts/Python/check_paths.py --files ${SAMPLE_NAMES_ARRAY[@]} --root_dir $ROOT_DIR --exp_name $EXP_NAME)
  echo "`date` call_moFF.sh $OK" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log

  if [ $OK -eq 1 ] 
  then
    echo "`date` call_moFF.sh SAMPLE_NAMES: ${SAMPLE_NAMES[*]}" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
    #SAMPLE_REPORTS=$(python $ROOT_DIR/scripts/Python/extract_files_path.py --sample $S --exp_design $EXP_DESIGN --prepend $ROOT_DIR/$EXP_NAME/$PS_OUT/reports --append "_Custom_PSM_Report.txt")
    #SAMPLE_MZML=$(python $ROOT_DIR/scripts/Python/extract_files_path.py --sample $S --exp_design $EXP_DESIGN --prepend $ROOT_DIR/$EXP_NAME/data/mzML --append ".mzML")
    echo $SAMPLE_REPORTS
    echo $SAMPLE_MZML
    echo "`date` call_moFF.sh Calling moFF match between runs module" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
    echo "`date` call_moFF.sh python $MOFF_PATH/moff_mbr.py --inputF $ROOT_DIR/$EXP_NAME/peptideShaker_out/PSM_reports --sample "$SAMPLE_NAMES" --ext txt --log_file_name mbr.log" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
    python $MOFF_PATH/moff_mbr.py --inputF $ROOT_DIR/$EXP_NAME/peptideShaker_out/PSM_reports --sample "$SAMPLE_NAMES" --ext txt --log_file_name mbr.log
  else
    echo "`date` call_moFF.sh Omitting sample $S" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
  fi
done

source deactivate moFF
Rscript --vanilla $ROOT_DIR/scripts/R/custom_PSM_report.R --root_dir $ROOT_DIR --exp_name $EXP_NAME --sample_names $SAMPLE_NAMES

