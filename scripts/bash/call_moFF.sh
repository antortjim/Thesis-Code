#! /bin/bash

source activate moFF

SOFT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SOFT_DIR/load_flags.sh $1 $2 > /dev/null 2>&1

######################################################################
## Perform match between runs and feature extraction using moFF
#######################################################################
EXP_DESIGN=$ROOT_DIR/$EXP_NAME/data/experimental_design.tsv
# The variable grouping samples into replicates i.e the variable where replicates have a shared value must be in the third column
GROUP_INDEX=2
# The variable giving the sample names must be in the first column
SAMPLE_NAMES_INDEX=1
EXPERIMENT_GROUPS=$(awk '{print $'$GROUP_INDEX'}' $EXP_DESIGN | tail -n +2 | sort | uniq | xargs)
REPORTS_DIR="PSM_reports"
echo "`date` call_moFF.sh EXPERIMENT_GROUPS: ${EXPERIMENT_GROUPS[*]}" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log


i=1
for S in ${EXPERIMENT_GROUPS[*]}
do
  #echo "`date` call_moFF.sh Running grep $S $EXP_DESIGN | awk '{print $4}' | xargs | sed -e 's/ /|/g'" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
  SAMPLE_NAMES=$(grep $S $EXP_DESIGN | awk '{print $'$SAMPLE_NAMES_INDEX'}' | xargs | sed -e 's/ /|/g')
  SAMPLE_NAMES_ARRAY=($(grep $S $EXP_DESIGN | awk '{print $'$SAMPLE_NAMES_INDEX'}' | xargs))
  CMD="python $ROOT_DIR/scripts/Python/check_reports_exist.py --files ${SAMPLE_NAMES_ARRAY[@]} --root_dir $ROOT_DIR --exp_name $EXP_NAME --reports_dir $REPORTS_DIR"
  echo "`date` call_moFF.sh Analyzing samples in condition $S " >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
  MOFF_OUTPUT=output_moff_RAW
  REPORTS_EXIST=$(python $ROOT_DIR/scripts/Python/check_reports_exist.py --files ${SAMPLE_NAMES_ARRAY[@]} --root_dir $ROOT_DIR --exp_name $EXP_NAME --reports_dir $REPORTS_DIR --suffix \"\")
  echo `date` call_moFF.sh python $ROOT_DIR/scripts/Python/check_reports_exist.py --files ${SAMPLE_NAMES_ARRAY[@]} --root_dir $ROOT_DIR --exp_name $EXP_NAME --reports_dir $REPORTS_DIR --suffix \"\" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
  MOFF_DONE=$(python $ROOT_DIR/scripts/Python/check_reports_exist.py --files ${SAMPLE_NAMES_ARRAY[@]} --root_dir $ROOT_DIR --exp_name $EXP_NAME --reports_dir  $REPORTS_DIR/$MOFF_OUTPUT --suffix "_match_moff_result")
  REPORTS_EXIST=1
  #CMD=python $ROOT_DIR/scripts/Python/check_reports_exist.py --files ${SAMPLE_NAMES_ARRAY[@]} --root_dir $ROOT_DIR --exp_name $EXP_NAME --reports_dir  $REPORTS_DIR/$MOFF_OUTPUT --suffix "_match_moff_result"
  #echo "`date` call_moFF.sh CMD $CMD" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
  #MOFF_DONE=$(eval $CMD)
  #MOFF_DONE=$(python $ROOT_DIR/scripts/Python/check_reports_exist.py --files ${SAMPLE_NAMES_ARRAY[@]} --root_dir $ROOT_DIR --exp_name $EXP_NAME --reports_dir  $MOFF_OUTPUT --suffix "_match_moff_result")
  echo "`date` call_moFF.sh Are the reports available $REPORTS_EXIST ?" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
  echo "`date` call_moFF.sh Is moFF done already? $MOFF_DONE" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
  RUN_MOFF=0
  if [ $REPORTS_EXIST -ne $MOFF_DONE ] && [ $REPORTS_EXIST -eq 1 ]; then RUN_MOFF=1; else RUN_MOFF=0; fi
  echo "`date` call_moFF.sh Calling moFF? $RUN_MOFF" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log

  if [ $RUN_MOFF -eq 1 ] 
  then
    echo "`date` call_moFF.sh SAMPLE_NAMES: ${SAMPLE_NAMES[*]}" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
    SAMPLE_REPORTS=$(python $ROOT_DIR/scripts/Python/extract_files_path.py --group $S --exp_design $EXP_DESIGN --prepend $ROOT_DIR/$EXP_NAME/$PS_OUT/$REPORTS_DIR --append ".txt")
    #SAMPLE_MZML=$(python $ROOT_DIR/scripts/Python/extract_files_path.py --group $S --exp_design $EXP_DESIGN --prepend $ROOT_DIR/$EXP_NAME/data/mzML --append ".mzML")
    SAMPLE_RAW=$(python $ROOT_DIR/scripts/Python/extract_files_path.py --group $S --exp_design $EXP_DESIGN --prepend $ROOT_DIR/$EXP_NAME/data/RAW --append ".RAW")
    echo "`date` call_moFF.sh Calling full moFF workflow (MBR and Apex)" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log

    CMD="python $MOFF_PATH/moff_all.py \
      --tol 10 \
      --output_folder $ROOT_DIR/$EXP_NAME/$PS_OUT/$REPORTS_DIR/$MOFF_OUTPUT \
      --inputraw $SAMPLE_RAW
      --inputtsv $SAMPLE_REPORTS
      --rt_w 3.0 --tag_pep_sum_file moFF_run --weight_comb 1 \
      --out_filt 1 --filt_width 1.5 --peptide_summary 1 --log_file_name $EXP_NAME"

    echo "`date` call_moFF.sh $CMD" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
    eval $CMD

  else
    echo "`date` call_moFF.sh Omitting sample $S" >> $ROOT_DIR/$EXP_NAME/log/pipeline.log
  fi

  #source deactivate moFF
  #Rscript --vanilla $ROOT_DIR/scripts/R/custom_PSM_report.R --root_dir $ROOT_DIR --exp_name $EXP_NAME --sample_names ${SAMPLE_NAMES[*]}
  #source activate moFF

done

## Run the Occam's razor
#$ROOT_DIR/scripts/bash/occam_razor.sh $ROOT_DIR $EXP_NAME $MOFF_OUTPUT


