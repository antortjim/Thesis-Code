#! /bin/bash

SOFT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SOFT_DIR/load_flags.sh > /dev/null 2>&1

######################################################################
## Perform match between runs and feature extraction using moFF
#######################################################################
EXP_DESIGN=$ROOT_DIR/$EXP_NAME/data/experimental_design.tsv
#echo $EXP_DESIGN
SAMPLE_INDEX=$(sed -n $'1s/\t/\\\n/gp' $EXP_DESIGN | grep -nx 'sample' | cut -f 1 -d:)
#echo $SAMPLE_INDEX
SAMPLES=$(cut -f $SAMPLE_INDEX -d$'\t' $EXP_DESIGN | tail -n +2 | sort | uniq)
cut -f $SAMPLE_INDEX -d$'\t' $EXP_DESIGN | tail -n +2 | sort | uniq | echo

#echo ${SAMPLES[*]}

i=1
for S in ${SAMPLES[*]}
do
  echo $S
  SAMPLE_REPORTS=$(python $ROOT_DIR/scripts/Python/extract_files_path.py --sample $S --exp_design $EXP_DESIGN --prepend $ROOT_DIR/$EXP_NAME/$PS_OUT/reports --append "_PS_EXTENSION")
  SAMPLE_MZML=$(python $ROOT_DIR/scripts/Python/extract_files_path.py --sample $S --exp_design $EXP_DESIGN --prepend $ROOT_DIR/$EXP_NAME/data/mzML --append ".mzML")
  echo $SAMPLE_REPORTS
  echo $SAMPLE_MZML
  #python $MOFF_PATH/moff_all.py --inputtsv $SAMPLE_REPORTS --inputraw $SAMPLE_MZML --tol 10 --ext txt --log_file_name all_workflow.log --output_folder output_moff --peptide_summary 1
  ((i++))
  if [ $i -ge 1 ]
  then
    break
  fi
done

#python moff_all.py --inputtsv $ROOT_DIR/$EXP_NAME/$PS_OUT/reports/20080311_CPTAC6_07_6A005.txt $ROOT_DIR/$EXP_NAME/$PS_OUT/reports/20080313_CPTAC6_07_6A005.txt $ROOT_DIR/$EXP_NAME/$PS_OUT/reports/20080315_CPTAC6_07_6A005.txt --inputraw input_raw/20080311_CPTAC6_07_6A005.mzML input_raw/20080313_CPTAC6_07_6A005.mzML input_raw/20080315_CPTAC6_07_6A005.mzML --tol 10 --ext txt --log_file_name all_workflow.log --output_folder output_moff --peptide_summary 1
