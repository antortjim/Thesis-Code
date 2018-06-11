ROOT_DIR=$1
EXP_NAME=$2

Rscript $ROOT_DIR/scripts/R/load_maxlfq_data.R
Rscript $ROOT_DIR/scripts/R/process_peptide_intensities.R --root_dir $ROOT_DIR --exp_name $EXP_NAME --input_dir "quantification" --output_dir "quantification" --pipeline FALSE
Rscript $ROOT_DIR/scripts/Python/LFQ.py --input_dir "quantification" --output_dir "quantification"
Rscript $ROOT_DIR/scripts/R/fold_changes.R  --plots_dir plots
