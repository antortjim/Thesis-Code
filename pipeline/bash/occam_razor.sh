#! /bin/bash
ROOT_DIR=$1
EXP_NAME=$2
OUTPUT_MOFF=$3

# Absolute path this script is in, thus /home/user/bin
SOFT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SOFT_DIR/load_flags.sh $ROOT_DIR $EXP_NAME > /dev/null 2>&1

Rscript $ROOT_DIR/scripts/R/smallestUniqueGroups.R --output_moff $OUTPUT_MOFF
