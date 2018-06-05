#! /bin/bash

ROOT_DIR=$1
EXP_NAME=$2
#SOFT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SOFT_DIR=$ROOT_DIR/scripts/bash
PIPELINE_SETTINGS="$ROOT_DIR/$EXP_NAME/pipeline_settings.txt"
NSETTINGS=$(wc -l $PIPELINE_SETTINGS | awk '{print $1}')
echo $NSETTINGS
i=0
while [ $i -le $NSETTINGS ]
do
  CMD="$(head -$i $PIPELINE_SETTINGS | tail -n 1 |cut -f1 -d:)=$(head -$i $PIPELINE_SETTINGS | tail -n 1 | cut -f2 -d:)"
  echo $CMD
  eval $CMD
  ((i++))
done

IFS=', ' read -r -a DATABASE_NAMES <<< "$DATABASE_NAMES"
IFS=', ' read -r -a SEARCH_ENGINES <<< "$SEARCH_ENGINES"
echo ${SEARCH_ENGINES[*]}

# From https://stackoverflow.com/questions/3685970/check-if-a-bash-array-contains-a-value
containsElement () {
  local e match="$1"
  shift
  for e; do [[ "$e" == "$match" ]] && return 0; done
  return 1
}

# https://www.linuxquestions.org/questions/programming-9/grep-till-second-occurance-of-pattern-770754/
catNOccurrences () {
  NUMBER_SPECTRA=$1
  MGF_FILE="$2"
  awk -vf=0 '/END IONS/{f=1;++d}d!='$NUMBER_SPECTRA'{print}f&&d=='$NUMBER_SPECTRA'{print;exit}' $MGF_FILE 
  return 1
}

main() {
  echo "Load flags sourced"
}
if [ "${3}" != "--source-only" ]; then
    main "${@}"
fi
