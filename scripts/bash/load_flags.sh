#! /bin/bash

SOFT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PIPELINE_SETTINGS="$SOFT_DIR/../../pipeline_settings.txt"
NSETTINGS=$(wc -l $PIPELINE_SETTINGS | awk '{print $1}')
echo $NSETTINGS
i=0
while [ $i -le $NSETTINGS ]
do
  eval "$(head -$i $PIPELINE_SETTINGS | tail -n 1 |cut -f1 -d:)=$(head -$i $PIPELINE_SETTINGS | tail -n 1 | cut -f2 -d:)"
  ((i++))
done
IFS=', ' read -r -a DATABASE_NAMES <<< "$DATABASE_NAMES"


