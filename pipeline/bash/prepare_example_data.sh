#! /bin/bash

SOFT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SOFT_DIR/load_flags.sh $1 $2 > /dev/null 2>&1
$SOFT_DIR/load_flags.sh $1 $2 --source-only


rm $ROOT_DIR/$EXP_NAME/data/mgf_input/*
NUMBER_SPECTRA=1000
for MGF_FILE in $ROOT_DIR/thp1/data/mgf/*.mgf
do
  echo $(basename ${MGF_FILE%.*})
  MGF_FILE_DEST=$ROOT_DIR/$EXP_NAME/data/mgf/$(basename $MGF_FILE)
  #echo $MGF_FILE_DEST
  catNOccurrences $NUMBER_SPECTRA $MGF_FILE > $MGF_FILE_DEST
  SYMBOLIC_LINK=$ROOT_DIR/$EXP_NAME/data/mgf_input/$(basename $MGF_FILE)
  #echo $SYMBOLIC_LINK
  ln -s $MGF_FILE_DEST $SYMBOLIC_LINK
  msconvert $SYMBOLIC_LINK --mzML  -o $ROOT_DIR/$EXP_NAME/data/mzML/ --outfile $(basename ${MGF_FILE%.*}).mzML > /dev/null
done

cd $ROOT_DIR/$EXP_NAME/data/mgf_input 
rm ./PD7502-GDTHP1-A_C1.mgf
rm ./PD7503-GDTHP1-A_C2.mgf
rm ./PD7504-GDTHP1-A_C1.mgf
rm ./PD7508-GDTHP1-A_C1.mgf
rm ./PD7509-GDTHP1-A_C2.mgf
rm ./PD7510-GDTHP1-A_C1.mgf
rm ./PD7511-GDTHP1-A_C2.mgf
rm ./PD7512-GDTHP1-A_C1.mgf
rm ./PD7513-GDTHP1-A_C2.mgf
cd -
