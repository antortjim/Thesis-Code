NSETTINGS=$(wc -l pipeline_settings.txt | awk '{print $1}')
echo $NSETTINGS
i=0
while [ $i -le $NSETTINGS ]
do
  eval "$(head -$i pipeline_settings.txt | tail -n 1 |cut -f1 -d:)=$(head -$i pipeline_settings.txt | tail -n 1 | cut -f2 -d:)"
  ((i++))
done
IFS=', ' read -r -a DATABASE_NAMES <<< "$DATABASE_NAMES"


