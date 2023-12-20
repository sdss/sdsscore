#!/bin/bash
echo " "
echo "--- Running telemetry ----"
date
echo " "

source /uufs/chpc.utah.edu/sys/modulefiles/scripts/module_init/module_init.sh
module load sdsscore/holtz

#APO
for file in `ls -t $APOGEE_TELEMETRY_N/APOGEE_PD* | head -5` 
do 
  echo $file
  $SDSSCORE_DIR/bin/telemetry --load $file -o apo --dir ' '
done
$SDSSCORE_DIR/bin/telemetry -o apo --html $APOGEE_TELEMETRY/apo.html --skip 12

#LCO
for file in `ls -t $APOGEE_TELEMETRY_S/AS_PD* | head -5` 
do 
  echo $file
  $SDSSCORE_DIR/bin/telemetry --load $file -o lco --dir ' '
done
$SDSSCORE_DIR/bin/telemetry -o lco --html $APOGEE_TELEMETRY/lco.html --skip 12
