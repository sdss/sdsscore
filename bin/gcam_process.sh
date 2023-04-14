#!/bin/bash
echo " "
echo "--- Running gcam_process ----"
date
echo " "

source /etc/profile.d/sdss5.sh
module load sdsscore
/home/sdss5/software/sdsscore/main/bin/gcam_process --plot
