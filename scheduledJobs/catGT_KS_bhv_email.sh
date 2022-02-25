#! /bin/bash

#HOME=/home/julie/bin
#ACCESS="${HOME}/access.log"
#ERROR="${HOME}/error.log"

ls ~/Dropbox/MATLAB/onPaths/JF_Scripts_CortexLab/scheduledJobs/KSerror.log >> $OUTPUT 2>&1

echo "yo bb" >> ~/Dropbox/MATLAB/onPaths/JF_Scripts_CortexLab/scheduledJobs/KS.log

DIRSTART="/home/netshare/zinu/JF067/"
DATE=$(date +'%Y-%m-%d')
DIREND="/ephys/site1/"
DIRFULL="$DIRSTART$DATE$DIREND"
#echo "$DIRFULL"

FILEEND="-JF067"
FILEDATE=$(date +'%Y_%m_%d')
FILEFULL="$FILEDATE$FILEEND"

/home/julie/Downloads/CatGTLnxApp/CatGT-linux/runit.sh '-dir=$DIRFULL -run=$FILEFULL -g=0 -t=0,Inf -ap -prb=0 -gblcar -no_run_fld'

matlab -desktop -r "sendLatestPlotsEmail"


