#! /bin/bash

#HOME=/home/julie/bin
#ACCESS="${HOME}/access.log"
#ERROR="${HOME}/error.log"

ls ~/Dropbox/MATLAB/onPaths/JF_Scripts_CortexLab/scheduledJobs/KSerror.log >> $OUTPUT 2>&1

echo "yo bb" >> ~/Dropbox/MATLAB/onPaths/JF_Scripts_CortexLab/scheduledJobs/KS.log


FILEDATE=$(date +'%Y_%m_%d')
DATE=$(date +'%Y-%m-%d')
DIRSTART="/home/netshare/zinu/JF067/"
SEP="/"
ASTERIX="*"
SPIKEGLXNAME=$DIRSTART$DATE$SEP$FILEDATE$ASTERIX
EPHYSNAME=$DIRSTART$DATE$SEP$"ephys"
mv $SPIKEGLXNAME $EPHYSNAME

SITE1SLX=$EPHYSNAME$SEP$FILEDATE$ASTERIX
SITE1NAME=$EPHYSNAME$SEP$"site1"
mv $SITE1SLX $SITE1NAME

DIREND="/ephys/site1/"
DIRFULL="$DIRSTART$DATE$DIREND"
#echo "$DIRFULL"

FILEEND="-JF067"
FILEDATE=$(date +'%Y_%m_%d')
FILEFULL="$FILEDATE$FILEEND"

/home/julie/Downloads/CatGTLnxApp/CatGT-linux/runit.sh '-dir=$DIRFULL -run=$FILEFULL -g=0 -t=0,Inf -ap -prb=0 -gblcar -no_run_fld'

matlab -desktop -r "sendLatestPlotsEmail"


