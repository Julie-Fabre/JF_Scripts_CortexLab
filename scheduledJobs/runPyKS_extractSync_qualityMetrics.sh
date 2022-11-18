#! /bin/bash

#HOME=/home/julie/bin
#ACCESS="${HOME}/access.log"
#ERROR="${HOME}/error.log"

ls ~/Dropbox/MATLAB/onPaths/JF_Scripts_CortexLab/scheduledJobs/KSerror.log >> $OUTPUT 2>&1

echo "yo bb" >> ~/Dropbox/MATLAB/onPaths/JF_Scripts_CortexLab/scheduledJobs/KS.log


FILEDATE=$(date +'%Y-%m-%d')
DATE=$(date +'%Y-%m-%d')

#conda activate pyks2

#python ~/Dropbox/Python/pyks/runPyKS.py 

#conda deactivate pyks2

#matlab -desktop -r "sendLatestPlotsEmail"
