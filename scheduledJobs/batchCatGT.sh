#! /bin/bash

#HOME=/home/julie/bin
#ACCESS="${HOME}/access.log"
#ERROR="${HOME}/error.log"

#ls ~/Dropbox/MATLAB/onPaths/JF_Scripts_CortexLab/scheduledJobs/KSerror.log >> $OUTPUT 2>&1

#echo "yo bb" >> ~/Dropbox/MATLAB/onPaths/JF_Scripts_CortexLab/scheduledJobs/KS.log


/home/julie/Downloads/CatGTLnxApp/CatGT-linux/runit.sh '-dir=/home/netshare/zaru/JF082/2022-07-28/ephys/site1/ -run=2022-07-28_JF082 -g=0 -t=0,Inf -ap -prb=0 -gblcar -no_run_fld'
/home/julie/Downloads/CatGTLnxApp/CatGT-linux/runit.sh '-dir=/home/netshare/zaru/JF082/2022-07-28/ephys/site2/ -run=2022-07-28_JF082 -g=0 -t=0,Inf -ap -prb=1 -gblcar -no_run_fld'



#matlab -desktop -r "sendLatestPlotsEmail"
