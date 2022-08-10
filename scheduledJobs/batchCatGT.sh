#! /bin/bash

#HOME=/home/julie/bin
#ACCESS="${HOME}/access.log"
#ERROR="${HOME}/error.log"

#ls ~/Dropbox/MATLAB/onPaths/JF_Scripts_CortexLab/scheduledJobs/KSerror.log >> $OUTPUT 2>&1

#echo "yo bb" >> ~/Dropbox/MATLAB/onPaths/JF_Scripts_CortexLab/scheduledJobs/KS.log


/home/julie/Downloads/CatGTLnxApp/CatGT-linux/runit.sh '-dir=/home/netshare/zaru/JF082/2022-08-09/ephys/site1/ -run=2022-08-09_JF082 -g=0 -t=0,Inf -ap -prb=0 -gblcar -no_run_fld'
/home/julie/Downloads/CatGTLnxApp/CatGT-linux/runit.sh '-dir=/home/netshare/zaru/JF082/2022-08-09/ephys/site2/ -run=2022-08-09_JF082 -g=0 -t=0,Inf -ap -prb=1 -gblcar -no_run_fld'

/home/julie/Downloads/CatGTLnxApp/CatGT-linux/runit.sh '-dir=/home/netshare/zaru/JF084/2022-08-09/ephys/site1/ -run=2022-08-09_JF084 -g=1 -t=0,Inf -ap -prb=0 -gblcar -no_run_fld'
/home/julie/Downloads/CatGTLnxApp/CatGT-linux/runit.sh '-dir=/home/netshare/zaru/JF084/2022-08-09/ephys/site2/ -run=2022-08-09_JF084 -g=1 -t=0,Inf -ap -prb=1 -gblcar -no_run_fld'
/home/julie/Downloads/CatGTLnxApp/CatGT-linux/runit.sh '-dir=/home/netshare/zaru/JF084/2022-08-09/ephys/site3/ -run=2022-08-09_JF084-2 -g=0 -t=0,Inf -ap -prb=0 -gblcar -no_run_fld'
/home/julie/Downloads/CatGTLnxApp/CatGT-linux/runit.sh '-dir=/home/netshare/zaru/JF084/2022-08-09/ephys/site4/ -run=2022-08-09_JF084-2 -g=0 -t=0,Inf -ap -prb=1 -gblcar -no_run_fld'

#/home/julie/Downloads/CatGTLnxApp/CatGT-linux/runit.sh '-dir=/home/netshare/zaru/JF072/2022-08-07/ephys/site1/ -run=2022-08-07_JF072 -g=0 -t=0,Inf -ap -prb=0 -gblcar -no_run_fld'
#/home/julie/Downloads/CatGTLnxApp/CatGT-linux/runit.sh '-dir=/home/netshare/zaru/JF072/2022-08-07/ephys/site2/ -run=2022-08-06_JF072-1 -g=1 -t=0,Inf -ap -prb=1 -gblcar -no_run_fld'

#matlab -desktop -r "sendLatestPlotsEmail"
