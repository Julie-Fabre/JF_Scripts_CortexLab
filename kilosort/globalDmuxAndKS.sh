
runit.sh '-dir=/ -run=run_name -no_run_fld -prb=0./runit.sh '-dir=/home/netshare/zinu/JF067/2022-02-16/ephys/site1/ -run=2022_02_16-JF067 -g=0 -t=0,Inf -ap -prb=0 -no_run_fld -gblcar'

matlab -nodisplay -nosplash -nodesktop -r "run('/home/julie/Dropbox/MATLAB/onPaths/bombcell/qualityMetrics/runKS.m');exit;" | tail -n +11
