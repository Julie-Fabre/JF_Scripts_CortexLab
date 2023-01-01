
rerunQM = 1;
JF_runBombcell('JF070', '2022-06-09', 1, [],rerunQM)
JF_runBombcell('JF070', '2022-06-10', 1, [],rerunQM)
JF_runBombcell('JF070', '2022-06-11', 1, [],rerunQM)
JF_runBombcell('JF070', '2022-06-12', 1, [],rerunQM)

rerunKS = 1;
rerunQM = 1;
JF_runKSandBombcell('JF070', '2022-06-09', 2, rerunKS,rerunQM);%
JF_runKSandBombcell('JF070', '2022-06-10', 2, rerunKS,rerunQM);
JF_runKSandBombcell('JF070', '2022-06-11', 2, rerunKS,rerunQM);%
JF_runKSandBombcell('JF070', '2022-06-13', 1, rerunKS,rerunQM);
JF_runKSandBombcell('JF070', '2022-06-16', 2, rerunKS,rerunQM);
JF_runKSandBombcell('JF070', '2022-06-17', 1, rerunKS,rerunQM);
JF_runKSandBombcell('JF070', '2022-06-18', 1, rerunKS,rerunQM);

% AP_preprocess_phase3_newOEJF('JF070','2022-06-11', 1,[],'neuropixPhase3A_kilosortChanMap.mat');
% AP_preprocess_phase3_newOEJF('JF070','2022-06-12', 1,[],'neuropixPhase3A_kilosortChanMap.mat');
