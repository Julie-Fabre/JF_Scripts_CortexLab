%% KS 3A 
%AP_preprocess_phase3_newOEJF('JF043','2021-06-25', 1)
%AP_preprocess_phase3_newOEJF('JF044','2021-06-25', 1)
%AP_preprocess_phase3_newOEJF('JF043','2021-06-27', 1)
%AP_preprocess_phase3_newOEJF('JF044','2021-06-27', 1)
%AP_preprocess_phase3_newOEJF('JF043','2021-06-28', 1)
JF_preprocess_NPX2('JF062', '2021-11-05', 'chanMapNP2_4Shank_bottRow_flipper0x2Emat_kilosortChanMap.mat', 1, 1, []);
JF_preprocess_NPX2('JF062', '2021-11-05', 'chanMapNP2_4Shank_bottRow_flipper0x2Emat_kilosortChanMap.mat', 1, 2, []);
JF_preprocess_NPX2('JF062', '2021-11-06', 'chanMapNP2_4Shank_bottRow_flipper0x2Emat_kilosortChanMap.mat', 1, 1, []);
JF_preprocess_NPX2('JF062', '2021-11-07', 'chanMapNP2_4Shank_bottRow_flipper0x2Emat_kilosortChanMap.mat', 1, 1, []);
JF_preprocess_NPX2('JF062', '2021-11-07', 'chanMapNP2_4Shank_bottRow_flipper0x2Emat_kilosortChanMap.mat', 1, 2, []);
JF_preprocess_NPX2('JF062', '2021-11-07', 'chanMapNP2_4Shank_bottRow_flipper0x2Emat_kilosortChanMap.mat', 1, 3, []);
JF_preprocess_NPX2('JF062', '2021-11-07', 'chanMapNP2_4Shank_bottRow_flipper0x2Emat_kilosortChanMap.mat', 1, 4, []);
% %% extract 2.0 
% [ephysAPfile,~] = AP_cortexlab_filenameJF('JF043','2021-06-27',2,'ephys_ap',2,[]);
% [ephysKSfile,~] = AP_cortexlab_filenameJF('JF043','2021-06-27',2,'ephys',2,[]);
% syncFT(ephysAPfile, 385, ephysKSfile)
% clear all;
% 
% [ephysAPfile,~] = AP_cortexlab_filenameJF('JF044','2021-06-27',2,'ephys_ap',2,[]);
% [ephysKSfile,~] = AP_cortexlab_filenameJF('JF044','2021-06-27',2,'ephys',2,[]);
% syncFT(ephysAPfile, 385, ephysKSfile)
% clear all;
% 
% [ephysAPfile,~] = AP_cortexlab_filenameJF('JF044','2021-06-28',2,'ephys_ap',2,[]);
% [ephysKSfile,~] = AP_cortexlab_filenameJF('JF044','2021-06-28',2,'ephys',2,[]);
% syncFT(ephysAPfile, 385, ephysKSfile)
% clear all;