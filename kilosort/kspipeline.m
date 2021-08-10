%% KS 3A 
%AP_preprocess_phase3_newOEJF('JF043','2021-06-25', 1)
%AP_preprocess_phase3_newOEJF('JF044','2021-06-25', 1)
%AP_preprocess_phase3_newOEJF('JF043','2021-06-27', 1)
%AP_preprocess_phase3_newOEJF('JF044','2021-06-27', 1)
%AP_preprocess_phase3_newOEJF('JF043','2021-06-28', 1)
AP_preprocess_phase3_newOEJF('JF044','2021-07-01', 1)

%% extract 2.0 
[ephysAPfile,~] = AP_cortexlab_filenameJF('JF043','2021-06-27',2,'ephys_ap',2,[]);
[ephysKSfile,~] = AP_cortexlab_filenameJF('JF043','2021-06-27',2,'ephys',2,[]);
syncFT(ephysAPfile, 385, ephysKSfile)
clear all;

[ephysAPfile,~] = AP_cortexlab_filenameJF('JF044','2021-06-27',2,'ephys_ap',2,[]);
[ephysKSfile,~] = AP_cortexlab_filenameJF('JF044','2021-06-27',2,'ephys',2,[]);
syncFT(ephysAPfile, 385, ephysKSfile)
clear all;

[ephysAPfile,~] = AP_cortexlab_filenameJF('JF044','2021-06-28',2,'ephys_ap',2,[]);
[ephysKSfile,~] = AP_cortexlab_filenameJF('JF044','2021-06-28',2,'ephys',2,[]);
syncFT(ephysAPfile, 385, ephysKSfile)
clear all;