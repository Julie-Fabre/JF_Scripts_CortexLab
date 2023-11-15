%% load + save data 
loadVids = 0;
[task_data, session_data] = cl_loadPerStimulusData('task', 2, loadVids); 
if loadVids
    save('/home/julie/Dropbox/MATLAB/task_data_goNogo_session.mat', '-struct', 'session_data', '-v7.3');
end
save('/home/julie/Dropbox/MATLAB/task_data_goNogo2.mat', '-struct', 'task_data', '-v7.3');
close all;

keep loadVids
[task_data_gogogo, session_data] = cl_loadPerStimulusData('taskGo', 2, loadVids); 
if loadVids
    save('/home/julie/Dropbox/MATLAB/task_data_gogogo_session.mat', '-struct', 'session_data', '-v7.3');
end
save('/home/julie/Dropbox/MATLAB/task_data_gogogo2.mat', '-struct', 'task_data_gogogo', '-v7.3');
close all;

keep loadVids
[task_data_passive, session_data] = cl_loadPerStimulusData('passive', 5, loadVids);  
if loadVids
    save('/home/julie/Dropbox/MATLAB/task_data_passive_session.mat', '-struct', 'session_data', '-v7.3');
end
save('/home/julie/Dropbox/MATLAB/task_data_passive_full2.mat', '-struct', 'task_data_passive', '-v7.3');
close all;

keep loadVids
[nat_passive, session_data] =  cl_loadPerStimulusData('passive', 6, loadVids); 
%save('/home/julie/Dropbox/MATLAB/nat_passive_full.mat', '-struct', 'nat_passive', '-v7.3');
close all;

keep loadVids
[gratings_passive, session_data] = cl_loadPerStimulusData('passive', 1, loadVids); 
%save('/home/julie/Dropbox/MATLAB/gratings_passive_full.mat', '-struct', 'gratings_passive', '-v7.3');


%% passive

% check where the snr and gpe stim resp come from -> partuicular regions?
% subregions? animals? 

%% task behavior 


%% task data
cl_plot_PSTHs; 

cl_visResp_summarize;
cl_dprime_summarize;

cl_decoding_summarize;
cl_dimReduction_summarize; % PCA + angle; tSNE
cl_modeling_summarize; 

% to do:
% - location
% - cell types 
% - reward omitted vs not 
% - does performance affect encoding? 
% - movement 
% - types of responses 

%% suppl: histology

%% suppl: bombcell 