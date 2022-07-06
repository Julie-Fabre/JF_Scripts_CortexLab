%% load ephys data, align to laser onsets 

%% ~~ loading info ~~
close all;
myPaths;
animal = 'XG006';
date = '2022-06-30';
site = 1;
recording = [];
experiment = 1;
verbose = true;
load_parts.cam=false;
load_parts.ephys=true;

%% ~~ extract sync channel ~~

[ephysAPfile,aa] = AP_cortexlab_filenameJF(animal,date,experiment,'ephys_ap',site,recording);
if size(ephysAPfile,2) ==2 %keep only ap
    ephysAPfile = ephysAPfile{1};
end
isSpikeGlx = contains(ephysAPfile, '_g');
if isSpikeGlx
     [ephysKSfile,~] = AP_cortexlab_filenameJF(animal,date,experiment,'ephys',site,recording);
    if isempty(dir([ephysKSfile filesep 'sync.mat']))
        sync = syncFT(ephysAPfile, 385, ephysKSfile);
    end
end


%% ~~ load data from experiment ~~ 
loadClusters = 1;% whether to load phy results
JF_loadExperiment_forXin;


%% ~~ plot data in GUI ~~ 
curr_shank=NaN;
AP_cellrasterJF_forXin({laser_on_flip_times,laser_flip_times(1:2:end),laser_on_flip_times}, ...
    {laserParamsAllLaserOn.Amp, laserParamsAllLaserOn.Freq, laserParamsAllLaserOn.Ramp});