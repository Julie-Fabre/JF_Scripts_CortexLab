
%% manual stuff
% allen coordinates
allen_ccf_npx;

%% AP080 - average response:
animals = {'AP080'};
curr_animal = 1; % (set which animal to use)
corona = 0;
animal = animals{curr_animal};

protocol = 'JF_Grating'; % (this is the name of the Signals protocol)
experiments = AP_find_experimentsJF(animal, protocol, true); %experiments = experiments([experiments.imaging] & [experiments.ephys]); % (use only experiments with both widefield + ephys)
experiments = experiments([experiments.ephys]);

curr_day = 3; % (set which day to use)

day = experiments(curr_day).day; % date
thisDay = experiments(curr_day).day; % date
thisDate = thisDay;
experiment = experiments(curr_day).experiment; % experiment number

verbose = true; % display load progress and some info figures
load_parts.cam = false;


load_parts.imaging = false;
load_parts.ephys = true;

%AP_load_experiment; % loads data from experiment
site = 2;
experiment = 2;
AP_load_experimentJF;

try
    AP_cellrasterJF({stimOn_times}, ...
        {trial_conditions(:, 1) + trial_conditions(:, 2) + trial_conditions(:, 3), ...
        trial_conditions(:, 4)});
catch
    try
        AP_cellrasterJF({stimOn_times}, ...
            {trial_conditions(:, 1) + trial_conditions(:, 2) + trial_conditions(:, 3), ...
            });
    catch
        AP_cellrasterJF({stimOn_times}, ...
            {stimIDs});
    end
end
dataGrating(1).spike_templates = spike_templates;
dataGrating(1).spike_times_timeline = spike_times_timeline;
dataGrating(1).stimOn = stimOn_times;
dataGrating(1).stimCond = trial_conditions;
dataGrating(1).depth = [1000, 3000];
dataGrating(1).spike_depths = spike_depths;

%% AP082 - average response:
animals = {'AP082'};
curr_animal = 1; % (set which animal to use)
corona = 0;
animal = animals{curr_animal};

protocol = 'JF_Grating'; % (this is the name of the Signals protocol)
experiments = AP_find_experimentsJF(animal, protocol, true); %experiments = experiments([experiments.imaging] & [experiments.ephys]); % (use only experiments with both widefield + ephys)
experiments = experiments([experiments.ephys]);

curr_day = 1; % (set which day to use)

day = experiments(curr_day).day; % date
thisDay = experiments(curr_day).day; % date
thisDate = thisDay;
experiment = experiments(curr_day).experiment; % experiment number

verbose = true; % display load progress and some info figures
load_parts.cam = false;
load_parts.imaging = false;
load_parts.ephys = true;

%AP_load_experiment; % loads data from experiment
site = 1;
experiment = 1;
AP_load_experimentJF;

try
    AP_cellrasterJF({stimOn_times}, ...
        {trial_conditions(:, 1) + trial_conditions(:, 2) + trial_conditions(:, 3), ...
        trial_conditions(:, 4)});
catch
    try
        AP_cellrasterJF({stimOn_times}, ...
            {trial_conditions(:, 1) + trial_conditions(:, 2) + trial_conditions(:, 3), ...
            });
    catch
        AP_cellrasterJF({stimOn_times}, ...
            {stimIDs});
    end
end
dataGrating(2).spike_templates = spike_templates;
dataGrating(2).spike_times_timeline = spike_times_timeline;
dataGrating(2).stimOn = stimOn_times;
dataGrating(2).stimCond = trial_conditions;
dataGrating(2).stimIDs = stimIDs;
dataGrating(2).depth = [1500, 3000];
dataGrating(2).spike_depths = spike_depths;

%% JF017 - average response:
animals = {'JF017'};
curr_animal = 1; % (set which animal to use)
corona = 0;
animal = animals{curr_animal};

protocol = 'JF_Grating'; % (this is the name of the Signals protocol)
experiments = AP_find_experimentsJF(animal, protocol, true); %experiments = experiments([experiments.imaging] & [experiments.ephys]); % (use only experiments with both widefield + ephys)
experiments = experiments([experiments.ephys]);

curr_day = 1; % (set which day to use)

day = experiments(curr_day).day; % date
thisDay = experiments(curr_day).day; % date
thisDate = thisDay;
experiment = experiments(curr_day).experiment; % experiment number

verbose = true; % display load progress and some info figures
load_parts.cam = false;
load_parts.imaging = false;
load_parts.ephys = true;

%AP_load_experiment; % loads data from experiment
site = 2;
experiment = 4;
AP_load_experimentJF;

try
    AP_cellrasterJF({stimOn_times}, ...
        {trial_conditions(:, 1) + trial_conditions(:, 2) + trial_conditions(:, 3), ...
        trial_conditions(:, 4)});
catch
    try
        AP_cellrasterJF({stimOn_times}, ...
            {trial_conditions(:, 1) + trial_conditions(:, 2) + trial_conditions(:, 3), ...
            });
    catch
        AP_cellrasterJF({stimOn_times}, ...
            {stimIDs});
    end
end

dataGrating(3).spike_templates = spike_templates;
dataGrating(3).spike_times_timeline = spike_times_timeline;
dataGrating(3).stimOn = stimOn_times;
dataGrating(3).stimCond = trial_conditions;
dataGrating(3).stimIDs = stimIDs;
dataGrating(3).depth = [1200, 2200];
dataGrating(3).spike_depths = spike_depths;

%% JF019-1 - average response:
animals = {'JF019'};
curr_animal = 1; % (set which animal to use)
corona = 0;
animal = animals{curr_animal};

protocol = 'JF_Grating'; % (this is the name of the Signals protocol)
experiments = AP_find_experimentsJF(animal, protocol, true); %experiments = experiments([experiments.imaging] & [experiments.ephys]); % (use only experiments with both widefield + ephys)
experiments = experiments([experiments.ephys]);

curr_day = 2; % (set which day to use)

day = experiments(curr_day).day; % date
thisDay = experiments(curr_day).day; % date
thisDate = thisDay;
experiment = experiments(curr_day).experiment; % experiment number

verbose = true; % display load progress and some info figures
load_parts.cam = false;
load_parts.imaging = false;
load_parts.ephys = true;

%AP_load_experiment; % loads data from experiment
site = 1;
experiment = 1;
AP_load_experimentJF;

try
    AP_cellrasterJF({stimOn_times}, ...
        {trial_conditions(:, 1) + trial_conditions(:, 2) + trial_conditions(:, 3), ...
        trial_conditions(:, 4)});
catch
    try
        AP_cellrasterJF({stimOn_times}, ...
            {trial_conditions(:, 1) + trial_conditions(:, 2) + trial_conditions(:, 3), ...
            });
    catch
        AP_cellrasterJF({stimOn_times}, ...
            {stimIDs});
    end
end

dataGrating(4).spike_templates = spike_templates;
dataGrating(4).spike_times_timeline = spike_times_timeline;
dataGrating(4).stimOn = stimOn_times;
dataGrating(4).stimCond = trial_conditions;
dataGrating(4).stimIDs = stimIDs;
dataGrating(4).depth = [2000, 4000];
dataGrating(4).spike_depths = spike_depths;

%% JF019 - 2
animals = {'JF019'};
curr_animal = 1; % (set which animal to use)
corona = 0;
animal = animals{curr_animal};

protocol = 'JF_Grating'; % (this is the name of the Signals protocol)
experiments = AP_find_experimentsJF(animal, protocol, true); %experiments = experiments([experiments.imaging] & [experiments.ephys]); % (use only experiments with both widefield + ephys)
experiments = experiments([experiments.ephys]);

curr_day = 1; % (set which day to use)

day = experiments(curr_day).day; % date
thisDay = experiments(curr_day).day; % date
thisDate = thisDay;
experiment = experiments(curr_day).experiment; % experiment number

verbose = true; % display load progress and some info figures
load_parts.cam = false;
load_parts.imaging = false;
load_parts.ephys = true;

%AP_load_experiment; % loads data from experiment
site = 3;
experiment = 7;
AP_load_experimentJF;

try
    AP_cellrasterJF({stimOn_times}, ...
        {trial_conditions(:, 1) + trial_conditions(:, 2) + trial_conditions(:, 3), ...
        trial_conditions(:, 4)});
catch
    try
        AP_cellrasterJF({stimOn_times}, ...
            {trial_conditions(:, 1) + trial_conditions(:, 2) + trial_conditions(:, 3), ...
            });
    catch
        AP_cellrasterJF({stimOn_times}, ...
            {stimIDs});
    end
end

dataGrating(5).spike_templates = spike_templates;
dataGrating(5).spike_times_timeline = spike_times_timeline;
dataGrating(5).stimOn = stimOn_times;
dataGrating(5).stimCond = trial_conditions;
dataGrating(5).stimIDs = stimIDs;
dataGrating(5).depth = [3500, 4000];
dataGrating(5).spike_depths = spike_depths;


%% plotting - average response all 
colorsT = {'r', 'k', 'b', 'g', 'm'}; 
figure();
thisWindow = [-0.2, 0.5];
psthBinSize=0.01;
for i=1:5
    theseSpikes =  dataGrating(i).spike_times_timeline(dataGrating(i).spike_depths > dataGrating(i).depth(1) & dataGrating(i).spike_depths < dataGrating(i).depth(2));
    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, dataGrating(i).stimOn, thisWindow, psthBinSize);
    [psthbaseline, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, dataGrating(i).stimOn, [-0.2,0], psthBinSize);
    psth = (psth - nanmean(psthbaseline)) / nanstd(psthbaseline);
    AP_errorfillJF([-0.2:0.01:0.5-0.01]',psth',std((binnedArray- nanmean(psthbaseline)) / nanstd(psthbaseline))*4,colorsT{i})
    makepretty;
    hold on;
end
ylabel('dFR/FR')
xlabel('time from Grating (s)')
set(gcf,'color','w');

%% plotting - response to orientation each 
colorsO = [rgb('DarkRed'); rgb('Red'); rgb('Orange'); rgb('Brown'); rgb('Purple');rgb('DarkBlue');rgb('Blue');rgb('SkyBlue')];
figure();
thisWindow = [-0.2, 0.5];
psthBinSize=0.01;
orien=unique(dataGrating(5).stimCond(:,3));
for i=1:5
subplot(5,1,i)
if i == 5
for iOr = 1:8
    ff=find(dataGrating(i).stimCond(:,3)==orien(iOr));
    theseSpikes =  dataGrating(i).spike_times_timeline(dataGrating(i).spike_depths > dataGrating(i).depth(1) & dataGrating(i).spike_depths < dataGrating(i).depth(2));
    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, dataGrating(i).stimOn(ff), thisWindow, psthBinSize);

    

    [psthbaseline, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, dataGrating(i).stimOn(ff), [-0.2,0], psthBinSize);
    psth = (psth - nanmean(psthbaseline)) / nanstd(psthbaseline);
    AP_errorfillJF([-0.2:0.01:0.5-0.01]',psth',std((binnedArray- nanmean(psthbaseline)) / nanstd(psthbaseline))*4,colorsO(iOr,:))
    makepretty;
    hold on;
end
else
for iOr = [1,3,5]
    ff=find(dataGrating(i).stimCond(:,3)==orien(iOr));
    theseSpikes =  dataGrating(i).spike_times_timeline(dataGrating(i).spike_depths > dataGrating(i).depth(1) & dataGrating(i).spike_depths < dataGrating(i).depth(2));
    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, dataGrating(i).stimOn(ff), thisWindow, psthBinSize);

    

    [psthbaseline, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, dataGrating(i).stimOn(ff), [-0.2,0], psthBinSize);
    psth = (psth - nanmean(psthbaseline)) / nanstd(psthbaseline);
    AP_errorfillJF([-0.2:0.01:0.5-0.01]',psth',std((binnedArray- nanmean(psthbaseline)) / nanstd(psthbaseline))*4,colorsO(iOr,:))
    makepretty;
    hold on;
end
end

end
ylabel('dFR/FR')
xlabel('time from Grating (s)')
set(gcf,'color','w');
%% plotting - response to sp. frequency each 
colorsO = [rgb('DarkRed'); rgb('Red'); rgb('Orange'); rgb('Brown'); rgb('Purple');rgb('DarkBlue');rgb('Blue');rgb('SkyBlue')];
figure();
thisWindow = [-0.2, 0.5];
psthBinSize=0.01;
orien=unique(dataGrating(5).stimCond(:,2));
for i=1:5
subplot(5,1,i)
if i == 5
for iOr = 1:4
    ff=find(dataGrating(i).stimCond(:,2)==orien(iOr));
    theseSpikes =  dataGrating(i).spike_times_timeline(dataGrating(i).spike_depths > dataGrating(i).depth(1) & dataGrating(i).spike_depths < dataGrating(i).depth(2));
    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, dataGrating(i).stimOn(ff), thisWindow, psthBinSize);

    

    [psthbaseline, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, dataGrating(i).stimOn(ff), [-0.2,0], psthBinSize);
    psth = (psth - nanmean(psthbaseline)) / nanstd(psthbaseline);
    AP_errorfillJF([-0.2:0.01:0.5-0.01]',psth',std((binnedArray- nanmean(psthbaseline)) / nanstd(psthbaseline))*4,colorsO(iOr,:))
    makepretty;
    hold on;
end
else
for iOr = 2:4
if iOr==4
uu=unique(dataGrating(i).stimCond(:,2));
ff=find(dataGrating(i).stimCond(:,2)==uu(end));
    
else
    ff=find(dataGrating(i).stimCond(:,2)==orien(iOr));
end
    theseSpikes =  dataGrating(i).spike_times_timeline(dataGrating(i).spike_depths > dataGrating(i).depth(1) & dataGrating(i).spike_depths < dataGrating(i).depth(2));
    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, dataGrating(i).stimOn(ff), thisWindow, psthBinSize);

    

    [psthbaseline, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, dataGrating(i).stimOn(ff), [-0.2,0], psthBinSize);
    psth = (psth - nanmean(psthbaseline)) / nanstd(psthbaseline);
    AP_errorfillJF([-0.2:0.01:0.5-0.01]',psth',std((binnedArray- nanmean(psthbaseline)) / nanstd(psthbaseline)).*4,colorsO(iOr,:))
    makepretty;
    hold on;
end
end

end
ylabel('dFR/FR')
xlabel('time from Grating (s)')
set(gcf,'color','w');