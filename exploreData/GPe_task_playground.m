% GPe task playground 
% ex days:
% GPe: JF070 21 
% STR: JF071 11
% SNr: JF070 32(?), JF071 32/3 
%% ~ load data ~
animals={'JF070'};
curr_animal = 1; % (set which animal to use)
protocol = []; % (this is the name of the Signals protocol)
curr_day = 2; % (set which day to use)
load_parts.cam=true;
load_parts.imaging=false;
load_parts.ephys=true;
site = 1;%1,1; 2,4; 3,7
recording = []; 
experiment = 2;
loadClusters = 0;

loadEphys; 
curr_shank = NaN;
AP_cellrasterJF({stimOn_times,wheel_move_time,signals_events.responseTimes(n_trials(1):n_trials(end))',stimOn_times}, ...
    {trial_conditions(:,1),trial_conditions(:,2), ...
    trial_conditions(:,3), movement_after200ms_and_type});

%% ~ 1. periods low baseline activity ~
%% ~~ drift artefact? ~~
theseUnits = [92, 86, 87, 82];
figure(); 
for iUnit = 1:length(theseUnits)
    thisUnit = theseUnits(iUnit);
    subplot(7,1,iUnit)
    theseSpikes = spike_templates == thisUnit;
    [raster_y,raster_x, curr_smoothed_psth, t] = ...
        jf_PSTH(spike_templates, thisUnit, spike_times_timeline, stimOn_times, ones(size(stimOn_times,1),1));
    %figure();
    %scatter(t(raster_x), raster_y)
    theseTrials = findgroups(raster_y(t(raster_x)<=0));
    nSpikes_baseline = grpstats(raster_y(t(raster_x)<=0),theseTrials ,@(x) length(x));
    [theseTrials_sorted, theseTrials_sortedIdx] = sort(unique(theseTrials));
    lowBaseline_trials = find(nSpikes_baseline < prctile(nSpikes_baseline(theseTrials_sortedIdx), 30));
    lowBaseline_times = stimOn_times(lowBaseline_trials);
    scatter(spike_times_timeline(theseSpikes),template_amplitudes(theseSpikes), 2, 'filled');hold on;
    for iLow = 1:size(lowBaseline_times,1)
        f= fill([lowBaseline_times(iLow)-0.5,lowBaseline_times(iLow)+2,lowBaseline_times(iLow)+2,lowBaseline_times(iLow)-0.5],...
            [0,0,40,40], 'red');
        f.FaceColor = 'red';
        f.EdgeColor = 'none';
        f.FaceAlpha = 0.5;

    end
    xlim([min(stimOn_times), max(stimOn_times)])
    
end
% -> not sure 

%% ~~ behavior ? ~~
%stim_to_move = wheel_starts(wheel_move_stim_idx) - stimOn_times(1:n_trials);
instHit_rate = movmean(double(~isnan(stim_to_move(ismember(stimIDs(1:length(n_trials)), [1, 2])))), 20);
instCR_rate = movmean(double(isnan(stim_to_move(ismember(stimIDs(1:length(n_trials)), [3])))), 20);
instReactionTime = movmean(double(stim_to_move(ismember(stimIDs(1:length(n_trials)), [1, 2]))), 20, 'omitnan');

subplot(7, 1, 5)
plot(stimOn_times(ismember(stimIDs(1:length(n_trials)), [1, 2])), instHit_rate); hold on;
plot(stimOn_times(ismember(stimIDs(1:length(n_trials)), [3])), instCR_rate); hold on;
plot(stimOn_times(ismember(stimIDs(1:length(n_trials)), [1, 2])), instReactionTime)
legend({'Hit rate', 'Correct No go rate', 'Reaction time'})
xlabel('time (s)')
makepretty;

%% ~~ state? ~~
subplot(7, 1, 6)
% pupil size, whisking motion energy 

% facemap 
face_proc = load('/home/netshare/zinu/JF070/2022-06-10/1/face_proc.mat');
eye_proc = load('/home/netshare/zinu/JF070/2022-06-10/1/eye_proc.mat');
%eye_proc_NPY = readNPY('/home/netshare/zinu/JF070/2022-06-10/1/eye_proc.npy');
startFrame = find(~isnan(eyecam_t), 1, 'first');
stopFrame = size(eyecam_t,1) - find(~isnan(eyecam_t(end:-1:1)), 1, 'first') + 1;
% figure();
% imagesc(face_proc.avgmotion_0)
% colormap(brewermap([],'*RdBu'))

yyaxis left; cla; 
plot(facecam_t, smoothdata(nanmean(face_proc.motSVD_0,2), 'movmean', [1, 1000])); hold on; 
ylabel('whisker motion energy')
makepretty;

yyaxis right; cla; 
plot(eyecam_t(startFrame:stopFrame),smoothdata(pupilArea_pixels(startFrame:end), 'movmedian', [1,1000]))
xlabel('time (s)')
xlim([eyecam_t(startFrame)+10, eyecam_t(stopFrame)-10])
ylabel('pupil size')
legend({'whisker motion energy', 'pupil size'})
makepretty;


subplot(7,1,7)
plot(Timeline.rawDAQTimestamps, movmean(abs(wheel_velocity), 10000)) 
ylabel('wheel velocity')
xlabel('time(s)')
makepretty;

%% ~ 2. difference passive / task ? ~
% is there a pupil size/whisking difference? mou
%% ~ 3. fraction cells respond ~ 

%% ~ 4. over learning ~ 

