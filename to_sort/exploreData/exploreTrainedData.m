%% exploring trained data 
% no JF070 striatum ?? 
% STR: JF070 21 (bottom after 3200 - posterior striatum ++) ? ; bottom 31 after 3600 top  1 1
% (until 1500)- after maybe went through white matter trackt?  ? 
% JF071 (/!\ azimuth passive) 
% JF078
% JF067 (/!\ azimuth passive)
% GPe : JF070 day 2,3,4 site 1 (?) 
% SNr : JF070

%% striatum, chaque bonne cellule, baseline go/no go action 

%% load all data 
allDataStruct = JF_loadAllData({'JF070'}, 'CP', 'stage', 1);
raster_window = [-0.5, 2];
psth_bin_size = 0.01;
allData_singleCellPSTH = JF_singleCellPSTH_allData(allDataStruct,'stim_id', 'stim', raster_window, psth_bin_size);


%% baseline in go action vs no go action 
% plot baseline in go vs no go trials depending on stim type 
baseline_goStim_goAction = nanmean(allData_singleCellPSTH.psth(:,1,56:71),3)./1.5;
baseline_goStim_noGoAction = nanmean(allData_singleCellPSTH.psth(:,4,56:71),3)./1.5;
baseline_noGoStim_goAction = nanmean(allData_singleCellPSTH.psth(:,3,56:71),3)./1.5;
baseline_noGoStim_noGoAction = nanmean(allData_singleCellPSTH.psth(:,5,56:71),3)./1.5;
figure('Color', 'white');
scatterHistDiff(baseline_goStim_goAction(baseline_goStim_goAction<=30), baseline_goStim_noGoAction(baseline_goStim_noGoAction<=30), [], [], 'blue',1)
xlim([0 30])
ylim([0 30])
xlabel(['baseline firing rate in' newline ...
    'go stimulus, go action trials'])
ylabel(['baseline firing rate in' newline ...
    'go stimulus, no go action trials'])

figure('Color', 'white');
scatterHistDiff(baseline_noGoStim_goAction(baseline_noGoStim_goAction<=30), baseline_noGoStim_noGoAction(baseline_noGoStim_noGoAction<=30), [], [], 'blue',1)
xlim([0 30])
ylim([0 30])
xlabel(['baseline firing rate in' newline ...
    'no go stimulus, go action trials'])
ylabel(['baseline firing rate in' newline ...
    'no go stimulus, no go action trials'])

JF_singleCellPSTH_flip_allData(allDataStruct,'stim_id', 'stim', raster_window, psth_bin_size);% db mode, run cells one by one currently - need to implement  flipping 
% single cells: inst FR IN BASELINE vs inst Go rate 
% recording 1, unit 2, 15, 56, 60
for iRecording = 1%:size(allDataStruct,2)
    uniqueUnits = unique(allDataStruct(iRecording).spike_templates);
    for iUnit = 1:length(uniqueUnits)
        
        %if strcmp(trialGroups, 'stim_id')
            align_group = allDataStruct(iRecording).trial_conditions(:,1) + 10*(2+allDataStruct(iRecording).trial_conditions(:,2));
      %  end
       % if strcmp(alignTo, 'stim')
            align_times = allDataStruct(iRecording).stimOn_times;
        %end

        [~ ,~, raster_x, raster_y, ~] = JF_raster_PSTH(allDataStruct(iRecording).spike_templates,...
            allDataStruct(iRecording).spike_times, ...
            uniqueUnits(iUnit), raster_window, psth_bin_size, align_times, [], [],[], 0, 1);
        for iTrial = 1:size(align_times,1)
        instBaseline_FR(iUnit, iTrial) = sum(raster_x(raster_y==iTrial)<=50 )*2; %qq hard coded
        end
        %allData_singleCellPSTH.psth(unitCount,:,:) = curr_smoothed_psth;

    end
end
iRecording =1;
[instHit_rate, instCR_rate, instGo_rate ] = JF_getBehavArousalMeasures(allDataStruct(iRecording).trial_conditions(:,1),...
    allDataStruct(iRecording).trial_conditions(:,2));

figure();
subplot(411)
yyaxis left;
plot(align_times, smoothdata(squeeze(nanmean(instBaseline_FR(:, :))),'movmedian', [1,20])); hold on;
yyaxis right;
plot(align_times, instGo_rate)
legend({'inst baseline FR population', 'inst Go rate'})
xlabel('time(s)')
xlim([0 2200])


figure();
subplot(411)
yyaxis left;
plot(align_times, smoothdata(squeeze(instBaseline_FR(2, :)),'movmedian', [1,20])); hold on;
yyaxis right;
plot(align_times, instGo_rate)
legend({'inst baseline FR unit', 'inst Go rate'})
xlabel('time(s)')
xlim([0 2200])

subplot(412)
yyaxis left;
plot(align_times, smoothdata(squeeze(instBaseline_FR(15, :)),'movmedian', [1,20])); hold on;
yyaxis right;
plot(align_times, instGo_rate)
xlim([0 2200])

subplot(413)
yyaxis left;
plot(align_times, smoothdata(squeeze(instBaseline_FR(56, :)),'movmedian', [1,20])); hold on;
yyaxis right;
plot(align_times, instGo_rate)
xlim([0 2200])

subplot(414)
yyaxis left;
plot(align_times, smoothdata(squeeze(instBaseline_FR(60, :)),'movmedian', [1,20])); hold on;
yyaxis right;
plot(align_times, instGo_rate)
xlim([0 2200])
% -> FSIs. show example
% -> MSNs. show example 

% single trial population raster 

% correlation no go rate with pupil/whisking 

% example cell 
[curr_smoothed_psth, curr_psth, raster_x, raster_y, curr_raster] = JF_raster_PSTH(spike_templates, spike_times_timeline, ...
    thisTemplate, raster_window, psth_bin_size, align_times, align_group, sort_by, color_by, plot_me, causal_smoothing);

JF_singleTrialPSTH_allData

%% raster map 
ops=struct;
clus = unique(allDataStruct(1).spike_templates);
ops.nC = 2;%length(clus);%, number of clusters to use 
ops.iPC = 1:length(clus);%, number of PCs to use 
ops.isort = [];%, initial sorting, otherwise will be the top PC sort
ops.useGPU = 0;%, whether to use the GPU
ops.upsamp = 100;%, upsampling factor for the embedding position
ops.sigUp = 1;%, % standard deviation for upsampling
[cell_traces, isort1, isort2, timeVector,removeMe ] = JF_runRasterMap(allDataStruct(1).spike_times, allDataStruct(1).spike_templates, ops);

figure();
subplot(5,1,[1:3])
imagesc(timeVector, [], cell_traces(isort1,:)); hold on;
xlabel('time (s)')
ylabel('neuron #')
makepretty;

subplot(5,1,4)
plot(Timeline.rawDAQTimestamps, smoothdata(abs(wheel_velocity),'movmean', 10000))
xlim([0 2000])
ylabel('wheel velocity')
makepretty;

subplot(5,1,5)
iRecording =1;
[instHit_rate, instCR_rate, instGo_rate ] = JF_getBehavArousalMeasures(allDataStruct(iRecording).trial_conditions(:,1),...
    allDataStruct(iRecording).trial_conditions(:,2));


plot(stimOn_times, smoothdata(instCR_rate, 'movmean', 10)); hold on;
plot(stimOn_times, smoothdata(instHit_rate, 'movmean', 10))
xlim([0 2000])
legend({'CR rate', 'Hit rate'})
makepretty;

figure();
imagesc(timeVector, [], cell_traces(:,:)); hold on;
xlabel('time (s)')
ylabel('neuron #')
makepretty;
                    
%     line([stimOn_times(iStim), stimOn_times(iStim)], [0.5, size(cell_traces,1)])
% end
 %plot movement at bottom 

%overlay stimOn times, rates ect 
% 
% 
% figure();
% imagesc(timeVector, [], cell_traces(isort1,isort2)); hold on;
% xlabel('time (s)')
% ylabel('neuron #')
% makepretty;

% single trial population raster (neuron x time for one trial: -0.5 to 2
% seconds
figure();

iTrial = iTrial+1;
all_templates=unique(allDataStruct(1).spike_templates);
for iUnit = 1:length(all_templates)
timeVector = stimOn_times(iTrial)-0.5:0.01:stimOn_times(iTrial)+2;
curr_raster_spike_times = allDataStruct(1).spike_times(allDataStruct(1).spike_templates == all_templates(iUnit));
curr_raster_spike_times(curr_raster_spike_times > timeVector(end))= [];
pop_raster(iUnit,:) = histcounts(curr_raster_spike_times, timeVector);

end
pop_raster(removeMe,:)=[];
clf;
imagesc(-0.5:0.01:2, [], pop_raster(isort1, :)); hold on;
colormap(brewermap([],'Greys'));
line([0, 0], [0.5, length(isort1)+0.5], 'Color', 'r');
line([stim_to_move(iTrial), stim_to_move(iTrial)], [0.5, length(isort1)+0.5], 'Color',rgb('Purple'));
line([stim_to_feedback(iTrial), stim_to_feedback(iTrial)], [0.5, length(isort1)+0.5], 'Color',rgb('Blue'));
title(['stim = ', num2str(stimIDs(iTrial))])