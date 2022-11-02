
%% ~~tracking cells over days playground~~
% dependancies: 
% - CCGBz 
% - makepretty

% load dataset from Celian:
% It contains the matched clusters IDs for 2 different days (d1 and d2). It should be hopefully straightforward to interpret:
% dat.d1 contains everything about day 1: the ephys raw and sorting folders & file,
% and in dat.d1.clu the cluster's IDs, depths and x-locations. Same in dat.d2.
% dat.pairCorr gives you an idea of "how well" these clusters were matched (that's a rough idea).
% clusters removed: noise, un-reliable response, no good match

%% load dataset info
load('/home/netshare/zserver-lab/Share/Celian/dataForJulie.mat')
d1_folder = strrep(dat.d1.ephysSortingFolder, '\', filesep); % windows -> linux
d2_folder = strrep(dat.d2.ephysSortingFolder, '\', filesep); % windows -> linux
d1_folder = strrep(d1_folder, '/128.40.224.65/Subjects/', 'home/netshare/tempserver/'); % folder location
d2_folder = strrep(d2_folder, '/128.40.224.65/Subjects/', 'home/netshare/tempserver/'); % folder location
d1_rawfolder = strrep(dat.d1.ephysFolder, '\', filesep); % windows -> linux
d2_rawfolder = strrep(dat.d2.ephysFolder, '\', filesep); % windows -> linux
d1_rawfolder = strrep(d1_rawfolder, '/128.40.224.65/Subjects/', 'home/netshare/tempserver/'); % folder location
d2_rawfolder = strrep(d2_rawfolder, '/128.40.224.65/Subjects/', 'home/netshare/tempserver/'); % folder location
dat.d1.clu.ID = dat.d1.clu.ID + 1; % 0-idx -> 1-idx
dat.d2.clu.ID = dat.d2.clu.ID + 1; % 0-idx -> 1-idx
foldersAll = {d1_folder; d2_folder};

%% load experimental data
for iDay = 1:size(foldersAll, 1)
    
    templates{iDay} = readNPY([foldersAll{iDay, 1}, filesep, 'templates.npy']);
    channel_positions{iDay} = readNPY([foldersAll{iDay, 1}, filesep, 'channel_positions.npy']);
    spike_times{iDay} = double(readNPY([foldersAll{iDay, 1}, filesep, 'spike_times.npy']))./30000; % sample rate hard-coded as 30000 - should load this in from params 
    spike_templates{iDay} = readNPY([foldersAll{iDay, 1}, filesep, 'spike_templates.npy']) + 1; % 0-idx -> 1-idx
    template_amplitude{iDay} = readNPY([foldersAll{iDay, 1}, filesep, 'amplitudes.npy']);
    spike_clusters{iDay} = readNPY([foldersAll{iDay, 1}, filesep, 'spike_clusters.npy']) + 1;
    %mean_template_ampltiude
end

%% re-do location plot (ground truth)
for iDay = 1:size(foldersAll, 1)
    
    % get unique clusters
    if iDay == 1
        unique_clu = dat.d1.clu.ID;
    else
        unique_clu = dat.d2.clu.ID;
    end
    n_clu = size(unique_clu, 1);
    
    % Get the waveform of all templates (channel with largest amplitude)
    [~, max_site] = max(max(abs(templates{iDay}(unique_clu, :, :)), [], 2), [], 3);
    templates_max = nan(n_clu, size(templates{iDay}, 2));
    
    for curr_template = 1:n_clu
        this_template = unique_clu(curr_template);
        templates_max(curr_template, :) = ...
            templates{iDay}(this_template, :, max_site(curr_template));
        template_positions(curr_template, :, iDay) = channel_positions{iDay}(max_site(curr_template), :); % template position 
    end
    waveforms{iDay} = templates_max; %template waveform on max channel

end

% get euclidean distance
euclidian_dist = abs(sum(abs(template_positions(:, :, 1)-template_positions(:, :, 2)), 2));

% plot
figure()
subplot(231)
scatter(1:n_clu, euclidian_dist, 'filled')
xlabel('cluster #')
ylabel({'euclidean distance', 'between best match'})
makepretty;
set(gcf, 'color', 'white')

% ear-mark units with  drift <= 100um as potentially good 
keep_clusters = find(euclidian_dist <= 100);
disp(['keeping ', num2str(length(keep_clusters)), ' out of ', num2str(n_clu+1), ' clusters = ', ...
    num2str(round(length(keep_clusters)/(n_clu + 1)*100)), ' % of good clusters']) 


subplot(232)
scatter(template_positions(:, 1, 1), euclidian_dist, 'filled')
xlabel('cluster shank position day 1')

makepretty;

subplot(233)
scatter(template_positions(:, 2, 1), euclidian_dist, 'filled')
xlabel('cluster depth day 1')

makepretty;

subplot(234)
spiking_stat_window = max(spike_times{1}) - min(spike_times{1});
spiking_stat_bins = [min(spike_times{1}), max(spike_times{1})]; % not uysing this right now, but change spiking stat bins (increase the number of bins)
% to exclude epochs where unit has drifted/there are no spikes 

% Get firing rate across the session
unique_clu = dat.d1.clu.ID;
n_clu = size(unique_clu, 1);
bin_spikes = nan(n_clu, ...
    length(spiking_stat_bins)-1);
for curr_template = 1:n_clu
    this_template = unique_clu(curr_template);
    bin_spikes(curr_template, :) = ...
        histcounts(spike_times{1}(spike_templates{1} == this_template), ...
        spiking_stat_bins);
end
min_spikes = 10;
use_spiking_stat_bins = bsxfun(@ge, bin_spikes, prctile(bin_spikes, 80, 2)) & bin_spikes > min_spikes; 
spike_rate = sum(bin_spikes, 2) ./ ...
    ( double(spiking_stat_window));

scatter(spike_rate, euclidian_dist, 'filled')
xlabel('spike rate day 1')
makepretty;

% amplitude (usually a good proxy for unit quality)
% subplot(235)
% scatter(nanmean(template_amplitude{iDay}(spike_template==unique_clu)), euclidian_dist, 'filled')
% xlabel('amplitude day 1')
% makepretty;

%Celian's pairCorr
subplot(236)
scatter(dat.pairCorr, euclidian_dist, 'filled')
xlabel('match index')
makepretty;



%% get ACGs and waveforms 

% examples, flip through -> how similar are they?
iCluster = 1;

iCluster=iCluster+1;
unique_clu = dat.d1.clu.ID;
theseSpikes1 = spike_times{1}(spike_templates{1} == unique_clu(iCluster)); 
thisWaveform1 = waveforms{1}(iCluster,:);

unique_clu = dat.d2.clu.ID;
theseSpikes2 = spike_times{2}(spike_templates{2} == unique_clu(iCluster)); 
thisWaveform2 = waveforms{2}(iCluster,:);

[ccg1, t] = CCGBz([double(theseSpikes1); double(theseSpikes1)], [ones(size(theseSpikes1, 1), 1); ...
        ones(size(theseSpikes1, 1), 1) * 2], 'binSize', 0.001, 'duration', 0.5, 'norm', 'rate');
[ccg, t] = CCGBz([double(theseSpikes2); double(theseSpikes2)], [ones(size(theseSpikes2, 1), 1); ...
        ones(size(theseSpikes2, 1), 1) * 2], 'binSize', 0.001, 'duration', 0.5, 'norm', 'rate');
    
wv_duration1 = find(thisWaveform1 == max(thisWaveform1 )) - find(thisWaveform1 == min(thisWaveform1 ));
wv_duration2 = find(thisWaveform2 == max(thisWaveform2 )) - find(thisWaveform2 == min(thisWaveform2 ));
acg_peak_end_ratio1 = max(smoothdata(ccg1(:, 1, 1), 'gaussian', 40))/nanmean(smoothdata(ccg1(450:500, 1, 1)));
%figure();plot(smoothdata(ccg1(:, 1, 1), 'gaussian', 40))
acg_peak_end_ratio2 = max(smoothdata(ccg(:, 1, 1), 'gaussian', 40))/nanmean(smoothdata(ccg(450:500, 1, 1)));
acg_pss1 = find(smoothdata(ccg1(:, 1, 1), 'gaussian', 40) >= nanmean(smoothdata(ccg1(450:500, 1, 1))),1,'first');
acg_pss2 = find(smoothdata(ccg(:, 1, 1), 'gaussian', 40) >= nanmean(smoothdata(ccg(450:500, 1, 1))),1,'first');
R2 = corrcoef(thisWaveform1 ,thisWaveform2);
R1 = corrcoef(smoothdata(ccg1(:, 1, 1), 'gaussian', 40) ,smoothdata(ccg(:, 1, 1), 'gaussian', 40));


clf; 
ccg1(t==0,1,1)=0; 
sgtitle(['pairCorr = ',num2str(dat.pairCorr(iCluster)), ...
    ' eucl dist = ' num2str(euclidian_dist(iCluster))]); hold on;
subplot(221)
title(['acg-rt = ' num2str(acg_peak_end_ratio1) ' acg-pss = ' num2str(acg_pss1) ' r = ' num2str(R1(2))]); hold on;
area(1:size(ccg1,1),ccg1(:, 1, 1));
xlim([round(size(ccg1,1)/2) size(ccg1,1)])

subplot(222)
title(['wf-dur = ' num2str(wv_duration1)  ' r = ' num2str(R2(2))]); hold on;
plot(thisWaveform1);

subplot(223)
ccg(t==0,1,1)=0;
title(['acg-rt = ' num2str(acg_peak_end_ratio2) ' acg-pss = ' num2str(acg_pss2)]); hold on;
area(1:size(ccg,1),ccg(:, 1, 1));
xlim([round(size(ccg,1)/2) size(ccg,1)])
disp(iCluster)

subplot(224)
title(['wf-dur = ' num2str(wv_duration1) ]); hold on;
plot(thisWaveform2);

%% ~~method 2: first use quality metrics to subselect units~~
qMetric = struct; 
for iUnit = 1:length(dat.d1.clu.ID)
%n spikes
theseSpikes = spike_times{iDay}(spike_templates{iDay} == iUnit);
theseAmplis = template_amplitude{iDay}(spike_templates{iDay} == iUnit);
if ~isempty(theseSpikes)
    qMetric.numSpikes(iUnit) = numel(theseSpikes);

    % % spikes missing in xx chunks
    %super hacky - you have to cd to directory where gaussFit.py lives for
    %this to work. - change to where your file lives on your computer
%     cd('/home/julie/Dropbox/MATLAB/onPaths/bombcell/helpers/qualityMetricHelpers')
%     try
%         [percent_missing_ndtrAll, ~] = ampli_fit_prc_missJF(theseAmplis, 0);
%     catch
%         percent_missing_ndtrAll = NaN;
%     end
%     qMetric.pMissing(iUnit) = percent_missing_ndtrAll;
     [p,mu,stdev,n,x] = gaussian_mssing(theseAmplis);
     qMetric.pMissing(iUnit) = p;
        waveformsTemp_mean = waveforms{1}(iUnit,:);
    minProminence = 0.2 * max(abs(squeeze(waveformsTemp_mean)));

    %figure();plot(qMetric.waveform(iUnit, :))
    [PKS, LOCS] = findpeaks(squeeze(waveformsTemp_mean), 'MinPeakProminence', minProminence);
    [TRS, LOCST] = findpeaks(squeeze(waveformsTemp_mean)*-1, 'MinPeakProminence', minProminence);
    if isempty(TRS)
        TRS = min(squeeze(waveformsTemp_mean));
        if numel(TRS) > 1
            TRS = TRS(1);
        end
        LOCST = find(squeeze(waveformsTemp_mean) == TRS);
    end
    if isempty(PKS)
        PKS = max(squeeze(waveformsTemp_mean));
        if numel(PKS) > 1
            PKS = PKS(1);
        end
        LOCS = find(squeeze(waveformsTemp_mean) == PKS);
    end

    peakLoc = LOCS(PKS == max(PKS));
    if numel(peakLoc) > 1
        peakLoc = peakLoc(1);

    end
    troughLoc = LOCST(TRS == max(TRS));
    if numel(troughLoc) > 1
        troughLoc = troughLoc(1);
    end

    if peakLoc > troughLoc
        qMetric.somatic(iUnit) = 1;
    else
        qMetric.somatic(iUnit) = 0;
    end

    %false positives
    [qMetric.fractionRPVchunk(iUnit), qMetric.numRPVchunk(iUnit)] = fractionRPviolationsJF( ...
        numel(theseSpikes), theseSpikes, 0.0010, 0.0002,  max(theseSpikes)- min(theseSpikes)); %method from Hill et al., 2011
end
    
end
keepUnits = qMetric.numSpikes > 300 & qMetric.somatic & qMetric.fractionRPVchunk <= 2 & qMetric.pMissing < 0.000001;
%% ~~method 1: find closest ACG and waveform match for each cell ~~
keepUnits = ones(numel(dat.d1.clu.ID),1);
% get ACGs on day 1+2
for iDay = 1:size(foldersAll, 1)
    % n clusters
    if iDay == 1
        unique_clu = dat.d1.clu.ID(keepUnits);
    else
        unique_clu = dat.d2.clu.ID(keepUnits);
    end
    n_clu = size(unique_clu, 1);
    % Get the waveform of all templates (channel with largest amplitude)
    for curr_template = 1:n_clu
       theseSpikes = spike_times{iDay}(spike_templates{iDay} == unique_clu(curr_template)); 
           [ccg, t] = CCGBz([double(theseSpikes); double(theseSpikes)], [ones(size(theseSpikes, 1), 1); ...
        ones(size(theseSpikes, 1), 1) * 2], 'binSize', 0.001, 'duration', 0.5, 'norm', 'rate');
    ccg(t==0,1,1)=0;
    all_ccg{iDay}(curr_template,:) = ccg(:,1,1);
    end

end

% correlations
clearvars all_corrs

unique_clu = dat.d1.clu.ID(keepUnits);%subset the @good@ ones
n_clu = size(unique_clu, 1);
% Get the waveform of all templates (channel with largest amplitude)
for curr_template = 1:n_clu
    unique_clu2 = dat.d2.clu.ID(keepUnits);
    n_clu2 = size(unique_clu, 1);
   for curr_template2 = 1:n_clu2
       
        smoothedACG1 = smoothdata(all_ccg{1}(curr_template, :), 'gaussian', 40);
        smoothedACG2 = smoothdata(all_ccg{2}(curr_template2, :), 'gaussian', 40);
        ccc = corrcoef(smoothedACG1, smoothedACG2);

        all_corrs(curr_template, curr_template2, 1) = ccc(2);
        % normalize waveforms to reduce effect of cell bring closer/further
        % away from probe
        normWF1 = (waveforms{1}(curr_template,25:end) - min(waveforms{1}(curr_template,25:end))) ./...
                (max(waveforms{1}(curr_template,25:end)) - min(waveforms{1}(curr_template,25:end)));
        normWF2 = (waveforms{2}(curr_template2,25:end) - min(waveforms{2}(curr_template2,25:end))) ./...
                (max(waveforms{2}(curr_template2,25:end)) - min(waveforms{2}(curr_template2,25:end)));
        
        cc = max(xcorr(normWF1 , ...
         normWF2 )); % take max of xcorr is simplest way to ensure there's no artefact of waveforms not being properly aligned 
        all_corrs(curr_template, curr_template2, 2) = cc;


   end

end

% super inefficent way of doing things 
clearvars single_best_match

single_best_match = zeros(size(all_corrs,1),1);


excludeThis = [];
all_corrs_excluded = sum(all_corrs,3);

 for iMax = 1:size(all_corrs,1)
     %iMax=iMax+1
     lin_corrs = reshape(all_corrs_excluded,1,[]);
     lin_corrs = lin_corrs(~isnan(lin_corrs));
     [~,sorted_corrs] = sort(lin_corrs);
     theseMaxCorrs = lin_corrs(sorted_corrs);
     thisMaxCorr = theseMaxCorrs(end);
    
     [thisMaxX, thisMaxY] = find(all_corrs_excluded == thisMaxCorr);
    
     if length(thisMaxX)>1
         disp('warning: 2 identical values')
     end
     single_best_match(thisMaxX) = thisMaxY;
     disp(single_best_match(thisMaxX))
     all_corrs_excluded(:,thisMaxY) = NaN;% ensure you more than one n day unit can't be matched with the same n+1 day unit
     all_corrs_excluded(thisMaxX,:) = NaN;% ensure you more than one n+1 day unit can't be matched with the same n day unit
 end

%% plot ACG/WF matching vs functional matching 
figure()
subplot(231)
keepUnits(single_best_match==0)=0;
single_best_match(single_best_match==0)=[];

euclidian_dist2 = abs(sum(template_positions(logical(keepUnits), :, 1)-template_positions(single_best_match, :, 2), 2));
scatter(find(keepUnits), euclidian_dist(logical(keepUnits)), 'filled')
hold on;
scatter(1:145, euclidian_dist, 'filled')
hold on;

scatter(find(keepUnits), euclidian_dist2, 'filled')
xlabel('cluster #')
ylabel({'euclidean distance', 'between best match'})
makepretty;
set(gcf, 'color', 'white')


subplot(232)
cla;
scatter(template_positions(:, 1, 1), euclidian_dist, 'filled')
hold on;
scatter(template_positions(logical(keepUnits), 1, 1), euclidian_dist2, 'filled')
xlabel('cluster shank position day 1')

makepretty;

subplot(233)
cla;
scatter(template_positions(:, 2, 1), euclidian_dist, 'filled')
hold on;
scatter(template_positions(logical(keepUnits), 2, 1), euclidian_dist2, 'filled')
xlabel('cluster depth day 1')

makepretty;

subplot(234)
cla;
spiking_stat_window = max(spike_times{1}) - min(spike_times{1});
spiking_stat_bins = [min(spike_times{1}), max(spike_times{1})]; % not uysing this right now, but change spiking stat bins (increase the number of bins)
% to exclude epochs where unit has drifted/there are no spikes 

% Get firing rate across the session
unique_clu = dat.d1.clu.ID;
n_clu = size(unique_clu, 1);
bin_spikes = nan(n_clu, ...
    length(spiking_stat_bins)-1);
for curr_template = 1:n_clu
    this_template = unique_clu(curr_template);
    bin_spikes(curr_template, :) = ...
        histcounts(spike_times{1}(spike_templates{1} == this_template), ...
        spiking_stat_bins);
end
min_spikes = 10;
use_spiking_stat_bins = bsxfun(@ge, bin_spikes, prctile(bin_spikes, 80, 2)) & bin_spikes > min_spikes; 
spike_rate = sum(bin_spikes, 2) ./ ...
    ( double(spiking_stat_window));

scatter(spike_rate, euclidian_dist, 'filled')
hold on;
scatter(spike_rate(logical(keepUnits)), euclidian_dist2, 'filled')
xlabel('spike rate day 1')
makepretty;

% amplitude (usually a good proxy for unit quality)
subplot(235)
scatter(template_amplitude{iDay}(unique_clu), euclidian_dist, 'filled');
hold on;
scatter(template_amplitude{iDay}(unique_clu(logical(keepUnits))), euclidian_dist2, 'filled');
xlabel('amplitude day 1')
makepretty;

%Celian's pairCorr
subplot(236)
scatter(dat.pairCorr, euclidian_dist, 'filled')
xlabel('match index')
makepretty;
%% 
n_correct  = sum(single_best_match' == 1:length(single_best_match));
euclidian_dist2 = abs(sum(template_positions(:, :, 1)-template_positions(single_best_match, :, 2), 2));

sum(euclidian_dist2 < 20)
sum(euclidian_dist1 < 20)
figure(); 
clf;
scatter(1:145, euclidian_dist1);
hold on; 
scatter(1:145, euclidian_dist2)% let's see examples of (false) best match - whjat is going wrong? 
curr_template = curr_template+1

figure(5)
clf;
subplot(121)

smoothedACG1= smoothdata(all_ccg{1}(curr_template, 250:end), 'gaussian', 40);
smoothedACG1 = (smoothedACG1- min(smoothedACG1)) ./ (max(smoothedACG1) - min(smoothedACG1));
plot(smoothedACG1)

hold on; 
smoothedACG2= smoothdata(all_ccg{2}(curr_template, 250:end), 'gaussian', 40);
smoothedACG2 = (smoothedACG2- min(smoothedACG2)) ./ (max(smoothedACG2) - min(smoothedACG2));

plot(smoothedACG2)
ss3= smoothdata(all_ccg{2}(single_best_match(curr_template), 250:end), 'gaussian', 40);
ss3 = (ss3- min(ss3)) ./ (max(ss3) - min(ss3));
plot(ss3)
legend({'day n unit', 'celian match' 'acg/wv match'})
makepretty; 
corrcoef(smoothedACG1,smoothedACG2)
corrcoef(smoothedACG1,ss3)
subplot(122)
plot((waveforms{1}(curr_template,25:end) - min(waveforms{1}(curr_template,25:end))) ./...
    (max(waveforms{1}(curr_template,25:end)) - min(waveforms{1}(curr_template,25:end))))

hold on;
plot((waveforms{2}(curr_template,25:end) - min(waveforms{2}(curr_template,25:end))) ./...
    (max(waveforms{2}(curr_template,25:end)) - min(waveforms{2}(curr_template,25:end))))

hold on; 
plot((waveforms{2}(single_best_match(curr_template),25:end) - min(waveforms{2}(single_best_match(curr_template),25:end))) ./...
    (max(waveforms{2}(single_best_match(curr_template),25:end)) - min(waveforms{2}(single_best_match(curr_template),25:end))))

max(xcorr(waveforms{1}(curr_template,25:end),waveforms{2}(curr_template,25:end)))
max(xcorr(waveforms{1}(curr_template,25:end),waveforms{2}(single_best_match(curr_template),25:end)))
legend({'', ['eucl. distance functional match = ' num2str(euclidian_dist1(curr_template))],...
    ['eucl. distance ACG/WV match = ' num2str(euclidian_dist2(curr_template))]})
makepretty;




