
%% compute quality metrics for certain quality metrics
function [qMetric, ephysParams] = qualityMetEphysParamJF(qMetric, ephysParams, ephysData, raw, param, allT, iUnit, str_templates, timeChunks)
%qMetric = struct;
%ephysParams = struct;
thisRawUnit = allT(iUnit);
qMetric.thisRawUnit(iUnit) = thisRawUnit;
if param.strOnly
    strUnits = find(str_templates);
    thisUnit = strUnits(iUnit);
else
    thisUnit = iUnit;
end
qMetric.thisUnit(iUnit) = thisUnit;
%%qq addeParams, add somatic

if ~isempty(timeChunks)
    theseSpikesIdx = ephysData.spike_templates == thisUnit & ephysData.spike_times_timeline >= timeChunks(1) ...
        & ephysData.spike_times_timeline < timeChunks(2);
    theseSpikes = ephysData.spike_times_timeline(theseSpikesIdx);
    theseAmplis = ephysData.template_amplitudes(theseSpikesIdx);
    pc_feature_ind = ephysData.pc_features_ind;

    allSpikesIdx = ephysData.spike_times_timeline >= timeChunks(1) ...
        & ephysData.spike_times_timeline < timeChunks(2);
    %pc_features = ephysData.pc_features(allSpikesIdx, :, :);
else
    theseSpikesIdx = ephysData.spike_templates == thisUnit;
    theseSpikes = ephysData.spike_times_timeline(theseSpikesIdx);
    theseAmplis = ephysData.template_amplitudes(theseSpikesIdx);
    pc_feature_ind = ephysData.pc_features_ind;

    allSpikesIdx = ephysData.spike_templates;
    %pc_features = ephysData.pc_features;
    timeChunks = [min(theseSpikes), max(theseSpikes)];
end

%% QUALITY METRICS

%% num spikes
qMetric.numSpikes(iUnit) = numel(theseSpikes);

if param.raw

    %% get raw data and raw ampli
    % gg = find(ephysData.good_templates);
    curr_template = thisRawUnit; %gg(iUnit);
    ns = find(ephysData.new_spike_idx == curr_template);

    curr_spikes_idx = find(ephysData.spike_templates_full == ns);
    if ~isempty(curr_spikes_idx)
        curr_pull_spikes = unique(round(linspace(1, length(curr_spikes_idx), raw.max_pull_spikes)));
        if curr_pull_spikes(1) == 0
            curr_pull_spikes(1) = [];
        end
        curr_spikeT = ephysData.spike_times_full(curr_spikes_idx(curr_pull_spikes));
        curr_spikeT_pull = double(curr_spikeT) + raw.pull_spikeT;

        out_of_bounds_spikes = any(curr_spikeT_pull < 1, 2) | ...
            any(curr_spikeT_pull > size(raw.ap_data.data.data, 2), 2);
        curr_spikeT_pull(out_of_bounds_spikes, :) = [];

        curr_spike_waveforms = reshape(raw.ap_data.data.data(:, reshape(curr_spikeT_pull', [], 1)), raw.n_channels, length(raw.pull_spikeT), []);
        if ~isempty(curr_spike_waveforms)
            curr_spike_waveforms_car = curr_spike_waveforms - nanmedian(curr_spike_waveforms, 1);
            curr_spike_waveforms_car_sub = curr_spike_waveforms_car - curr_spike_waveforms_car(:, 1, :);

            waveforms_mean(curr_template, :, :) = ...
                permute(nanmean(curr_spike_waveforms_car_sub(raw.used_channels_idx, :, :), 3), [3, 2, 1]) * raw.microVoltscaling;


            thisChannelRaw = find(squeeze(max(waveforms_mean(curr_template, :, :))) == ...
                max(squeeze(max(waveforms_mean(curr_template, :, :)))));
            if numel(thisChannelRaw) > 1
                thisOne = find(max(abs(waveforms_mean(curr_template, :, thisChannelRaw))) == max(max(abs(waveforms_mean(curr_template, :, thisChannelRaw)))));
                if numel(thisOne) > 1
                    thisOne = thisOne(1);
                end
                thisChannelRaw = thisChannelRaw(thisOne);
            end
            qMetric.waveforms_mean(iUnit, :, :) = waveforms_mean(curr_template, :, :);
            qMetric.thisChannelRaw(iUnit) = thisChannelRaw;

            wavefPeak = max(waveforms_mean(curr_template, :, thisChannelRaw));
            if numel(wavefPeak) > 1
                wavefPeak = mean(wavefPeak);
            end
            wavefTrough = abs(min(waveforms_mean(curr_template, :, thisChannelRaw)));
            if numel(wavefTrough) > 1
                wavefTrough = mean(wavefTrough);
            end
            qMetric.waveformRawAmpli(iUnit) = wavefPeak + ...
                wavefTrough;
            wavefPeakLoc = (find(waveforms_mean(curr_template, :, thisChannelRaw) == ...
                max(waveforms_mean(curr_template, :, thisChannelRaw))));
            if numel(wavefPeakLoc) > 1
                wavefPeakLoc = mean(wavefPeakLoc);
            end
            wavefTroughLoc = (find(waveforms_mean(curr_template, :, thisChannelRaw) == ...
                min(waveforms_mean(curr_template, :, thisChannelRaw))));
            if numel(wavefTroughLoc) > 1
                wavefTroughLoc = mean(wavefTroughLoc);
            end
            qMetric.waveformRawDurationUs(iUnit) = (wavefPeakLoc - wavefTroughLoc) / ephysData.ephys_sample_rate * 1e6;
        else
            qMetric.waveformRawAmpli(iUnit) = NaN;
            qMetric.waveformRawDuration(iUnit) = NaN;
        end
    else
        qMetric.waveformRawAmpli(iUnit) = NaN;
        qMetric.waveformRawDuration(iUnit) = NaN;

    end

    %% number of peaks of raw mean waveform

    if ~isempty(curr_spikes_idx) && ~isempty(curr_spike_waveforms)
        %disp(iUnit)
        minProminence = 0.2 * max(abs(squeeze(qMetric.waveforms_mean(iUnit, :, qMetric.thisChannelRaw(iUnit)))));

        [PKS, LOCS] = findpeaks(squeeze(qMetric.waveforms_mean(iUnit, :, qMetric.thisChannelRaw(iUnit))), 'MinPeakProminence', minProminence);
        [TRS, LOCST] = findpeaks(squeeze(qMetric.waveforms_mean(iUnit, :, qMetric.thisChannelRaw(iUnit)))*-1, 'MinPeakProminence', minProminence);
        if isempty(TRS)
            TRS = min(squeeze(waveforms_mean(:, :, qMetric.thisChannelRaw(iUnit))));
            if numel(TRS) > 1
                TRS = TRS(1);
            end
            LOCST = find(squeeze(waveforms_mean(:, :, qMetric.thisChannelRaw(iUnit))) == TRS);
        end
        if isempty(PKS)
            PKS = max(squeeze(waveforms_mean(:, :, qMetric.thisChannelRaw(iUnit))));
            if numel(PKS) > 1
                PKS = PKS(1);
            end
            LOCS = find(squeeze(waveforms_mean(:, :, qMetric.thisChannelRaw(iUnit))) == PKS);
        end

        qMetric.numPeaksTroughs(iUnit) = numel(PKS) + numel(TRS);
        peakLoc = LOCS(PKS == max(PKS));
        if numel(peakLoc) > 1
            peakLoc = peakLoc(1);
        end
        troughLoc = LOCST(TRS == max(TRS));
        if numel(troughLoc) > 1
            troughLoc = troughLoc(1);
        end
        ephysParams.templateDuration(iUnit) = (peakLoc - troughLoc) / ephysData.ephys_sample_rate * 1e6;

        ephysParams.rawDuration(iUnit) = (peakLoc - troughLoc) / ephysData.ephys_sample_rate * 1e6;
        clearvars PKS TRS
    else
        qMetric.numPeaksTroughs(iUnit) = NaN;
    end

    %% spatial decay of raw mean waveform
    if ~isempty(curr_spikes_idx) && ~isempty(curr_spike_waveforms)
        maxXC = ephysData.channel_positions(thisChannelRaw, 1);
        maxYC = ephysData.channel_positions(thisChannelRaw, 2);
        chanDistances = ((ephysData.channel_positions(:, 1) - maxXC).^2 ...
            +(ephysData.channel_positions(:, 2) - maxYC).^2).^0.5;
        chansToPlot = find(chanDistances < param.chanDisMax);


        distancesOrder = chanDistances(chanDistances < param.chanDisMax);
        [~, sortIdx] = sort(distancesOrder);
        thisO = chansToPlot(sortIdx);

        maxChan = nan(size(chansToPlot, 1), 1);
        for iChanToPlot = 1:size(chansToPlot, 1)
            maxChan(iChanToPlot) = min(squeeze(qMetric.waveforms_mean(iUnit, :, thisO(iChanToPlot))));
        end

        if size(maxChan, 2) > 3
            qMetric.spatialDecay(iUnit, :) = [maxChan(1) / mean(maxChan(2:3)), mean(maxChan(2:3)) / mean(maxChan(4:end))];
        else
            qMetric.spatialDecay(iUnit, :) = [maxChan(1) / mean(maxChan(2:end)), NaN];
        end
    else
        qMetric.spatialDecay(iUnit, :) = [NaN, NaN];
    end


end
if param.drift == 0

    %%  template ampli
    %     thisChannelTempl = find(squeeze(max(ephysData.templates(thisUnit, :))) == ...
    %         max(squeeze(max(ephysData.templates(thisUnit, :, :)))));
    %     qMetric.thisChannelTempl(iUnit) = thisChannelTempl;
    waveformsTemp_mean = ephysData.template_waveforms(thisUnit, :);
    qMetric.waveformsTemp_mean(iUnit, :) = waveformsTemp_mean;

    wavefTPeak = max(waveformsTemp_mean);
    if numel(wavefTPeak) > 1
        wavefTPeak = mean(wavefTPeak);
    end
    wavefTTrough = abs(min(waveformsTemp_mean));
    if numel(wavefTTrough) > 1
        wavefTTrough = mean(wavefTTrough);
    end
    qMetric.waveformTemplAmpli(iUnit) = wavefTPeak + ...
        wavefTTrough;

    %% number of peaks of template

    minProminence = 0.2 * max(abs(squeeze(waveformsTemp_mean)));
    qMetric.waveformUnit(iUnit, :) = squeeze(waveformsTemp_mean);
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
    qMetric.numPeaksTroughsTemp(iUnit) = numel(PKS) + numel(TRS);
    %is somatic?
    peakLoc = LOCS(PKS == max(PKS));
    if numel(peakLoc) > 1
        peakLoc = peakLoc(1);

    end
    troughLoc = LOCST(TRS == max(TRS));
    if numel(troughLoc) > 1
        troughLoc = troughLoc(1);
    end

    if peakLoc > troughLoc
        ephysParams.somatic(iUnit) = 1;
    else
        ephysParams.somatic(iUnit) = 0;
    end
    if max(PKS) < max(TRS)
        ephysParams.trouP(iUnit) = 1;
    else
        ephysParams.trouP(iUnit) = 0;
    end
    if numel(PKS) == 1 && numel(TRS) == 2
        if LOCST(1) < LOCS(1) < LOCST(2)
            ephysParams.troughPeaktrough(iUnit) = 1;
        else
            ephysParams.troughPeaktrough(iUnit) = 1;
        end
    end
    ephysParams.templateDuration(iUnit) = (peakLoc - troughLoc) / ephysData.ephys_sample_rate * 1e6;

    %% spatial decay of template
    qMetric.thisChannelTempl(iUnit) = find(squeeze(max(ephysData.templates(thisUnit, :, :))) == ...
        max(squeeze(max(ephysData.templates(thisUnit, :, :)))));
    maxXC = ephysData.channel_positions(qMetric.thisChannelTempl(iUnit), 1);
    maxYC = ephysData.channel_positions(qMetric.thisChannelTempl(iUnit), 2);
    chanDistances = ((ephysData.channel_positions(:, 1) - maxXC).^2 ...
        +(ephysData.channel_positions(:, 2) - maxYC).^2).^0.5;
    chansToPlot = find(chanDistances < param.chanDisMax);


    distancesOrder = chanDistances(chanDistances < param.chanDisMax);
    [~, sortIdx] = sort(distancesOrder);
    thisO = chansToPlot(sortIdx);
    % figure();
    for iChanToPlot = 1:size(chansToPlot, 1)
        %    subplot(2, 15, iChanToPlot)
        %    plot(squeeze(ephysData.templates(thisUnit, :, thisO(iChanToPlot))))
        maxChan(iChanToPlot) = min(squeeze(ephysData.templates(thisUnit, :, thisO(iChanToPlot))));
    end

    %endP = size(chansToPlot, 1);
    midP = size(chansToPlot, 1) / 2;
    if size(chansToPlot, 1) > 3
        qMetric.spatialDecayTemp1(iUnit) = maxChan(1) / mean(maxChan(2:floor(midP)));
        qMetric.spatialDecayTemp2(iUnit) = mean(maxChan(2:floor(midP))) / mean(maxChan(floor(midP)+1:end));
    else
        qMetric.spatialDecayTemp1(iUnit) = maxChan(1) / mean(maxChan(2:end));
        qMetric.spatialDecayTemp2(iUnit) = NaN;
    end

    %% distance metrics: isolation distance, l-ratio, silhouette score %could add NN, d-prime
    if param.dist
        [isoD, Lratio, silhouetteScore, ~, ~, ~] = getDistanceMetricsJF(ephysData.pc_features, pc_feature_ind, thisUnit, ...
            numel(theseSpikes), theseSpikesIdx, allSpikesIdx, param.nChannelsIsoDist); % QQ degree of freedom parameter
        %QQ add need trough in at least on CCG.
        qMetric.isoD(iUnit) = isoD;
        qMetric.Lratio(iUnit) = Lratio;
        qMetric.silhouetteScore(iUnit) = silhouetteScore;
    end
end

%% contamination (rpvs)
[qMetric.fractionRPVchunk(iUnit), qMetric.numRPVchunk(iUnit)] = fractionRPviolationsJF( ...
    numel(theseSpikes), theseSpikes, param.tauR, param.tauC, timeChunks(end)-timeChunks(1)); %method from Hill et al., 2011

%% percentage spikes missing
% using gaussian %%QQ add test for this

try
    [percent_missing_ndtrAll, ~] = ampli_fit_prc_missJF(theseAmplis, 0);
catch
    percent_missing_ndtrAll = NaN;
end
qMetric.pMissing(iUnit) = percent_missing_ndtrAll;


qMetric.symSpikesMissing(iUnit) = prctgMissingSymetry(theseAmplis);

%% EPHYS PARAMS

%% acg
[ccg, t] = CCGBz([double(theseSpikes); double(theseSpikes)], [ones(size(theseSpikes, 1), 1); ...
    ones(size(theseSpikes, 1), 1) * 2], 'binSize', param.ACGbinSize, 'duration', param.ACGduration, 'norm', 'rate'); %function
%from the Zugaro lab mod. by Buzsaki lab-way faster than my own!
thisACG = ccg(:, 1, 1);
ephysParams.ACG(iUnit, :) = thisACG;
ephysParams.ACGtime(iUnit, :) = t;

% figure();bar(t*1000,ccg(:,1,1));
% xlabel('time(ms)');ylabel('sp/s');makepretty;

%% f.r.
ephysParams.spike_rateSimple(iUnit) = numel(theseSpikes) / (max(theseSpikes) - min(theseSpikes));

spiking_stat_window = max(theseSpikes) - min(theseSpikes);
spiking_stat_bins = [min(theseSpikes), max(theseSpikes)];

% Get firing rate across the session
bin_spikes = ...
    histcounts(theseSpikes, ...
    spiking_stat_bins);

min_spikes = 10;
use_spiking_stat_bins = bsxfun(@ge, bin_spikes, prctile(bin_spikes, 80, 2)) & bin_spikes > min_spikes;
spike_rate = sum(bin_spikes.*use_spiking_stat_bins, 2) ./ ...
    (sum(use_spiking_stat_bins, 2) * spiking_stat_window);
ephysParams.spike_rateAP(iUnit) = spike_rate;
ephysParams.acgFR(iUnit) = nanmean(ephysParams.ACG(iUnit, 190:200));

%% post spike suppression
acgfr = find(ephysParams.ACG(iUnit, 500:1000) >= ...
    nanmean(ephysParams.ACG(iUnit, 900:1000)));
%QQ smooth

if ~isempty(acgfr)
    acgfr = acgfr(1);
else
    acgfr = NaN;
end
ephysParams.postSpikeSuppression(iUnit) = acgfr;
% no longer used below (buggy function) 
% postS = find(ephysParams.ACG(iUnit, 500:1000) >= ...
%     spike_rate);
% [thisACG, ~, ~] = crosscorrelogram(theseSpikes, theseSpikes, [-param.ACGduration, param.ACGduration]); % 50 ms time window on either side, 0.1ms bins
% thisACG(thisACG == 0) = [];
% ephysParams.ACGBf(iUnit, :) = histcounts(thisACG, param.histBins); %  1 ms time bins
% 
% acgf = find(ephysParams.ACGBf(iUnit, 500:1000) >= ...
%     nanmean(ephysParams.ACGBf(iUnit, 600:900)));
% 
% if ~isempty(acgf)
%     acgf = acgf(1);
% else
%     acgf = NaN;
% end
% ephysParams.postSpikeSuppressionBf(iUnit) = acgf;

%% prop isi

long_isi_total = 0;
%isi_ratios = [];
for curr_bin = find(use_spiking_stat_bins)
    curr_spike_times = theseSpikes( ...
        theseSpikes > spiking_stat_bins(curr_bin) & ...
        theseSpikes < spiking_stat_bins(curr_bin+1));
    curr_isi = diff(curr_spike_times);

    long_isi_total = long_isi_total + sum(curr_isi(curr_isi > 2));

    %isi_ratios = [isi_ratios; (2 * abs(curr_isi(2:end)-curr_isi(1:end-1))) ./ ...
    %    (curr_isi(2:end) + curr_isi(1:end-1))]; %WRONG, see Holt 1996
end

ephysParams.prop_long_isi(iUnit) = long_isi_total / ...
    (sum(use_spiking_stat_bins(:)) * spiking_stat_window);
%cv2 = nanmean(isi_ratios);

%% cv
ephysParams.cv(iUnit) = nanstd(diff(theseSpikes)) / nanmean(diff(theseSpikes));

%% cv2
theseISIs = diff(theseSpikes);
cv2 = [];
for iISI = 1:numel(theseISIs) - 1
    cv2 = [cv2, (sqrt(2) * abs(theseISIs(iISI)-theseISIs(iISI+1))) ./ ...
        (nanmean(theseISIs(iISI)+theseISIs(iISI+1)))];
end
ephysParams.cv2(iUnit) = nanmean(cv2);
%Standard deviation of all adjacent ISIs, divide by their mean (*?2 so Poisson process has mean of 1)

%% Fano factor

%% skewISI
ephysParams.ISIskew(iUnit) = skewness(theseISIs);

%% max firing rate
[n, ~] = histcounts(theseSpikes, min(theseSpikes):param.maxFRbin:max(theseSpikes)); % firing rate over each xseconds
ephysParams.maxFR(iUnit) = max(n);

%% bursting things
[~, ~, numberBurstSpikes] = detectBurstsJF(theseSpikes);
ephysParams.numberBursts(iUnit) = numel(numberBurstSpikes);
ephysParams.avBurstSpikes(iUnit) = nanmean(numberBurstSpikes);
ephysParams.maxBurstSpikes(iUnit) = max(numberBurstSpikes);
end
