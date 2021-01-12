
%% Natural images decoding

%% 1. load all data (DMS, PS), check alignement (plot, + look at alignement to movement) -> dont manually
toLoad = readtable('C:/Users/Julie/Dropbox/NaturalImagesDatasetsForDecodingGrant.ods'); % get which animals and probes to load from each animal, and depth corr if necessary

theseAnimals = unique(toLoad.Animal);
ephysData = struct;
thisCount = 1;
for iAnimal = 1:size(theseAnimals)

    animal = theseAnimals{iAnimal};
    theseProbes = toLoad.Probe(find(strcmp(toLoad.Animal, animal))); %get probes to load
    theseDepthsCorr = toLoad.Depth(find(strcmp(toLoad.Animal, animal))); %get probes to load
    theseLocationsCorr = toLoad.Location(find(strcmp(toLoad.Animal, animal))); %get probes to load

    protString = 'mage';
    histoFile = AP_cortexlab_filenameJF(animal, [], [], 'histo', [], []);
    load(histoFile)
    probe2ephysFile = AP_cortexlab_filenameJF(animal, [], [], 'probe2ephys', [], []);
    load(probe2ephysFile)
    recInfoFile = AP_cortexlab_filenameJF(animal, [], [], 'acuteRecInfo', [], []);
    recInfo = readtable(recInfoFile);
    allenAt = loadStructureTreeJF(['C:\Users\Julie\Dropbox\Atlas\allenCCF\structure_tree_safe_2017.csv']);

    locations = {'CP'};
    %load probe
    uniqueD = unique(recInfo.Date);
    if isa(uniqueD, 'datetime')
        [dn, idx] = sort(uniqueD, 1, 'ascend');
    else
        uniqueD = uniqueD(contains(uniqueD, '/') | contains(uniqueD, '-')); %keep only dates (ie have a /) QQ
        try
            [dn, idx] = sort(datenum(uniqueD, 'dd/mm/yyyy'), 1, 'ascend'); %sort by ascending order
        catch
            [dn, idx] = sort(datenum(uniqueD, 'dd-mm-yyyy'), 1, 'ascend'); %sort by ascending order
        end
    end
    uniqueD = uniqueD(idx);

    for iP = 1:length(theseProbes)
        iProbe = theseProbes(iP);
        curr_day = probe2ephys(iProbe).day;
        curr_site = probe2ephys(iProbe).site;
        if isa(uniqueD, 'cell')
            cday = uniqueD{curr_day};
            dayBreaks = find(cday == '/' | cday == '-');
            if dayBreaks(1) == 2 && dayBreaks(2) == 4
                cday = strcat(cday(dayBreaks(2)+1:end), '-0', cday(dayBreaks(1)+1:dayBreaks(2)-1), ...
                    '-0', cday(1:dayBreaks(1)-1)); %in standard format
            elseif dayBreaks(1) == 2
                cday = strcat(cday(dayBreaks(2)+1:end), '-', cday(dayBreaks(1)+1:dayBreaks(2)-1), ...
                    '-0', cday(1:dayBreaks(1)-1)); %in standard format
            elseif dayBreaks(2) == 5
                cday = strcat(cday(dayBreaks(2)+1:end), '-0', cday(dayBreaks(1)+1:dayBreaks(2)-1), ...
                    '-', cday(1:dayBreaks(1)-1)); %in standard format
            else
                cday = strcat(cday(dayBreaks(2)+1:end), '-', cday(dayBreaks(1)+1:dayBreaks(2)-1), ...
                    '-', cday(1:dayBreaks(1)-1)); %in standard format
            end
            theseExpDates = strcmp(recInfo.Date, uniqueD{curr_day});
        else
            cday = uniqueD(curr_day);
            day = cday;
            theseExpDates = isbetween(recInfo.Date, uniqueD(curr_day), uniqueD(curr_day));
        end


        if ~iscell(recInfo.Site)
            theseExpSites = recInfo.Site == curr_site;
        else
            a = recInfo.Site;
            b = cellfun(@numel, recInfo.Site);
            a(b > 1 | b == 0) = {'0'};
            c = str2num(cell2mat(a));
            theseExpSites = c == curr_site;
        end
        if exist('recInfo.error', 'var')
            theseExpErrors = contains(recInfo.error, 'es') | ~isnan(recInfo.ephys_not_started) ...
                | ~isnan(recInfo.timeline_not_started);
        elseif isa(recInfo.ephys_not_started, 'string') && isa(recInfo.timeline_not_started, 'string')
            theseExpErrors = strlength(recInfo.ephys_not_started) > 0 ...
                | strlength(recInfo.timeline_not_started) > 0;
        elseif isa(recInfo.ephys_not_started, 'double') && isa(recInfo.timeline_not_started, 'double')
            theseExpErrors = ~isnan(recInfo.ephys_not_started) ...
                | ~isnan(recInfo.timeline_not_started);
        elseif isa(recInfo.ephys_not_started, 'cell') && isa(recInfo.timeline_not_started, 'double')
            theseExpErrors = strlength(recInfo.ephys_not_started) > 0 ...
                | ~isnan(recInfo.timeline_not_started);
        elseif isa(recInfo.ephys_not_started, 'cell') && isa(recInfo.timeline_not_started, 'cell')
            theseExpErrors = strlength(recInfo.ephys_not_started) > 0 ...
                | strlength(recInfo.timeline_not_started) > 0;
        elseif isa(recInfo.ephys_not_started, 'double') && isa(recInfo.timeline_not_started, 'string')
            theseExpErrors = ~isnan(recInfo.ephys_not_started) ...
                | strlength(recInfo.timeline_not_started) > 0;
        end
        theseExp = find(theseExpDates & theseExpSites & ~theseExpErrors);
        theseProtocols = recInfo.Protocol(theseExp);

        experiment = theseExp(contains(theseProtocols, protString));
        protocol = 'atural';
        experiments = AP_find_experimentsJF(animal, protocol, true); %experiments = experiments([experiments.imaging] & [experiments.ephys]); % (use only experiments with both widefield + ephys)
        %experiments = experiments([experiments.imaging] & [experiments.ephys]);
        experiments = experiments([experiments.ephys]);

        site = recInfo.Site(experiment);

        day = experiments(curr_day).day; % date
        experiment = recInfo.Protocol_number(experiment);
        verbose = false; % display load progress and some info figures
        load_parts.cam = false;
        load_parts.imaging = false;
        load_parts.ephys = true;

        lfp_channel = 'all';
        close all;
        AP_load_experimentJF_WIP;
        if dontAnalyze == 1 %try another method
            AP_load_experimentJF;
        end
        % plot alignement stuffs
        if isempty(theseDepthsCorr{iP, :})
            this_ccf = probe_ccf(iProbe);
            theseLocations = allenAt.acronym(this_ccf.trajectory_areas);
            theseLocationsInterest = contains(theseLocations, locations{1});
            theseDepths = this_ccf.probe_depths(theseLocationsInterest);
        else
            theseDepths_op = strsplit(theseDepthsCorr{iP, :}, '-');
            clearvars theseDepths
            theseDepths(1) = str2num(theseDepths_op{1, 1});
            theseDepths(2) = str2num(theseDepths_op{1, 2});
        end
        if ~isempty(theseDepths)
            %save data in structure
            theseTemplates = find(template_depths >= min(theseDepths) & template_depths <= max(theseDepths)); %correct depth units
            %find closest depth and get closest location
            if isempty(theseDepthsCorr{iP, :})
                A = repmat(this_ccf.probe_depths, [1, length(template_depths)]);
                [minValue, closestIndex] = min(abs(A-template_depths'));
                closestLocation = this_ccf.trajectory_coords(closestIndex, :);
            else
                closestLocation = NaN;
            end
            %raster general and for each trial condition
            theseSpikes = ismember(spike_templates, theseTemplates);
            ephysData(thisCount).spike_times_timeline = spike_times_timeline(theseSpikes);
            ephysData(thisCount).spike_templates = spike_templates(theseSpikes);
            ephysData(thisCount).spike_amplitudes = template_amplitudes(theseSpikes);
            ephysData(thisCount).templatesOI = theseTemplates;
            ephysData(thisCount).template_depths = template_depths(theseTemplates);
            ephysData(thisCount).template_wf = waveforms(theseTemplates, :);
            ephysData(thisCount).template_wfs = templates(theseTemplates, :, :);
            ephysData(thisCount).template_location = closestLocation;
            ephysData(thisCount).protocol = expDef;
            ephysData(thisCount).animal = animal;
            ephysData(thisCount).date = day;
            ephysData(thisCount).site = iProbe;
            ephysData(thisCount).experiment = experiment;
            ephysData(thisCount).location = theseLocationsCorr(iP);
            ephysData(thisCount).stimOn_times = stimOn_times;
            ephysData(thisCount).stimIDs = stimIDs;
            ephysData(thisCount).trial_conditions = trial_conditions;
            thisCount = thisCount + 1;
        end
        keep toLoad  theseAnimals  ephysData  thisCount iP iAnimal animal theseProbes theseDepthsCorr theseLocationsCorr protString ...
            probe_ccf probe2ephys recInfo allenAt locations uniqueD

    end
end
%save data on dropbox

%% sanity check: average responses
colorsT = lines(11);
figure();
thisWindow = [-0.2, 0.5];
psthBinSize = 0.01;
for i = 1:thisCount - 1
    theseSpikes = ephysData(i).spike_times_timeline;
    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, ephysData(i).stimOn_times, thisWindow, psthBinSize);
    [psthbaseline, bins, rasterX, rasterY, spikeCounts, binnedArrayBase] = psthAndBA(theseSpikes, ephysData(i).stimOn_times, [-0.2, 0], psthBinSize);
    psth = (psth - nanmean(psthbaseline)) / nanstd(psthbaseline);
    w = gausswin(10);
    y = filter(w, 1, (binnedArray ./ 0.01 - -nanmean(psthbaseline))/nanstd(psthbaseline));
    AP_errorfillJF([-0.2:0.01:0.5 - 0.01]', psth', std(y)./sqrt(size(binnedArray, 1)), colorsT(i, :), 0.1);
    makepretty;
    hold on;
end

ylabel('dFR/FR')
xlabel('time from nat. img (s)')
set(gcf, 'color', 'w');

%% 2. for each cell FR 1/2 trials against 1/2 -> r2, plot distribution
param.tauR = 0.0010; %refractory period time (s)
param.tauC = 0.0002; %censored period time (s)
param.maxNumPeak = 4;
param.maxRPV = 5;
param.maxPercMissing = 30;
unitCount = 0;
for i = 1:thisCount - 1
    theseUnits = unique(ephysData(i).spike_templates);
    ephysData(i).goodUnit = zeros(length(theseUnits), 1);
    for ii = 1:length(theseUnits)
        unitCount = unitCount + 1;
        theseSpikes = ephysData(i).spike_times_timeline(ephysData(i).spike_templates == theseUnits(ii));
        theseAmplis = ephysData(i).spike_amplitudes(ephysData(i).spike_templates == theseUnits(ii));
        thisWaveform = ephysData(i).template_wf(ii, :);
        theseWvs = ephysData(i).template_wfs(ii, :, :);
        [~, max_site] = max(max(abs(ephysData(i).template_wfs(ii, :, :)), [], 2), [], 3);
        if max_site < 5
            nearest_sites = [max_site:max_site + 7];
        elseif max_site > size(ephysData(i).template_wfs(ii, :, :), 3) - 7
            nearest_sites = [max_site - 7:max_site];
        else
            nearest_sites = [max_site - 4:max_site + 3];
        end

        %% quality metrics
        %check rfp
        [fractionRPVchunk, numRPVchunk] = fractionRPviolationsJF( ...
            numel(theseSpikes), theseSpikes, param.tauR, param.tauC, theseSpikes(end)-theseSpikes(1)); %method from Hill et al., 2011
        [ccg, t] = CCGBz([double(theseSpikes); double(theseSpikes)], [ones(size(theseSpikes, 1), 1); ...
            ones(size(theseSpikes, 1), 1) * 2], 'binSize', 0.01, 'duration', 0.5, 'norm', 'rate'); %function
        %from the Zugaro lab mod. by Buzsaki lab-way faster than my own!
        thisACG = ccg(:, 1, 1);
        %         figure(1);
        %         clf;
        %         plot(thisACG)
        %         title(num2str(fractionRPVchunk))
        %check amplitudes
        try
            [percent_missing_ndtrAll, ~] = ampli_fit_prc_missJF(theseAmplis, 0);
        catch
            percent_missing_ndtrAll = NaN;
        end


        %check number of peaks
        minProminence = 0.2 * max(abs(squeeze(thisWaveform)));

        %figure();plot(qMetric.waveform(iUnit, :))
        [PKS, LOCS] = findpeaks(squeeze(thisWaveform), 'MinPeakProminence', minProminence);
        [TRS, LOCST] = findpeaks(squeeze(thisWaveform)*-1, 'MinPeakProminence', minProminence);
        if isempty(TRS)
            TRS = min(squeeze(thisWaveform));
            if numel(TRS) > 1
                TRS = TRS(1);
            end
            LOCST = find(squeeze(thisWaveform) == TRS);
        end
        if isempty(PKS)
            PKS = max(squeeze(thisWaveform));
            if numel(PKS) > 1
                PKS = PKS(1);
            end
            LOCS = find(squeeze(thisWaveform) == PKS);
        end
        numPeaksTroughsTemp = numel(PKS) + numel(TRS);


        %is somatic?
        peakLoc = LOCS;
        if numel(peakLoc) > 1
            peakLoc = peakLoc(end);

        end
        troughLoc = LOCST(TRS == max(TRS));
        if numel(troughLoc) > 1
            troughLoc = troughLoc(1);
        end


        %check correlation/spatial decay
        [rho, pval] = corr(squeeze(theseWvs(:, :, nearest_sites)), squeeze(theseWvs(:, :, nearest_sites)));
        troughVals = min(squeeze(theseWvs(:, :, nearest_sites)));
        %         figure(2);
        %         clf;
        %         subplot(2, 5, 1)
        %         plot(thisWaveform)
        %         for iSite = 1:length(nearest_sites)
        %             subplot(2, 5, 2+iSite)
        %             plot(squeeze(theseWvs(:, :, nearest_sites(iSite))))
        %
        %             subplot(2, 5, 2)
        %             plot(squeeze(theseWvs(:, :, nearest_sites(iSite))))
        %             hold on;
        %         end
        %         title([num2str(min(troughVals)) num2str(max(troughVals))])

        if min(troughVals) < max(troughVals) * 2 && numPeaksTroughsTemp < param.maxNumPeak && peakLoc > troughLoc && ...
                fractionRPVchunk <= param.maxRPV && numel(theseSpikes) > 300 %&& percent_missing_ndtrAll < param.maxPercMissing
            %disp('goodUnit')
            ephysData(i).goodUnit(ii) = 1;
        else
            ephysData(i).goodUnit(ii) = 0;
        end

        %% FR 1/2 trials
        if ephysData(i).goodUnit(ii) == 1
            thisWindow = [0.05, 0.2];
            psthBinSize = 0.01;
            if max(ephysData(i).stimIDs) > 29
                BA=[];
                stim =[];
                PSTH=[];
                for iStim = 1:30
                    theseTrials = find(ephysData(i).stimIDs == iStim);

                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials(1:2:end)), thisWindow, psthBinSize); %psth aligned
                    FRunits(unitCount, iStim, 1) = nanmean(psth);
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials(2:2:end)), thisWindow, psthBinSize); %psth aligned
                    FRunits(unitCount, iStim, 2) = nanmean(psth);
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials), thisWindow, psthBinSize); %psth aligned
                    
                    BA = [BA, nanmean(binnedArray,2)];
                    stim = [stim, iStim*ones(size(binnedArray,1),1)];
                    PSTH = [PSTH; nanmean(binnedArray,1)];

                end
                FRunitsCorr(unitCount) = corr(squeeze(FRunits(unitCount, :, 1))', squeeze(FRunits(unitCount, :, 2))');
                kk=kstest(BA(1,:));
                load carsmall
                bb=BA';
                yv = bb(:);
                %p =  vartestn(yv,repmat(1:30,1,size(BA,1)),'TestType','LeveneAbsolute');                                                                                                                                         ,Model_Year,'TestType','LeveneAbsolute')
                FRunitsP(unitCount) = anova1(BA') ,[],'off');
                FRunitsPSTH(unitCount,:) = nanmean(PSTH,1);
                %FRunitsVAR(unitCount,:) = p;
                FRunitaWA(unitCount) = welchanova([yv, repmat(1:30,1,size(BA,1))'],0.05);
                
                %[p,tbl,stats] = anova1(BA)
                %         figure();
                %         scatter(squeeze(FRunits(unitCount,:,1)), squeeze(FRunits(unitCount,:,2)))
            else
                FRunits(unitCount, 1:10, 2) = NaN;
                FRunits(unitCount, 1:30, 1) = NaN;
                FRunitsCorr(unitCount) = NaN;
                FRunitsPSTH(unitCount,1:15) = NaN;
                
            end
        else
            FRunits(unitCount, 1:30, 2) = NaN;
            FRunits(unitCount, 1:30, 1) = NaN;
            FRunitsCorr(unitCount) = NaN;
            FRunitsPSTH(unitCount,1:15) = NaN;
        end


    end
    keep FRunits ii theseUnits unitCount i ephysData thisCount param FRunitsCorr FRunitsP BA
end
%summary
figure();
hist(FRunitsCorr, 20)
xlabel('r2 coeff.')
ylabel('# of units')
makepretty;
prctile(FRunitsCorr,50)
%example cells 
ff=find(FRunitsCorr<-0.3&FRunitsCorr>-0.5);%>0.5; <0.1 > -0.1 ;
for iW=ff
    figure();
    %subplot(221)
    scatter(squeeze(FRunits(iW,:,1)), squeeze(FRunits(iW,:,2)))
    xlabel('FR 1/2 trials')
    ylabel('FR 1/2 trials')
    makepretty;
    %subplot(222) 
    %subplot(223)
    
    %subplot(224)
end
FRunitsCorr(ff(21))
FRunitsCorr(ff(35))
FRunitsCorr(ff(28))
FRunitsCorr(ff(37))
FRunitsCorr(ff(35))
%% 3. one-way ANOVA: -way anova for each cell: is there a main effect "stim" on firing rate  ? pvalue for each cell. can then plot the districbution of pvalues and see which ones are signif.
figure();
del=round(FRunitsP*10000)==0;
aa=FRunitsP;
aa(del)=[];
hist(aa, 40)
ylabel('unit #')
xlabel('one-way ANOVA p-value')
makepretty;
%% 4. decoder
%%Training Side
% random three class data with target matrix -- [9X3] 9 observation with 3 features
data = FRunits(:,:,1)';
target = 1:30;
% create a naive bayes model
% data must not have zero variance
% var(data(target==1,:)) for checking variance for class 1
nb = fitcnb(data,target); % nb is the trained model. save it at end for doing testing
save('nb.mat','nb');
% train performance
label = predict(nb,data);
perf=sum(label==target)/size(label,1); % performance in the range of 0 to 1
%%Testing Side
% for testing load the trained model
load('nb.mat');
testdata = [10 0 0]; % take 1 new unknown observation and give to trained model
Group = predict(nb,testdata);

%% TO DO: EXAMPLE 