
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
thisCount=11;
       % FRunitsPSTHtrainStim = nan(1,30);
       %         FRunitsPSTHtestStim=nan(1,30);
                %FRunitsPSTHtrain=nan(1,15);
                %FRunitsPSTHtest=nan(1,15);
 unitIdx=[]; 
 siteIdx=[];
 goodCount = 1;
for i = 1:thisCount - 1
    theseUnits = unique(ephysData(i).spike_templates);
    ephysData(i).goodUnit = zeros(length(theseUnits), 1);
    siteIdx = [siteIdx; i*ones(length(theseUnits),1)];
    unitIdx = [unitIdx, [1:length(theseUnits)]];
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
                stimtrain =[];
                stimtest =[];
                PSTHtrain=[];
                PSTHtest=[];
                trainData1=[];
                trainData2=[];
                trainData3=[];
                trainData4=[];
                trainData5=[];
                trainData6=[];
                trainData7=[];
                trainData8=[];
                trainData9=[];
                trainData10=[];

                for iStim = 1:30
                    theseTrials = find(ephysData(i).stimIDs == iStim);

                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray1] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials(1:2:end)), thisWindow, psthBinSize); %psth aligned
                    FRunits(unitCount, iStim, 1) = nanmean(psth);
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray2] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials(2:2:end)), thisWindow, psthBinSize); %psth aligned
                    FRunits(unitCount, iStim, 2) = nanmean(psth);
                    FRunitsCorr(goodCount) = corr(squeeze(FRunits(unitCount, :, 1))', squeeze(FRunits(unitCount, :, 2))');
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials), thisWindow, psthBinSize); %psth aligned
                    
                    BA = [BA, nanmean(binnedArray,2)];
                    stimtrain = [stimtrain; iStim];
                    stimtest = [stimtest; iStim];
                    PSTHtrain = [PSTHtrain; nanmean(binnedArray1,1)];
                    PSTHtest = [PSTHtest; nanmean(binnedArray2,1)];
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray1] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials(1:10:end)), thisWindow, psthBinSize); %psth aligned
                    
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray2] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials(2:10:end)), thisWindow, psthBinSize); %psth aligned
                   
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray3] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials(3:10:end)), thisWindow, psthBinSize); %psth aligned
                   
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray4] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials(4:10:end)), thisWindow, psthBinSize); %psth aligned
                  
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray5] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials(5:10:end)), thisWindow, psthBinSize); %psth aligned
                    
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray6] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials(6:10:end)), thisWindow, psthBinSize); %psth aligned
                    
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray7] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials(7:10:end)), thisWindow, psthBinSize); %psth aligned
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray8] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials(8:10:end)), thisWindow, psthBinSize); %psth aligned
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray9] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials(9:10:end)), thisWindow, psthBinSize); %psth aligned
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray10] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials(10:10:end)), thisWindow, psthBinSize); %psth aligned
                    trainData1= [trainData1; nanmean(nanmean(binnedArray1))]; 
                    trainData2= [trainData2; nanmean(nanmean(binnedArray2))];
                    trainData3= [trainData3; nanmean(nanmean(binnedArray3))];
                    trainData4= [trainData4; nanmean(nanmean(binnedArray4))];
                    trainData5= [trainData5; nanmean(nanmean(binnedArray5))];
                    trainData6= [trainData6; nanmean(nanmean(binnedArray6))];
                   trainData7= [trainData7; nanmean(nanmean(binnedArray7))];
                    trainData8= [trainData8; nanmean(nanmean(binnedArray8))];
                    trainData9= [trainData9; nanmean(nanmean(binnedArray9))];
                    trainData10= [trainData10; nanmean(nanmean(binnedArray10))];
                    %PSTHtest = [PSTHtest; nanmean(binnedArray2,1)];
                    

                end
                 kk=kstest(BA(1,:));
                load carsmall
                bb=BA';
                yv = bb(:);
                %p =  vartestn(yv,repmat(1:30,1,size(BA,1)),'TestType','LeveneAbsolute');                                                                                                                                         ,Model_Year,'TestType','LeveneAbsolute')
                FRunitsP(unitCount) = anova1(BA ,[],'off');
                FRunitsPKS(unitCount) = kruskalwallis(BA, [], 'off');
                if FRunitsPKS(unitCount)<0.05
                    FRunitsCorr(goodCount) = corr(squeeze(FRunits(unitCount, :, 1))', squeeze(FRunits(unitCount, :, 2))');
               
                FRunitsPSTHtrain((goodCount-1)*30+1:(goodCount)*30,1) = trainData1;
                FRunitsPSTHtrain((goodCount-1)*30+1:(goodCount)*30,2) = trainData2;
                FRunitsPSTHtrain((goodCount-1)*30+1:(goodCount)*30,3) = trainData3;
                FRunitsPSTHtrain((goodCount-1)*30+1:(goodCount)*30,4) = trainData4;
                FRunitsPSTHtrain((goodCount-1)*30+1:(goodCount)*30,5) = trainData5;
                FRunitsPSTHtrain((goodCount-1)*30+1:(goodCount)*30,6) = trainData6;
                FRunitsPSTHtrain((goodCount-1)*30+1:(goodCount)*30,7) = trainData7;
                FRunitsPSTHtrain((goodCount-1)*30+1:(goodCount)*30,8) = trainData8;
                FRunitsPSTHtrain((goodCount-1)*30+1:(goodCount)*30,9) = trainData9;
                FRunitsPSTHtrain((goodCount-1)*30+1:(goodCount)*30,10) = trainData10;
                FRunitsPSTHtrainStim((goodCount-1)*30+1:(goodCount)*30) = stimtrain;
                goodCount = goodCount+1;
                end
                %FRunitsPSTHtest((unitCount-1)*30+1:(unitCount)*30,:) = PSTHtest;
                %FRunitsVAR(unitCount,:) = p;
                FRunitsWA(unitCount) = welchanova([yv, repmat(1:30,1,size(BA,1))'],0.05);
                
                FRunitsPSTHtestStim((unitCount-1)*30+1:(unitCount)*30) = stimtest;
                %[p,tbl,stats] = anova1(BA)
                %         figure();
                %         scatter(squeeze(FRunits(unitCount,:,1)), squeeze(FRunits(unitCount,:,2)))
            else
                FRunits(unitCount, 1:10, 2) = NaN;
                FRunits(unitCount, 1:30, 1) = NaN;
                FRunitsPKS(unitCount) = NaN;
                %FRunitsCorr(unitCount) = NaN;
%                 FRunitsPSTHtrain(unitCount,1:15) = NaN;
%                 FRunitsPSTHtest(unitCount,1:15) = NaN;
                FRunitsWA(unitCount)=NaN;
               
                
            end
        else
            FRunits(unitCount, 1:30, 2) = NaN;
            FRunits(unitCount, 1:30, 1) = NaN;
            FRunitsPKS(unitCount) = NaN;
%             FRunitsCorr(unitCount) = NaN;
%             FRunitsPSTHtrain(unitCount,1:15) = NaN;
            FRunitsWA(unitCount)=NaN;
        end


    end
    keep FRunitsPKS FRunits ii theseUnits unitCount i ephysData thisCount param FRunitsCorr FRunitsP BA FRunitsWA FRunitsPSTHtrain ...
        FRunitsPSTHtestStim FRunitsPSTHtrainStim siteIdx unitIdx goodCount
end
size(FRunitsCorr)%1031
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
hist(aa, 20)
ylabel('unit #')
xlabel('one-way ANOVA p-value')
makepretty;

figure();
del=round(FRunitsWA*10000)==0;
aa=FRunitsWA;
aa(del)=[];
hist(aa,20)
ylabel('unit #')
xlabel('welch''s ANOVA p-value')
makepretty;

figure();
del=round(FRunitsPKS*10000)==0;
aa=FRunitsPKS;
aa(del)=[];
hist(aa,30)
ylabel('unit #')
xlabel('Kruskal-Wallis p-value')
makepretty;
%% 4. decoder - MAP
%%Training Side
% random three class data with target matrix -- [9X3] 9 observation with 3 features
data = FRunitsPSTHtrain(:,:,1)';
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

%% 5. decoder - SVM one vs all 
%%Training Side
% AND gate data with target matrix -- [4X2] 4 observation with 2 features
data = [0 0;
        0 1;
        1 0;
        1 1];
target = [1;1;1;2]; % [0 0 0 1] is representd as [1 1 1 2] 

data = FRunitsPSTHtrain; 
target = FRunitsPSTHtrainStim; 
target(target==0)=NaN;
%train an binary SVM model
SVMstructc = fitcsvm(data,target,'Standardize',true); % SVMstructc is the trained model. save it at end for doing testing
save('SVMstructc.mat','SVMstructc');
% train performance
Group = predict(SVMstructc,data) % give the data to model for checking its training level
perf=sum(Group==target)/size(Group,1) % performance in the range of 0 to 1
%%Testing Side
% for testing load the trained model
load('SVMstructc.mat');
testdata = [0 0]; % take 1 new unknown observation and give to trained model
Group = predict(SVMstructc,testdata)


%% 6. decoder - KNN
GreatUnits = find(FRunitsPKS<0.05  );
for iStim = 1:30
    data1(iStim, :)= FRunitsPSTHtrain(find( FRunitsPSTHtrainStim(1:size(FRunitsPSTHtrain,1))==iStim),1);
    data2(iStim, :)= FRunitsPSTHtrain(find( FRunitsPSTHtrainStim(1:size(FRunitsPSTHtrain,1))==iStim),2);
    data3(iStim, :)= FRunitsPSTHtrain(find( FRunitsPSTHtrainStim(1:size(FRunitsPSTHtrain,1))==iStim),3);
    data4(iStim, :)= FRunitsPSTHtrain(find( FRunitsPSTHtrainStim(1:size(FRunitsPSTHtrain,1))==iStim),4);
    data5(iStim, :)= FRunitsPSTHtrain(find( FRunitsPSTHtrainStim(1:size(FRunitsPSTHtrain,1))==iStim),5);
    data6(iStim, :)= FRunitsPSTHtrain(find( FRunitsPSTHtrainStim(1:size(FRunitsPSTHtrain,1))==iStim),6);
    data7(iStim, :)= FRunitsPSTHtrain(find( FRunitsPSTHtrainStim(1:size(FRunitsPSTHtrain,1))==iStim),7);
    data8(iStim, :)= FRunitsPSTHtrain(find( FRunitsPSTHtrainStim(1:size(FRunitsPSTHtrain,1))==iStim),8);
    data9(iStim, :)= FRunitsPSTHtrain(find( FRunitsPSTHtrainStim(1:size(FRunitsPSTHtrain,1))==iStim),9);
    data10(iStim, :)= FRunitsPSTHtrain(find( FRunitsPSTHtrainStim(1:size(FRunitsPSTHtrain,1))==iStim),10);
    
end
kNNeigh = 7; %try 5, 7
allData = [data1; data2; data3; data4; data5;data6; data7; data8; data9;data10];
for iCV = 1:5
    nums = [(iCV - 1)*30*2+1, iCV * 2*30]; 
data =  allData([1:nums(1)-1, nums(2)+1:end],:); 
target = repmat(1:30,[1,8]); 
kNNModel = fitcknn(data,target,'NumNeighbors',kNNeigh); % kNNModel is the trained model. save it at end for doing testing
save('kNNModel.mat','kNNModel');
% train performance
label = predict(kNNModel,data);
perf=sum(label==target)/size(label,1); % performance in the range of 0 to 1
%%Testing Side
% for testing load the trained model
load('kNNModel.mat');
testdata = allData(nums(1):nums(2),:); % take 1 new unknown observation and give to trained model
Group = predict(kNNModel,testdata);

pCorr(iCV) = numel(find(Group'-repmat([1:30],[1,2]) ==0))/60*100;

end
chanceLevel = (1/30)*100;

figure();
%boxplot(pCorr)
hold on;
scatter(ones(size(pCorr)).*(1+(rand(size(pCorr))-0.5)/2),pCorr,'b','filled')
hold on;
line([0 2], [chanceLevel, chanceLevel], 'Color','r');
ylabel('CV-split classification accuracy')
xlabel('KNN classifier');
makepretty;
%how many times isd the label right? 
%% TO DO: EXAMPLES
GreatUnits = find(FRunitsPKS<0.05  );

colorsO = lines(30); 
for iGreatUnit = 1:length(GreatUnits)
    iGreatUnit = iGreatUnit + 1;
    i=siteIdx(GreatUnits(iGreatUnit));
    ii=unitIdx(GreatUnits(iGreatUnit));
    theseUnits = unique(ephysData(i).spike_templates);
    theseSpikes = ephysData(i).spike_times_timeline(ephysData(i).spike_templates == theseUnits(ii));
    thisWindow = [-0.2, 0.5];
    psthBinSize = 0.01;
    binnFull = [];
    stim=[];
    for iStim = 1:30
        theseTrials = find(ephysData(i).stimIDs == iStim);
        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials), thisWindow, psthBinSize); %psth aligned
        binnFull = [binnFull; binnedArray];
        meanP(iStim) = nanmean(psth(11:31)); 
        stim =[stim; iStim*ones(size(binnedArray,1),1)];
    end
    %[ss, sidx]=sort(meanP);
    sidx=1:30;
    clf; 
    binnFull(binnFull > 1) = 1;
    sortedBA=[];
    for iStim = sidx
        sortedBA = [sortedBA; binnFull(stim==iStim,:)];
    end
    %h = imagesc(-0.1:0.01:0.3-0.01, [], 1-binnFull);
    [yPoints,xPoints] = find(sortedBA==1);
    %xPoints = 0.0xPoints  -0.1; 
    tt=repmat(-0.2:0.01:0.5-0.1, [1,600]);
    plot(tt(xPoints),yPoints,'.k');
    xlim([-0.2, 0.5])
    colormap(gray)
    yy = ylim;
    for iOr = 1:30
        line([0, 0], [1 + round(yy(2)/30) * (iOr - 1), round(yy(2)/30) * (iOr)], 'Color', colorsO(iOr, :))
        makepretty;
    end
    xlabel('time from stim onset (s)')
    ylabel('trial # (sorted by image)')
    makeprettyLite;
    %site and unit # 
    %get binned Array for each stim 
    
end
%% luminosity of each natural image 
load('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\naturalImages\img11')
sum(sum(img))
% -> pretty comparable 

%% cells are selective same images ?
%%get data aligned to stim 
param.tauR = 0.0010; %refractory period time (s)
param.tauC = 0.0002; %censored period time (s)
param.maxNumPeak = 4;
param.maxRPV = 5;
param.maxPercMissing = 30;
unitCount = 0;
thisCount=11;
       % FRunitsPSTHtrainStim = nan(1,30);
       %         FRunitsPSTHtestStim=nan(1,30);
                %FRunitsPSTHtrain=nan(1,15);
                %FRunitsPSTHtest=nan(1,15);
 unitIdx=[]; 
 siteIdx=[];
 goodCount = 1;
 clearvars FRunitsTime
 warning off
for i = 1:thisCount - 1
    theseUnits = unique(ephysData(i).spike_templates);
    ephysData(i).goodUnit = zeros(length(theseUnits), 1);
    siteIdx = [siteIdx; i*ones(length(theseUnits),1)];
    unitIdx = [unitIdx, [1:length(theseUnits)]];
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
                stimtrain =[];
                stimtest =[];
                PSTHtrain=[];
                PSTHtest=[];
                allTrialsdata=[];
                allTrialsdata1=[];
                allTrialsdata2=[];
                allTrialTimedata=[];
               

                for iStim = 1:30
                    theseTrials = find(ephysData(i).stimIDs == iStim);

                    
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials), thisWindow, psthBinSize); %psth aligned
                    
                    BA = [BA, nanmean(binnedArray,2)];
                    stimtrain = [stimtrain; iStim];
                    stimtest = [stimtest; iStim];
                    allTrialsdata= [allTrialsdata; nanmean(nanmean(binnedArray))];
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray1] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials(1:2:end)), thisWindow, psthBinSize); %psth aligned
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray2] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials(2:2:end)), thisWindow, psthBinSize); %psth aligned
                    allTrialsdata1= [allTrialsdata1; nanmean(nanmean(binnedArray1))];
                    allTrialsdata2= [allTrialsdata2; nanmean(nanmean(binnedArray2))];
                   [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials), [-0.2,0.5], psthBinSize); %psth aligned
                    
                    allTrialTimedata = [allTrialTimedata; nanmean(binnedArray)];
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray1] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials(1:2:end)), thisWindow, psthBinSize); %psth aligned
                    FRunits(unitCount, iStim, 1) = nanmean(psth);
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray2] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials(2:2:end)), thisWindow, psthBinSize); %psth aligned
                    FRunits(unitCount, iStim, 2) = nanmean(psth);
                   
                    %PSTHtest = [PSTHtest; nanmean(binnedArray2,1)];
                    

                end
                 FRunitsCorr(goodCount) = corr(squeeze(FRunits(unitCount, :, 1))', squeeze(FRunits(unitCount, :, 2))');
                FRunitsPSTHtrain(goodCount, :) = allTrialsdata;
                FRunitsTime(goodCount,:,:)=allTrialTimedata; 
                goodCount = goodCount+1;
                
                %FRunitsPSTHtest((unitCount-1)*30+1:(unitCount)*30,:) = PSTHtest;
                %FRunitsVAR(unitCount,:) = p;
               
                %[p,tbl,stats] = anova1(BA)
                %         figure();
                %         scatter(squeeze(FRunits(unitCount,:,1)), squeeze(FRunits(unitCount,:,2)))
            else
               
               
                
            end
        else
            
        end


    end
    keep  ii theseUnits unitCount i ephysData thisCount param  BA FRunitsWA FRunitsPSTHtrain ...
        siteIdx unitIdx goodCount FRunitsTime FRunitsPSTHtrain1 FRunitsPSTHtrain2
end
size( FRunitsCorr)
%%MUA time * stim matrix 
figure(); 
imagesc(-0.2:0.05:0.5-0.05,[],squeeze(nanmean(FRunitsTime)))
colormap(brewermap([],'*RdBu'))
xlabel('time from image onset (s)')
ylabel('image #')
makepretty; 

%%cell * stim matrix (average)
figure(); 
%get max 

%ops.nCall = [30,2];

[val,ii]=max(zscore(FRunitsPSTHtrain,[],2),[],2);
[sV, sI] =sort(ii);
zz=zscore(FRunitsPSTHtrain,[],2);
imagesc(zz(sI,:))
colormap(brewermap([],'*RdBu'))
ylabel('unit # (sorted by max response)')
xlabel('image #')
makepretty; 

%%cell* stim matrix, sorted with rtaster map

%get max 
ops=struct;
%ops.nCall = [30,2];
[isort1, isort2, Sm] = mapTmap(zscore(FRunitsPSTHtrain,[],2), ops);

zz=zscore(FRunitsPSTHtrain,[],2);
figure();
imagesc(zz(isort1,isort2))
colormap(brewermap([],'*RdBu'))
ylabel('unit # (sorted by max response)')
xlabel('image #')
makepretty; 

%%MSN vs TAN vs FSI 
param.tauR = 0.0010; %refractory period time (s)
param.tauC = 0.0002; %censored period time (s)
param.maxNumPeak = 4;
param.maxRPV = 5;
param.maxPercMissing = 30;
unitCount = 0;
thisCount=11;
       % FRunitsPSTHtrainStim = nan(1,30);
       %         FRunitsPSTHtestStim=nan(1,30);
                %FRunitsPSTHtrain=nan(1,15);
    acgA=[];
    wvA=[];%FRunitsPSTHtest=nan(1,15);
 unitIdx=[]; 
 siteIdx=[];
 cellT=[];
 goodCount = 1;
 clearvars FRunitsTime
for i = 1:thisCount - 1
    theseUnits = unique(ephysData(i).spike_templates);
    ephysData(i).goodUnit = zeros(length(theseUnits), 1);
    siteIdx = [siteIdx; i*ones(length(theseUnits),1)];
    unitIdx = [unitIdx, [1:length(theseUnits)]];
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
            ones(size(theseSpikes, 1), 1) * 2], 'binSize', 0.001, 'duration', 1, 'norm', 'rate'); %function
        %from the Zugaro lab mod. by Buzsaki lab-way faster than my own!
        thisACG = ccg(:, 1, 1);
        acgf = find(thisACG( 500:1000) >= ...
    nanmean(thisACG( 600:900)));
ephysData(i).acg(ii,:)=thisACG;
if ~isempty(acgf)
    acgf = acgf(1);
else
    acgf = NaN;
end
postSpikeSuppressionBf = acgf;
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
            if (-abs(troughLoc) + abs(peakLoc)) * 1e6 /30000 <= 400
                ephysData(i).celltype(ii) = 2;%fsi
            elseif postSpikeSuppressionBf < 50 && (-abs(troughLoc) + abs(peakLoc)) * 1e6 /30000 > 500
                ephysData(i).celltype(ii) = 1;%msn
            elseif postSpikeSuppressionBf>=50
                ephysData(i).celltype(ii) = 3;%tan
            end
        else
            ephysData(i).goodUnit(ii) = 0;
        end

        %% FR 1/2 trials
        if ephysData(i).goodUnit(ii) == 1
            thisWindow = [0.05, 0.2];
            psthBinSize = 0.01;
            if max(ephysData(i).stimIDs) > 29
                BA=[];
                stimtrain =[];
                stimtest =[];
                PSTHtrain=[];
                PSTHtest=[];
                allTrialsdata=[];
                allTrialsdata1=[];
                allTrialsdata2=[];
                allTrialTimedata=[];
               

                for iStim = 1:30
                    theseTrials = find(ephysData(i).stimIDs == iStim);

                    
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials), thisWindow, psthBinSize); %psth aligned
                    
                    BA = [BA, nanmean(binnedArray,2)];
                    stimtrain = [stimtrain; iStim];
                    stimtest = [stimtest; iStim];
                    allTrialsdata= [allTrialsdata; nanmean(nanmean(binnedArray))];
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray1] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials(1:2:end)), thisWindow, psthBinSize); %psth aligned
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray2] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials(2:2:end)), thisWindow, psthBinSize); %psth aligned
                    allTrialsdata1= [allTrialsdata1; nanmean(nanmean(binnedArray1))];
                    allTrialsdata2= [allTrialsdata2; nanmean(nanmean(binnedArray2))];
                   [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials), [-0.2,0.5], psthBinSize); %psth aligned
                    
                    allTrialTimedata = [allTrialTimedata; nanmean(binnedArray)];
                     [psth, bins, rasterX, rasterY, spikeCounts, binnedArray1] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials(1:2:end)), thisWindow, psthBinSize); %psth aligned
                    FRunits( iStim, 1) = nanmean(psth);
                    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray2] = psthAndBA(theseSpikes, ephysData(i).stimOn_times(theseTrials(2:2:end)), thisWindow, psthBinSize); %psth aligned
                    FRunits(iStim, 2) = nanmean(psth);
                   
                    %PSTHtest = [PSTHtest; nanmean(binnedArray2,1)];
                    

                end
                cellT = [cellT;  ephysData(i).celltype(ii)];
                acgA = [acgA, thisACG];
                wvA = [wvA; thisWaveform];
                FRunitsPSTHtrain(goodCount, :) = allTrialsdata;
                FRunitsPSTHtrain1(goodCount, :) = allTrialsdata1;
                FRunitsPSTHtrain2(goodCount, :) = allTrialsdata2;
                FRunitsTime(goodCount,:,:)=allTrialTimedata; 
                
                FRunitsCorr(goodCount) = corr(squeeze(FRunits( :, 1)), squeeze(FRunits( :, 2)));
               
                goodCount = goodCount+1;
                
                %FRunitsPSTHtest((unitCount-1)*30+1:(unitCount)*30,:) = PSTHtest;
                %FRunitsVAR(unitCount,:) = p;
               
                %[p,tbl,stats] = anova1(BA)
                %         figure();
                %         scatter(squeeze(FRunits(unitCount,:,1)), squeeze(FRunits(unitCount,:,2)))
            else
               
               
                
            end
        else
            
        end


    end
    keep  ii theseUnits unitCount i ephysData thisCount param  BA FRunitsWA FRunitsPSTHtrain ...
        siteIdx unitIdx goodCount FRunitsTime cellT acgA wvA FRunitsPSTHtrain1  FRunitsPSTHtrain2 FRunitsCorr
end
size(FRunitsCorr)
%%ACG and WV check
figure(); 
subplot(231)
plot(nanmean(acgA(:,cellT==1),2),'r')
xlim([0 1001])
ylabel('FR')
makepretty;
subplot(232)
plot(nanmean(acgA(:,cellT==2),2),'b')
xlim([0 1001])
xlabel('time (ms)')
makepretty;
subplot(233)
plot(nanmean(acgA(:,cellT==3),2),'g')
xlim([0 1001])
makepretty;
subplot(234)
t_wv=[1:size(wvA,2)] * 1e3/ 30000;
plot(t_wv, nanmean(wvA(cellT==1,:)),'r')
xlim([t_wv(1), t_wv(end)])
makepretty;
subplot(235)
plot(t_wv,nanmean(wvA(cellT==2,:)),'b')
xlim([t_wv(1), t_wv(end)])
xlabel('time (ms)')
makepretty;
subplot(236)
plot(t_wv,nanmean(wvA(cellT==3,:)),'g')
xlim([t_wv(1), t_wv(end)])
makepretty;
%%matrix mua 
figure(); 
subplot(131)
imagesc(-0.2:0.05:0.5-0.05,[],squeeze(nanmean(FRunitsTime(cellT==1,:,:))))
colormap(brewermap([],'*RdBu'))
title('MSN')
%xlabel('time from image onset (s)')
ylabel('image #') 
makepretty;
subplot(132)
imagesc(-0.2:0.05:0.5-0.05,[],squeeze(nanmean(FRunitsTime(cellT==2,:,:))))
colormap(brewermap([],'*RdBu'))
title('FSI')
xlabel('time from image onset (s)')
%ylabel('image #')
makepretty; 
subplot(133)
imagesc(-0.2:0.05:0.5-0.05,[],squeeze(nanmean(FRunitsTime(cellT==3,:,:))))
colormap(brewermap([],'*RdBu'))
title('TAN')
%xlabel('time from image onset (s)')
%ylabel('image #')
makepretty; 
%%
figure();
subplot(131)
[val,ii]=max(zscore(FRunitsPSTHtrain(cellT==1,:),[],2),[],2);
[sV, sI] =sort(ii);
zz=zscore(FRunitsPSTHtrain(cellT==1,:),[],2);
imagesc(zz(sI,:))
colormap(brewermap([],'*RdBu'))
ylabel('unit # (sorted by max response)')
xlabel('image #')
title('MSN')
makepretty; 

subplot(132)
[val,ii]=max(zscore(FRunitsPSTHtrain(cellT==2,:),[],2),[],2);
[sV, sI] =sort(ii);
zz=zscore(FRunitsPSTHtrain(cellT==2,:),[],2);
imagesc(zz(sI,:))
colormap(brewermap([],'*RdBu'))
ylabel('unit # (sorted by max response)')
xlabel('image #')
title('FSI')
makepretty; 

subplot(133)
[val,ii]=max(zscore(FRunitsPSTHtrain(cellT==3,:),[],2),[],2);
[sV, sI] =sort(ii);
zz=zscore(FRunitsPSTHtrain(cellT==3,:),[],2);
imagesc(zz(sI,:))
colormap(brewermap([],'*RdBu'))
ylabel('unit # (sorted by max response)')
xlabel('image #')
title('TAN')
makepretty; 
%% MAX, 1/2 trials 
frclean = FRunitsCorr(~isnan(FRunitsCorr))
theseCells = FRunitsCorr >= 0.5; 
% figure();
% H=dendrogram(Z);
figure();
subplot(131)
ops=struct;
%ops.nCall = [30,2];
[isort1, isort2, Sm] = mapTmap(zscore(FRunitsPSTHtrain(theseCells' & cellT==1,:),[],2), ops);
[val,ii]=max(zscore(FRunitsPSTHtrain1(cellT==1,:),[],2),[],2);
[sV, sI] =sort(ii);
tree= linkage(FRunitsPSTHtrain(cellT==1,:));
D = pdist(FRunitsPSTHtrain(cellT==1,:));
leafOrder = optimalleaforder(tree,D)
zz=zscore(FRunitsPSTHtrain(theseCells' &cellT==1,:),[],2);
imagesc(zz(isort1,isort2))
colormap(brewermap([],'*RdBu'))
ylabel('unit # (sorted by max response)')
xlabel('image #')
title('MSN')
makepretty; 

subplot(132)
[isort1, isort2, Sm] = mapTmap(zscore(FRunitsPSTHtrain(theseCells' &cellT==2,:),[],2), ops);
[val,ii]=max(zscore(FRunitsPSTHtrain1(cellT==2,:),[],2),[],2);
[sV, sI] =sort(ii);
tree= linkage(FRunitsPSTHtrain(cellT==2,:));
D = pdist(FRunitsPSTHtrain(cellT==2,:));
leafOrder = optimalleaforder(tree,D)
zz=zscore(FRunitsPSTHtrain(theseCells' &cellT==2,:),[],2);
imagesc(zz(isort1,isort2))
colormap(brewermap([],'*RdBu'))
ylabel('unit # (sorted by max response)')
xlabel('image #')
title('FSI')
makepretty; 

subplot(133)
[isort1, isort2, Sm] = mapTmap(zscore(FRunitsPSTHtrain(theseCells' &cellT==3,:),[],2), ops);
[val,ii]=max(zscore(FRunitsPSTHtrain1(cellT==3,:),[],2),[],2);
[sV, sI] =sort(ii);
tree= linkage(FRunitsPSTHtrain(cellT==3,:));
D = pdist(FRunitsPSTHtrain(cellT==3,:));
leafOrder = optimalleaforder(tree,D)
zz=zscore(FRunitsPSTHtrain(theseCells' &cellT==3,:),[],2);
imagesc(zz(isort1,isort2))
colormap(brewermap([],'*RdBu'))
ylabel('unit # (sorted by max response)')
xlabel('image #')
title('TAN')
makepretty; 
%% max 
figure();
subplot(131)
ops=struct;
%ops.nCall = [30,2];
ops.isort = [];
[isort1, isort2, Sm] = mapTmap(zscore(FRunitsPSTHtrain1(theseCells' & cellT==1,:),[],2), ops);
[val,ii]=max(zscore(FRunitsPSTHtrain1(theseCells' &cellT==1,:),[],2),[],2);
[sV, sI] =sort(ii);
tree= linkage(FRunitsPSTHtrain(cellT==1,:));
D = pdist(FRunitsPSTHtrain(cellT==1,:));
leafOrder = optimalleaforder(tree,D)
zz=zscore(FRunitsPSTHtrain2(theseCells' &cellT==1,:),[],2);
imagesc(zz(sI,:))
colormap(brewermap([],'*RdBu'))
ylabel('unit # (sorted by max response)')
xlabel('image #')
title('MSN')
makepretty; 

subplot(132)
[isort1, isort2, Sm] = mapTmap(zscore(FRunitsPSTHtrain1(theseCells' &cellT==2,:),[],2), ops);
[val,ii]=max(zscore(FRunitsPSTHtrain1(theseCells' &cellT==2,:),[],2),[],2);
[sV, sI] =sort(ii);
tree= linkage(FRunitsPSTHtrain(cellT==2,:));
D = pdist(FRunitsPSTHtrain(cellT==2,:));
leafOrder = optimalleaforder(tree,D)
zz=zscore(FRunitsPSTHtrain2(theseCells' &cellT==2,:),[],2);
imagesc(zz(sI,:))
colormap(brewermap([],'*RdBu'))
ylabel('unit # (sorted by max response)')
xlabel('image #')
title('FSI')
makepretty; 

subplot(133)
[isort1, isort2, Sm] = mapTmap(zscore(FRunitsPSTHtrain(theseCells' &cellT==3,:),[],2), ops);
[val,ii]=max(zscore(FRunitsPSTHtrain1(theseCells' &cellT==3,:),[],2),[],2);
[sV, sI] =sort(ii);
tree= linkage(FRunitsPSTHtrain(cellT==3,:));
D = pdist(FRunitsPSTHtrain(cellT==3,:));
leafOrder = optimalleaforder(tree,D)
zz=zscore(FRunitsPSTHtrain2(theseCells' &cellT==3,:),[],2);
imagesc(zz(sI,:))
colormap(brewermap([],'*RdBu'))
ylabel('unit # (sorted by max response)')
xlabel('image #')
title('TAN')
makepretty; 
%% tSNE 
Y=tsne(zscore(FRunitsPSTHtrain(theseCells' & cellT==1,:),[],2));
figure();
subplot(131)
scatter(Y(:, 1), Y(:, 2),7, 'r', 'filled'); hold on;
[ss,ssi]=sort(Y(:,1));
ylabel('dim 2')
xlabel('dim 1')
title('tSNE')
makepretty;
subplot(132)
zz=zscore(FRunitsPSTHtrain(theseCells' & cellT==1,:),[],2);
imagesc(zz(ssi,:))

colormap(brewermap([],'*RdBu'))
[ss,ssi]=sort(Y(:,2));
ylabel('neuron (sorted by dim 1)')
xlabel('stim #')
makepretty;
subplot(133)
zz=zscore(FRunitsPSTHtrain(theseCells' & cellT==1,:),[],2);
imagesc(zz(ssi,:))
ylabel('neuron (sorted by dim 2)')
xlabel('stim #')
makepretty;
colormap(brewermap([],'*RdBu'))

ff=find(theseCells');
Y=tsne(zscore(FRunitsPSTHtrain(theseCells',:),[],2));
figure();
scatter(Y(cellT(ff)==1, 1), Y(cellT(ff)==1, 2),7, 'r', 'filled'); hold on;
scatter(Y(cellT(ff)==2, 1), Y(cellT(ff)==2, 2), 7,'b','filled')
%scatter(Y(cellT==3, 1), Y(cellT==3, 2), 7,'g','filled')
ylabel('Dim. 1')
xlabel('Dim. 2')
makepretty;

Y=tsne(zscore(squeeze(nanmean(FRunitsTime(:,:,:))),[],2));
figure();
scatter(Y(:, 1), Y(:, 2),7, 'r', 'filled'); hold on;
ylabel('Dim. 1')
xlabel('Dim. 2')
makepretty;
%% Shuffle test r2 

%% Selectivity in space 