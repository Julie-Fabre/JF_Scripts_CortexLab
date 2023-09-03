
%% Where to save data - CHANGE THESE PATHS
SaveDir = '/home/netshare/zinu/JF067'; % Folder where to store the results
tmpdatafolder = '/media/julie/ExtraHD/data_temp'; % temporary folder for temporary decompression of data 

%% Information on mice and recording types - CHANGE THESE MOUSE NAMES AND RECORDING TYPES
MiceOpt = {'JF067', 'JF078', 'JF051', 'AL035', 'JF082', 'JF084'}; % Add all mice you want to analyze
RecordingType(ismember(MiceOpt, {'JF067', 'JF078', 'JF051', 'AL035', 'JF082', 'JF084'}))={'Chronic'}; % 'Acute' or 'Chronic'
miceDays = [1,29; 1,16; 1,17; 1,14; 1,31];% 1:11
% 1:16
% 1:end
% 1:end
% 1:31
% diff saving (celian run)
saveJF=1;
for iMouse = 2:size(MiceOpt,2)-1 %JF084 - saving deal with 

    %% get all raw ephys and kilosort directories - CHANGE THESE PATHS
    site = 1;
    experiments = AP_find_experimentsJF(MiceOpt{iMouse}, '', true, site); % find all experiments for this mouse
    experiments = experiments([experiments.ephys]); % keep only experiments that have ephys

    ephys_dirs = {experiments.ephys_raw_paths};
    kilosort_dirs = arrayfun(@(x) {[experiments(x).ephys_ks_paths, filesep, 'site', num2str(site)]}, 1:size(experiments,1));
    thisRecordingType = RecordingType{iMouse};
    mouseName = MiceOpt{iMouse};

    %% subselect two first for testing 
    ephys_dirs =  ephys_dirs(miceDays(iMouse,1):miceDays(iMouse,2));
    kilosort_dirs =  kilosort_dirs(miceDays(iMouse,1):miceDays(iMouse,2));

    %% run unit match
    set(0,'DefaultFigureVisible','off')
    [PrepareClusInfoparams, UMparam, UniqueIDConversion, MatchTable, WaveformInfo] = um_runUnitMatch(kilosort_dirs, ephys_dirs, SaveDir, tmpdatafolder, thisRecordingType, mouseName);
    set(0,'DefaultFigureVisible','on')

    %% check output
    EvaluatingUnitMatch([SaveDir, filesep, 'UnitMatch']);
    ComputeFunctionalScores([SaveDir, filesep, 'UnitMatch'], saveJF)

    DrawBlind = 0; %1 for blind drawing (for manual judging of pairs)
    DrawPairsUnitMatch([SaveDir, filesep, 'UnitMatch'], DrawBlind, saveJF);
    
    % Key presses:
    %   Right arrow: next pair
    %   Left arrow: previous pair
    %   Up arrow: label as match
    %   Down arrow: label as non-match
    %   o: label as I don't know (=uncurated)
    %   m: go to first uncurated pair
    %   p: select pair number
    %   s: SAVE
    recompute = 1;
    FigureFlick([SaveDir, filesep, 'UnitMatch'],'julie',recompute, saveJF)
    % your labels are saved in MatchTable.<user>
end

%% get matches across 10 first days 
thisThreshold = 0.5;
% % find cells with pair-wise matches
% all_matches = nan(max(MatchTable.RecSes1)-1, 1);
% matchCount = 1;
% for iRecording = 1:max(MatchTable.RecSes1)-1
%     clearvars match_pairID
%     %theseCandidates = MatchTable.RecSes1==iRecording & MatchTable.MatchProb > 0.5;
%     sessionMatches = find(MatchTable.RecSes1==iRecording &...
%         MatchTable.MatchProb > thisThreshold);
%     match_pairID = [double(MatchTable.ID1(sessionMatches)), double(MatchTable.ID2(sessionMatches)),...
%         double(MatchTable.UID1(sessionMatches)), double(MatchTable.UID2(sessionMatches)),...
%         MatchTable.MatchProb(sessionMatches)];
%     excludeThese = match_pairID(:,1) == ...
%         match_pairID(:,2);
%     uniqueMatches = unique(match_pairID(~excludeThese,1)); % get unique matches
%     for iMatch = 1:size(uniqueMatches, 1)
%         % find the best match QQ
%         %best_match = find(match_pairID(:,3) == max(match_pairID(match_pairID(:,1) ==...
%         %    uniqueMatches(iMatch),3)));
%         idx = ~excludeThese & match_pairID(:,1) == uniqueMatches(iMatch);
%         best_match = find(match_pairID(:,3) == match_pairID(:,4));
%         all_matches(iRecording, matchCount) = MatchTable.ID1(sessionMatches(best_match));
%         all_matches(iRecording, matchCount) = MatchTable.ID2(sessionMatches(best_match));
%         matchCount = matchCount + 1;
%     end
% end

% find all UIDs equal, ands remove the others 
excludeThese =  MatchTable.ID1 == MatchTable.ID2;
all_matches = table;
all_matches.uid1 = MatchTable.UID1(find(MatchTable.UID1(~excludeThese) == MatchTable.UID2(~excludeThese)...
    & MatchTable.MatchProb(~excludeThese) > thisThreshold));
all_matches.day1 = MatchTable.RecSes1(find(MatchTable.UID1(~excludeThese) == MatchTable.UID2(~excludeThese)...
    & MatchTable.MatchProb(~excludeThese) > thisThreshold));
all_matches.day2 = MatchTable.RecSes2(find(MatchTable.UID1(~excludeThese) == MatchTable.UID2(~excludeThese)...
    & MatchTable.MatchProb(~excludeThese) > thisThreshold));
% get unique days 
matchCounts =histc(all_matches,unique(all_matches));
find(matchCounts>1)




%% plot visual activity over time  and waveform and ACG 
%pairs = [379,2754,2772,1929,514,746,913,939,968,1676,1780,1826,2138,2159,2341,2451,2555,2749,2773,2866,2871,2922,3029];
clearvars num_recs
for iPair = 1:max(MatchTable.UID1)
    %thisPair = MatchTable.UID1 
    % find all recording + unit labels
    %excludeThese =  MatchTable.ID1 == MatchTable.ID2;
    these_indices = find(MatchTable.UID1==iPair & MatchTable.UID2==iPair);
    num_recs(iPair) = numel(unique([MatchTable.RecSes1(these_indices); MatchTable.RecSes2(these_indices)]));
end

animal = MiceOpt{iMouse};
protocol = 'JF_choiceworldStimuli'; % (this is the name of the Signals protocol)
experiments = AP_find_experimentsJF(animal, protocol, true);
experiments = experiments([experiments.ephys]);
raster_window = [-0.5, 1];
psth_bin_size = 0.001;

curr_psth_all = nan(size(thesePairs, 2), size(recs, 1), 3, 1500);
thesePairs = find(num_recs >= 4);
for iPair = 1:size(thesePairs, 2)
    thisPair = thesePairs(iPair);
    these_indices = MatchTable.UID1 == thisPair & MatchTable.UID2 == thisPair;
    recs = unique([MatchTable.RecSes1(these_indices); MatchTable.RecSes2(these_indices)]);

    for iRec = 1:size(recs, 1)
        
        % load recording
        thisRecording = recs(iRec); 
        site = 1;
        recording = [];
        day = experiments(thisRecording).day;
        experiment = experiments(thisRecording).experiment(1);

        JF_load_experiment;

        % get visually aligned activity for unit(s) 
        theseUnits = unique([MatchTable.ID1(these_indices & MatchTable.RecSes1 == thisRecording);...
            MatchTable.ID2(these_indices & MatchTable.RecSes2 == thisRecording)]);
        
        [align_group_a, align_group_b] = ismember(trial_conditions(:,2), unique(trial_conditions(:,2)));
        [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(spike_templates, spike_times_timeline, ...
            theseUnits, raster_window, psth_bin_size, stimOn_times, align_group_b);

        % store visual activity 
        curr_psth_all(iPair, iRec, 1:size(unique(trial_conditions(:,2)),1), :) = curr_psth; 
         
    end
end

