
%% Where to save data - CHANGE THESE PATHS
SaveDir = '/home/netshare/zinu/JF067/'; % Folder where to store the results
tmpdatafolder = '/media/julie/ExtraHD/data_temp'; % temporary folder for temporary decompression of data 

%% Information on mice and recording types - CHANGE THESE MOUSE NAMES AND RECORDING TYPES
MiceOpt = {'JF067'}; % Add all mice you want to analyze
RecordingType(ismember(MiceOpt,{'JF067'}))={'Chronic'}; % 'Acute' or 'Chronic'

for iMouse = 1:size(MiceOpt,2)

    %% get all raw ephys and kilosort directories - CHANGE THESE PATHS
    site = 1;
    experiments = AP_find_experimentsJF(MiceOpt{iMouse}, '', true, site); % find all experiments for this mouse
    experiments = experiments([experiments.ephys]); % keep only experiments that have ephys

    ephys_dirs = {experiments.ephys_raw_paths};
    kilosort_dirs = arrayfun(@(x) {[experiments(x).ephys_ks_paths, filesep, 'site', num2str(site)]}, 1:size(experiments,1));
    thisRecordingType = RecordingType{iMouse};
    mouseName = MiceOpt{iMouse};

    %% subselect two first for testing 
    ephys_dirs =  ephys_dirs(1:29);
    kilosort_dirs =  kilosort_dirs(1:29);

    %% run unit match
    set(0,'DefaultFigureVisible','off')
    [PrepareClusInfoparams, UMparam, UniqueIDConversion, MatchTable, WaveformInfo] = um_runUnitMatch(kilosort_dirs, ephys_dirs, SaveDir, tmpdatafolder, thisRecordingType, mouseName);
    set(0,'DefaultFigureVisible','on')

    %% check output
    EvaluatingUnitMatch([SaveDir, filesep, 'UnitMatch']);
    ComputeFunctionalScores([SaveDir, filesep, 'UnitMatch'])

    DrawBlind = 0; %1 for blind drawing (for manual judging of pairs)
    DrawPairsUnitMatch([SaveDir, filesep, 'UnitMatch'], DrawBlind);
    
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
    FigureFlick(SaveDir,'julie',recompute)
    % your labels are saved in MatchTable.<user>
end