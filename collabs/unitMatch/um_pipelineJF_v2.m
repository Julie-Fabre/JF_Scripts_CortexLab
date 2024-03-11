
%% Where to save data - CHANGE THESE PATHS
clear all;
close all;
cl_myPaths;
tmpdatafolder = extraHDPath; %'/media/julie/ExtraHD/data_temp'; % temporary folder for temporary decompression of data

%% Information on mice and recording types - CHANGE THESE MOUSE NAMES AND RECORDING TYPES
MiceOpt = {'JF067', 'JF084', 'JF078', 'JF082','JFAL035'}; % Add all mice you want to analyze. 51 = ventricule mostly, don't include
RecordingType(ismember(MiceOpt, {'JF067', 'JF078', 'JF_AL035', 'JF082', 'JF084'})) = {'Chronic'}; % 'Acute' or 'Chronic'
miceDays = [1, 16;  3, 15; 1, 10;  1, 12; 1, 20]; % 1:11
sites = [1, 1; 2, 4; 1, 3; 1, 2; 1, 1];
runMe = 1;
saveJF = 1;
for iMouse = 1:size(MiceOpt, 2)

    %% get all raw ephys and kilosort directories - CHANGE THESE PATHS
    theseSites = sites(iMouse, 1):sites(iMouse, 2);
    for iSite = 1:size(theseSites, 2)
        site = theseSites(iSite);
        experiments = cl_find_experiments(MiceOpt{iMouse}, '', true, site); % find all experiments for this mouse
        experiments = experiments([experiments.ephys]); % keep only experiments that have ephys

        ephys_dirs = {experiments.ephys_raw_paths};
        kilosort_dirs = arrayfun(@(x) {[experiments(x).ephys_ks_paths, filesep, 'site', num2str(site)]}, 1:size(experiments, 1));
        raw_data_paths = arrayfun(@(x) {experiments(x).ephys_raw_paths}, 1:size(experiments, 1));

        thisRecordingType = RecordingType{iMouse};
        mouseName = MiceOpt{iMouse};

        [filename, file_exists] = cl_cortexlab_filename(MiceOpt{iMouse}, experiments(end).day, '', 'histo_folder', '', '', '');
        SaveDir = [fileparts(filename), filesep, 'UnitMatch', filesep, 'site', num2str(site)];

        if runMe

            %% get location info
            ephys_dirs = ephys_dirs(miceDays(iMouse, 1):miceDays(iMouse, 2));
            kilosort_dirs = kilosort_dirs(miceDays(iMouse, 1):miceDays(iMouse, 2));
            raw_data_paths = raw_data_paths(miceDays(iMouse, 1):miceDays(iMouse, 2));
            day_str = {experiments.day};
            day_str = day_str(miceDays(iMouse, 1):miceDays(iMouse, 2));

            % remove any empty ones
            nonEmptyCells = ~cellfun(@isempty, ephys_dirs) & ~cellfun(@isempty, kilosort_dirs) & ...
                arrayfun(@(x) length(kilosort_dirs{x}) > 6, 1:size(kilosort_dirs, 2));
            ephys_dirs = ephys_dirs(nonEmptyCells);
            kilosort_dirs = kilosort_dirs(nonEmptyCells);
            raw_data_paths = raw_data_paths(nonEmptyCells);
            day_str = day_str(nonEmptyCells);

            % store in UMparam structure
            UMparam.SaveDir = SaveDir;
            UMparam.KSDir = kilosort_dirs;
            UMparam.RawDataPaths = raw_data_paths;

            %% pre-process
            set(0, 'DefaultFigureVisible', 'off')
            UMparam = ExtractKilosortData(UMparam.KSDir, UMparam, UMparam.RawDataPaths); % Extract KS data and do some noise removal, optionally decompresses cbin to bin data and uses BOMBCELL quality metric to define good single units
            clusinfo = getClusinfo(UMparam.KSDir); % prepare clusinfo struct

            set(0, 'DefaultFigureVisible', 'on')

            %% Load default parameters
            UMparam = DefaultParametersUnitMatch(UMparam);

            %% UnitMatch algorithm:
            [UniqueIDConversion, MatchTable, WaveformInfo, UMparam] = UnitMatch(clusinfo, UMparam);
            if UMparam.AssignUniqueID
                [UniqueIDConversion, MatchTable] = AssignUniqueID(UMparam.SaveDir);
            end

            %% Visualization
            PlotUnitsOnProbe(clusinfo, UMparam, UniqueIDConversion, WaveformInfo)

            %% Automatic evaluation:
            EvaluatingUnitMatch(UMparam.SaveDir); % Within session cross-validation
            ComputeFunctionalScores(UMparam.SaveDir, 1) % Only works when having access to Kilosort output (e.g. spike times etc.)

            %% Curation:
            if UMparam.MakePlotsOfPairs
                DrawPairsUnitMatch(UMparam.SaveDir);
                if UMparam.GUI
                    FigureFlick(UMparam.SaveDir)
                    pause
                end
            end

            %% Further evaluation - only works in combination with Bombcell
            QualityMetricsROCs(UMparam.SaveDir); % Only works in combination with BOMBCELL (and is added to path!!)

            %% Summary plot
            UMFiles = fullfile(UMparam.SaveDir, 'UnitMatch.mat');
            %summaryFunctionalPlots(UMFiles, 'Corr')
            % summaryFunctionalPlots_Part2(UMFiles, groupvec)
            %summaryMatchingPlots(UMFiles)

        else
            load([SaveDir, filesep, 'UnitMatch', filesep, 'UnitMatch.mat'])
        end
    end
end
