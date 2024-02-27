
MiceOpt = {'JF067'}; % Add all mice you want to analyze
SaveDir = '/home/netshare/zinu/';
%% Parameters and settings
PrepareClusInfoparams.RunPyKSChronicStitched = 0; % Default 0. if 1, run PyKS chronic recordings stitched when same IMRO table was used
PrepareClusInfoparams.CopyToTmpFirst = 1; % If 1, copy data to local first, don't run from server (= advised!)
PrepareClusInfoparams.DecompressLocal = 1; % If 1, uncompress data first if it's currently compressed (= necessary for unitmatch and faster for QualityMetrics)
PrepareClusInfoparams.deNoise = 0;
PrepareClusInfoparams.nSavedChans = 385;
PrepareClusInfoparams.nSyncChans = 1;

% Storing preprocessed data?
PrepareClusInfoparams.ReLoadAlways = 1; % If 1, SP & Clusinfo are always loaded from KS output
PrepareClusInfoparams.saveSp = 1; % Save SP struct for easy loading of preprocessed data
PrepareClusInfoparams.binsz = 0.01; %Bin size for PSTHs in seconds
PrepareClusInfoparams.deNoise = 1; %whether to try and denoise data

% Quality Metrics
PrepareClusInfoparams.RunQualityMetrics = 1; % If 1, Run the quality metrics (Bombcell @JulieFabre)
PrepareClusInfoparams.RedoQM = 0; %if 1, redo quality metrics if it already exists
PrepareClusInfoparams.InspectQualityMetrics = 0; % If 1, Inspect the quality matrix/data set using the GUI (manual inspection)
PrepareClusInfoparams.loadPCs = 0; % Only necessary when computiong isoluation metrics/drift in QM. You save a lot of time keeping this at 0

% UnitMatch
PrepareClusInfoparams.UnitMatch = 1; % If 1, find identical units across sessions or oversplits in a fast and flexible way
PrepareClusInfoparams.RedoUnitMatch = 1; % if 1, Redo unitmatch
PrepareClusInfoparams.separateIMRO = 0; % Run for every IMRO separately (for memory reasons or when having multiple probes this might be a good idea)

% UnitMatch Parameters:
% All parameters to choose from: {'AmplitudeSim','spatialdecaySim','WavformMSE','WVCorr','CentroidDist','CentroidVar','CentroidDistRecentered','TrajAngleSim','TrajDistSim'};
% WavformSim is average of WVCorr and WavformMSE
% CentroidOverlord is average of CentroidDistRecentered and CentroidVar
PrepareClusInfoparams.Scores2Include = {'CentroidDist', 'WavformSim', 'CentroidOverlord', 'spatialdecaySim', 'AmplitudeSim', 'TrajAngleSim'}; %{'AmplitudeSim','spatialdecayfitSim','WavformSim','CentroidDist','CentroidVar','TrajAngleSim'}; %
PrepareClusInfoparams.ApplyExistingBayesModel = 0; %If 1, use probability distributions made available by us -
PrepareClusInfoparams.MakePlotsOfPairs = 0; % Plots pairs for inspection (UnitMatch)
PrepareClusInfoparams.AssignUniqueID = 1; % Assign UniqueID
PrepareClusInfoparams.GoodUnitsOnly = 1; % Include only good untis in the UnitMatch analysis - faster and more sensical
PrepareClusInfoparams.extractSync = 0;
PrepareClusInfoparams.plot = 0;


PrepareClusInfoparams.loadMATsToSave = 1;
SummarizeAcrossMice_JF;