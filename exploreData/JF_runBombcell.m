function JF_runBombcell(animal, date, site, recording, rerunQM )
experiment = [];
%% run quality metrics
ephysPath = AP_cortexlab_filenameJF(animal, date, experiment, 'ephys', site, recording);
rawFile = AP_cortexlab_filenameJF(animal, date, experiment, 'ephys_ap', site, recording);

ephysDirPath = AP_cortexlab_filenameJF(animal, date, experiment, 'ephys_dir',site, recording);
if contains(rawFile, 'continuous.dat')
    metaFolder  = fileparts(fileparts(fileparts(rawFile)));
    ephysMetaDir = dir([metaFolder, '/*.oebin']);
else 
    ephysMetaDir = dir([ephysDirPath, '/*.meta']);
end
savePath = fullfile(ephysDirPath, 'qMetrics');
qMetricsExist = dir(fullfile(savePath, 'qMetric*.mat'));
if isempty(qMetricsExist) || rerunQM 
    [spikeTimes_samples, spikeTemplates, ...
        templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions] = bc_loadEphysData(ephysPath);

    bc_qualityParamValues; 
  

   bc_runAllQualityMetrics(param, spikeTimes_samples, spikeTemplates, ...
    templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions, savePath)
end

end