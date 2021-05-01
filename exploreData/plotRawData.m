 [ephys_path, protocol_exists] = AP_cortexlab_filenameJF('AP084', '2020-12-11', experiment, 'ephys_ap', 1);
        raw.n_channels = 384; % number of channels
        raw.ephys_datatype = 'int16';
        raw.kp = strcat(ephys_path(1:end) ,filesep, '..');
        raw.dat = dir([raw.kp, filesep, '*.dat']);
        raw.ephys_ap_filename = [raw.dat(1).folder, filesep, raw.dat(1).name];
        raw.ap_dat_dir = dir(raw.ephys_ap_filename);
        raw.pull_spikeT = -40:41; % number of points to pull for each waveform
        raw.microVoltscaling = 0.19499999284744263; %in structure.oebin for openephys, this never changed so hard-coded here-not loading it in.

        raw.dataTypeNBytes = numel(typecast(cast(0, raw.ephys_datatype), 'uint8'));
        raw.n_samples = raw.ap_dat_dir.bytes / (raw.n_channels * raw.dataTypeNBytes);
        raw.ap_data = memmapfile(raw.ephys_ap_filename, 'Format', {raw.ephys_datatype, [raw.n_channels, raw.n_samples], 'data'});
        raw.max_pull_spikes = 200; % number of spikes to pull
        
        
 cCount=cumsum(repmat(1000,size(1:384,1),1),1);
    

   % t=timeChunkStart:timeChunkStop;
    LinePlotReducer(@plot, t, double(raw.ap_data.data.data(1:384, :))+double(cCount),'k');
   
    LinePlotExplorer(gcf);
    %overlay the spikes