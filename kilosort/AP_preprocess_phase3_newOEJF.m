function AP_preprocess_phase3_newOEJF(animal,day, thisSite,t_range,chanMap)
% AP_preprocess_phase3(animal,day,t_range)
%
% t_range = specify data time range to use
% (now using kilosort2, putting into a 'kilosort2' folder)
%
% Copied from AP_preprocess_phase3_newOE:
% doesn't accomodate old OE, done when started multiple experiments
% (stopped and started ephys recordings, which creates multiple experiment
% folders)

%% Get paths and filenames

[ephys_path,ephys_exists] = AP_cortexlab_filenameJF(animal,day,[],'ephys_dir');
if ~exist('chanMap','var')%default to 3A probes of no channel map provided. 
    chanMap = '/home/julie/Dropbox/MATLAB/onPaths/JF_Scripts_CortexLab/kilosort/forPRBimecP3opt3.mat';
end
if ~ephys_exists
    error([animal ' ' day ': No ephys data found']);
end

save_paths = {[ephys_path filesep 'kilosort2']};
data_paths = {ephys_path};

% Check for multiple sites (assume sites are marked as site#)
data_path_dir = dir([data_paths{1} filesep 'site*']);
if ~isempty(data_path_dir)
    data_paths = cellfun(@(x) [data_paths{1} filesep x],{data_path_dir.name},'uni',false);
    save_paths = cellfun(@(x) [save_paths{1} filesep x],{data_path_dir.name},'uni',false);
end
if ~exist('thisSite','var')
    thisSite = 1:length(data_paths);
end
for curr_site = thisSite 
    
    % Get experiments (if turned off between)
    curr_data_path = data_paths{curr_site};   
    ephys_exp_paths = dir([curr_data_path filesep 'experiment*']);
    if isempty(ephys_exp_paths)
        ephysType = 2;
        ephys_exp_paths = dir([curr_data_path filesep '*g0']);
    else
        ephysType=1;
    end
       
    for curr_exp = 1:length(ephys_exp_paths)
        
        % Update save path with experiment (only if more than one, legacy)
        if length(ephys_exp_paths) == 1
            curr_save_path = save_paths{curr_site};
        elseif length(ephys_exp_paths) > 1
            curr_save_path = [save_paths{curr_site} filesep ephys_exp_paths(curr_exp).name];
        end
        
        % Make save path
        if ~exist(curr_save_path,'dir')
            mkdir(curr_save_path)
        end
              
        % Get OE filenames (check for multiple experiments, do separately)
        if ephysType == 2
                exp_rec_dir = [ ephys_exp_paths(curr_exp).name];
        if ~isempty(t_range)
            ap_data_filename = dir([curr_data_path filesep exp_rec_dir filesep '*tRange*.bin']);

        else
             ap_data_filename = dir([curr_data_path filesep exp_rec_dir filesep '*.bin']);
        end
        lfp_data_filename =  dir([curr_data_path filesep exp_rec_dir filesep '*.bin']);
       % sync_filename = [curr_data_path filesep filesep exp_rec_dir filesep 'events' filesep 'Neuropix-3a-100.0' filesep 'TTL_1' filesep 'channel_states.npy' ];
       % sync_timestamps_filename = [curr_data_path filesep exp_rec_dir filesep 'events' filesep 'Neuropix-3a-100.0' filesep 'TTL_1' filesep 'timestamps.npy' ];
       % messages_filename = [curr_data_path filesep exp_rec_dir filesep 'sync_messages.txt'];
        settings_filename = dir([curr_data_path filesep exp_rec_dir filesep '*.meta']);
        
        else
            
        exp_rec_dir = [ ephys_exp_paths(curr_exp).name filesep 'recording1'];
        if ~isempty(t_range)

                ap_data_filename = [curr_data_path filesep exp_rec_dir filesep 'continuous' filesep 'Neuropix-3a-100.0' filesep 'continuous*tRange*.dat'];
        
  
        else
             ap_data_filename = [curr_data_path filesep exp_rec_dir filesep 'continuous' filesep 'Neuropix-3a-100.0' filesep 'continuous.dat'];
        
        end
        lfp_data_filename = [curr_data_path filesep exp_rec_dir filesep 'continuous' filesep 'Neuropix-3a-100.1' filesep 'continuous.dat'];
        sync_filename = [curr_data_path filesep filesep exp_rec_dir filesep 'events' filesep 'Neuropix-3a-100.0' filesep 'TTL_1' filesep 'channel_states.npy' ];
        sync_timestamps_filename = [curr_data_path filesep exp_rec_dir filesep 'events' filesep 'Neuropix-3a-100.0' filesep 'TTL_1' filesep 'timestamps.npy' ];
        messages_filename = [curr_data_path filesep exp_rec_dir filesep 'sync_messages.txt'];
        settings_filename = [curr_data_path filesep exp_rec_dir filesep 'structure.oebin'];
        end
        
        %% Get and save recording parameters
        if ephysType==1
        % The gains and filter cuts aren't recorded anymore?!
        ap_gain = {500};
        lfp_gain = {125};
        filter_cut = {300};
        
        % (0.195x for int16 to uV? how's this change with gain, just another x?)
        
        % Hard-coded parameters
        n_channels = 384;
        ap_sample_rate = 30000;
        lfp_sample_rate = 2500;
        
        params = {'raw_path',['''' curr_data_path '''']; ...
            'n_channels',num2str(n_channels); ...
            'ap_sample_rate',num2str(ap_sample_rate); ... % this should be 30000 AP, 2500 LFP
            'lfp_sample_rate',num2str(lfp_sample_rate);
            'ap_gain',num2str(ap_gain{1}); ...
            'lfp_gain',num2str(lfp_gain{1})
            'filter_cutoff',num2str(filter_cut{1})};
        
        param_filename = [curr_save_path filesep 'dat_params.txt'];
        
        formatSpec = '%s = %s\r\n';
        fid = fopen(param_filename,'w');
        for curr_param = 1:size(params,1)
            fprintf(fid,formatSpec,params{curr_param,:});
        end
        fclose(fid);
        
        
        %% Get/save digital input events
        
        % Get processor start information 
        % (assume 1 = software, 2 = AP, 3 = LFP)
        % This way of getting start times is garbage too, assume AP
        messages_id = fopen(messages_filename);
        messages_text = textscan(messages_id,'%*s %d@%dHz','delimiter',{'time: '});
        fclose(messages_id);
        
        % Processor start time: to be subtracted from sync so start = 0
        % (this is usually 0 but not if hit record directly from preview,
        % so either does nothing or corrects offset)
        start_time_sec = double(messages_text{1}(2)/messages_text{2}(2));
        
        % Get/save digital input event times,
        sync_data = readNPY(sync_filename);
        sync_timestamps_int = readNPY(sync_timestamps_filename);
        sync_timestamps = double(sync_timestamps_int)/ap_sample_rate;
        
        sync_channels = unique(abs(sync_data));
        sync = struct('timestamps',cell(size(sync_channels)),'values',cell(size(sync_channels)));
        for curr_sync = 1:length(sync_channels)
            sync_events = abs(sync_data) == (sync_channels(curr_sync));
            sync(curr_sync).timestamps = sync_timestamps(sync_events) - start_time_sec;
            sync(curr_sync).values = sign(sync_data(sync_events)) == 1;
        end
        
        sync_save_filename = [curr_save_path filesep 'sync.mat'];
        save(sync_save_filename,'sync');
        end
        %% Run kilosort
        
        % Set up local directory and clear out
        ssd_kilosort_path = '/media/julie/ExtraHD/data_temp/kilosort';
        
        % Clear out local kilosort directories
        if exist(ssd_kilosort_path,'dir')
            rmdir(ssd_kilosort_path,'s')
        end
        mkdir(ssd_kilosort_path);
        
        % Copy AP data locally
        disp('Copying AP data to local drive...')
        if ephysType==1
        ap_temp_filename = [ssd_kilosort_path filesep animal '_' day  '_' 'ephys_apband.dat'];
      
        else
              ap_temp_filename = [ssd_kilosort_path filesep animal '_' day  '_' 'ephys_apband.bin'];
                ap_data_filename = [curr_data_path filesep exp_rec_dir filesep ap_data_filename.name];
        end
        %copyfile(ap_data_filename,ap_temp_filename);
        disp('Done');
        
        % Clean AP data of artifacts
        disp('Cleaning AP data...')
        if ephysType==1
            ap_clean_filename = [ssd_kilosort_path filesep animal '_' day '_' 'ephys_apband_clean.dat'];
        else
            ap_clean_filename = [ssd_kilosort_path filesep animal '_' day '_' 'ephys_apband_clean.bin'];
        end
        
        if ephysType==1
            ttl_path = fileparts(sync_filename);
            %AP_clean_datJF(ap_temp_filename,n_channels,ttl_path,ap_clean_filename);
            %delete(ap_temp_filename);
        end
        % Delete local raw data
        
        
        % Kilosort 1 (old)
        %     AP_run_kilosort(ap_temp_car_filename,ap_sample_rate);
        
        % Kilosort 2
        % Set default time range to be [0,inf]
        if ~exist('t_range','var')
            t_range = [0,inf];
        end
        if ephysType==1
            dirAP = dir(ap_data_filename);
            AP_run_kilosort2JF([dirAP.folder, filesep, dirAP.name],ap_sample_rate,dirAP.folder,[],[],chanMap);
        else
            ap_sample_rate = 30000;
            AP_run_kilosort2JF(ap_temp_filename,ap_sample_rate,ssd_kilosort_path,[],[],chanMap);
        end
        %% Copy kilosort results to server
        
        disp('Copying sorted data to server...');
        ks_results_path = [ssd_kilosort_path filesep 'results'];
        copyfile(ks_results_path,curr_save_path);
        
        %% Copy kilosort results and raw data to phy folder for clustering
        
        local_phy_path = ['//media/julie/ExtraHD/data_temp/phy_staging' filesep animal '_' day];
        mkdir(local_phy_path)
        
        % Move the cleaned data into the phy directory
        [~,ap_file,ap_ext] = fileparts(ap_clean_filename);
        movefile(ap_clean_filename,[local_phy_path filesep ap_file ap_ext])
        
        % Move the results into the phy directory
        movefile([ks_results_path filesep '*'],local_phy_path)
        
        %% Delete all temporary local data
        rmdir(ssd_kilosort_path,'s');
        mkdir(ssd_kilosort_path);
        
    end
    
end

disp('Done processing phase 3 data.');


