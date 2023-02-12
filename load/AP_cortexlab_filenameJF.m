function [filename,file_exists] = AP_cortexlab_filenameJF(animal,date,experiment,file,site,recording, shank)
myPaths;
% [filename,file_exists] = AP_cortexlab_filename(animal,date,experiment,file,site,recording)
%
% This is an absolute mess because of lab-wide inconsistency
%
% file types:
% expInfo
% timeline
% block
% parameters
% protocol
% eyecam
% eyecam_processed
% facecam
% facecam_processed
% facecam_movement
% hardware
% imaging
% ephys
% ephys_ks1
% ephys_dir
% ephys_ap

% Make inputs strings if they're numbers
if isnumeric(experiment)
    experiment = num2str(experiment);
end

if isnumeric(date)
    experiment = num2str(date);
end

% Site = multiple probes
if exist('site','var') && ~isempty(site)
    if isnumeric(site)
        site = num2str(site);
    end
    site_dir = [filesep 'site' site];
else
    site = [];
    site_dir = [];
end

% Recording = separated recordings from single probe
if exist('recording','var') && ~isempty(recording)
    if isnumeric(recording)
        recording = num2str(recording);
    end
    recording_dir = [filesep 'experiment' recording];
else
    recording = [];
    recording_dir = [];
end

% List servers
server1 = zinuPath;
server2 = zaruPath;
server3 = znasPath;
server4 = zserverPath;
server5 = localExtHdPath;
% Check that servers are accessible (login needed on restart)
% if ~exist([server1])
%     error('Zserver not available');
% end
% if ~exist([server2])
%     error('Zubjects not available');
% end
% if ~exist([server3])
%     error('Znas not available');
% end

% List all folders to check
server_location = cell(0);
server_location{end+1} = [server3 ];
server_location{end+1} = [server2 ];
server_location{end+1} = [server1 ];
server_location{end+1} = [server4 ];
server_location{end+1} = [server5 ];


switch file
    
    case 'expInfo'
        filepattern = [animal filesep date filesep];
        [filename,file_exists] = check_locations(filepattern,server_location);
        if file_exists
            filename = fileparts(filename{1});
        end
    case 'syncMess'
        filepattern = [animal filesep date filesep ...
                'ephys' site_dir filesep 'experiment*' filesep 'recording*' ...
                filesep 'sync_messages.txt'];
            [filename,file_exists] = check_locations(filepattern,server_location);
    case 'histo'
        filepattern = [animal filesep '*istology/processed/slices/probe_ccf.mat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        if ~file_exists
             filepattern = [animal filesep '*istology/025_micron/brainreg/probe_ccf.mat'];
            [filename,file_exists] = check_locations(filepattern,server_location);
        end
        if ~file_exists
             filepattern = [animal filesep '*istology/slices/probe_ccf.mat'];
            [filename,file_exists] = check_locations(filepattern,server_location);
        end
    case 'acuteRecInfo'
        filepattern = [animal filesep 'Acute_rec_' animal '.csv'];
        [filename,file_exists] = check_locations(filepattern,server_location);
    case 'probe2ephys'
        filepattern = [animal filesep 'Histology/processed/probe2ephys.mat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        if ~file_exists
             filepattern = [animal filesep '*istology/025_micron/brainreg/probe2ephys.mat'];
            [filename,file_exists] = check_locations(filepattern,server_location);
        end
    case 'mainfolder'
        filepattern = [animal filesep date];
        [filename,file_exists] = check_locations(filepattern,server_location);
    case 'timeline'
        filepattern = [animal filesep date filesep experiment ...
            filesep date '_' experiment '_' animal '_Timeline.mat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        if ~file_exists %server down, had to use 'default'
             filepattern = [animal filesep date filesep experiment ...
            filesep date '_' experiment '_default_Timeline.mat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        end
        
    case 'block'
        filepattern = [animal filesep date filesep experiment ...
            filesep date '_' experiment '_' animal '_Block.mat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        if ~file_exists %server down, had to use 'default'
            filepattern = [animal filesep date filesep experiment ...
            filesep date '_' experiment '_default_Block.mat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        end
            
    case 'parameters'
        filepattern = [animal filesep date filesep experiment ...
            filesep date '_' experiment '_' animal '_parameters.mat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        if ~file_exists %server down, had to use 'default'
            filepattern = [animal filesep date filesep experiment ...
            filesep date '_' experiment '_default_parameters.mat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        end
    case 'protocol'
        filepattern = [animal filesep date filesep experiment filesep 'Protocol.mat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        
        
    case 'eyecam'
        filepattern = [animal filesep date filesep experiment filesep 'eye.mj2'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        
    case 'facecam'
        filepattern = [animal filesep date filesep experiment filesep 'face.mj2'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        
    case 'eyecam_processed'
        filepattern = [animal filesep date filesep experiment filesep 'eye_proc.mat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        
    case 'facecam_processed'
        filepattern = [animal filesep date filesep experiment filesep 'face_proc.mat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        
    case 'facecam_movement'
        % (output from AP_mouse_movie_movement)
        filepattern = [animal filesep date filesep experiment filesep 'facecam_movement.mat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        
    case 'hardware'
        filepattern = [animal filesep date filesep experiment ...
            filesep date '_' experiment '_' animal '_hardwareInfo.mat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        if ~file_exists %server down, had to use 'default'
            filepattern = [animal filesep date filesep experiment ...
            filesep date '_' experiment '_default_hardwareInfo.mat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        end
    case 'imaging'
        filepattern = [animal filesep date filesep filesep 'svd*'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        if file_exists && iscell(filename)
            filename = fileparts(filename{1});
        elseif file_exists && isstr(filename)
            filename = fileparts(filename);
        end
        
    case 'ephys_dir'
        % (the path where the ephys data is kept)
        filepattern = [animal filesep date filesep 'ephys' site_dir];
        [filename,file_exists] = check_locations(filepattern,server_location);
        if ~file_exists
            filepattern = [animal filesep date filesep 'ephys' filesep 'recording' num2str(recording) site_dir];
        [filename,file_exists] = check_locations(filepattern,server_location);
      
        end
        if file_exists
            filename = fileparts(filename{1});
        end
    case 'ephys_parentDir'
        % (the path where the ephys data is kept)
        filepattern = [animal filesep date filesep 'ephys'];
        [filename,file_exists] = check_locations(filepattern,server_location);
       
        if file_exists
            filename = fileparts(filename{1});
        end
        
        
        
    case 'ephys_ap'
        % (the raw action potential band data file)
        
        % Old open ephys
        filepattern = [animal filesep date filesep 'ephys' site_dir filesep 'experiment*_10*-0_0.dat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        
        % New open ephys
        if ~file_exists
            filepattern = [animal filesep date filesep ...
                'ephys' site_dir filesep 'experiment*' filesep 'recording*' ...
                filesep 'continuous'  filesep 'Neuropix-3a-100.0' filesep 'continuous.dat'];
            [filename,file_exists] = check_locations(filepattern,server_location);
        end
        if ~file_exists
            filepattern = [animal filesep date filesep ...
                'ephys' site_dir filesep 'recording*' ...
                filesep 'continuous'  filesep 'Neuropix-3a-100.0' filesep 'continuous.dat'];
            [filename,file_exists] = check_locations(filepattern,server_location);
        end
        if ~file_exists
            filepattern = [animal filesep date filesep ...
                'ephys' site_dir  ...
                filesep 'continuous'  filesep 'Neuropix-3a-100.0' filesep 'continuous.dat'];
            [filename,file_exists] = check_locations(filepattern,server_location);
        end
        %spike glx
        if ~file_exists
            filepattern = [animal filesep date filesep ...
                'ephys' site_dir filesep recording_dir filesep '*ap.bin'];
            [filename,file_exists] = check_locations(filepattern,server_location);
        end
        if ~file_exists
            filepattern = [animal filesep date filesep ...
                'ephys' site_dir filesep 'experiment*' filesep '*ap.bin'];
            [filename,file_exists] = check_locations(filepattern,server_location);
        end
        if ~file_exists
            filepattern = [animal filesep date filesep ...
                'ephys' site_dir filesep 'experiment*' filesep 'recording' num2str(recording) filesep '*ap.bin'];
            [filename,file_exists] = check_locations(filepattern,server_location);
        end
        if ~file_exists
            filepattern = [animal filesep date filesep ...
                'ephys'  filesep  'recording' num2str(recording)  site_dir filesep '*ap.bin'];
            [filename,file_exists] = check_locations(filepattern,server_location);
        end
        %spike_glx local 
        if ~file_exists
            filepattern = [animal filesep date filesep ...
                 site_dir filesep '*ap.bin'];
            [filename,file_exists] = check_locations(filepattern,server_location);
        end
    case 'ephys_includingCompressed'
         % New open ephys
        
            filepattern = [animal filesep date filesep ...
                'ephys' site_dir filesep 'experiment*' filesep 'recording*' ...
                filesep 'continuous'  filesep 'Neuropix-3a-100.0' filesep 'continuous.dat'];
            [filename,file_exists] = check_locations(filepattern,server_location);
        
        if ~file_exists
            filepattern = [animal filesep date filesep ...
                'ephys' site_dir filesep 'recording*' ...
                filesep 'continuous'  filesep 'Neuropix-3a-100.0' filesep 'continuous.dat'];
            [filename,file_exists] = check_locations(filepattern,server_location);
        end
        if ~file_exists
            filepattern = [animal filesep date filesep ...
                'ephys' site_dir  ...
                filesep 'continuous'  filesep 'Neuropix-3a-100.0' filesep 'continuous.dat'];
            [filename,file_exists] = check_locations(filepattern,server_location);
        end
        
          if ~file_exists
             
            filepattern = [animal filesep date filesep ...
                'ephys' site_dir filesep recording_dir filesep '*ap.*bin'];
            [filename,file_exists] = check_locations(filepattern,server_location);
          end
     
        if ~file_exists
            filepattern = [animal filesep date filesep ...
                'ephys' site_dir filesep 'experiment*' filesep '*ap.*bin'];
            [filename,file_exists] = check_locations(filepattern,server_location);
        end
        if ~file_exists
            filepattern = [animal filesep date filesep ...
                'ephys' site_dir filesep 'experiment*' filesep 'recording' num2str(recording) filesep '*ap.bin'];
            [filename,file_exists] = check_locations(filepattern,server_location);
        end
        if ~file_exists
            filepattern = [animal filesep date filesep ...
                'ephys'  filesep  'recording' num2str(recording)  site_dir filesep '*ap.*bin'];
            [filename,file_exists] = check_locations(filepattern,server_location);
        end
        %spike_glx local 
        if ~file_exists
            filepattern = [animal filesep date filesep ...
                 'ephys' filesep site_dir filesep '*ap.*bin'];
            [filename,file_exists] = check_locations(filepattern,server_location);
        end
     case 'ephys_histology'
        
            site_dir = ['site' num2str(site), '-*' num2str(shank-1), '*'];
            filepattern = [animal filesep date filesep ...
                'ephys' filesep site_dir filesep recording_dir filesep '*ap.*bin'];
            [filename,file_exists] = check_locations(filepattern,server_location);
     
        if ~file_exists
            filepattern = [animal filesep date filesep ...
                'ephys' filesep site_dir filesep 'experiment*' filesep '*ap.*bin'];
            [filename,file_exists] = check_locations(filepattern,server_location);
        end
        if ~file_exists
            filepattern = [animal filesep date filesep ...
                'ephys' filesep site_dir filesep 'experiment*' filesep 'recording' num2str(recording) filesep '*ap.bin'];
            [filename,file_exists] = check_locations(filepattern,server_location);
        end
        if ~file_exists
            filepattern = [animal filesep date filesep ...
                'ephys'  filesep  'recording' num2str(recording) filesep site_dir filesep '*ap.*bin'];
            [filename,file_exists] = check_locations(filepattern,server_location);
        end
        %spike_glx local 
        if ~file_exists
            filepattern = [animal filesep date filesep ...
                 'ephys' filesep site_dir filesep '*ap.*bin'];
            [filename,file_exists] = check_locations(filepattern,server_location);
        end
    case 'ephys'
        % (folder with kilosort/phy outputs)
        
        kilosort_version = 2; % (kilosort 2 by default)
        
        % Drop the kilosort version in the base workspace
        assignin('base','kilosort_version',kilosort_version);
         filepattern = [animal filesep date filesep ...
                'ephys'  filesep  'kilosort2' filesep 'recording' num2str(recording)  site_dir ];
            [filename,file_exists] = check_locations(filepattern,server_location);
         
            if ~file_exists
        filepattern = [animal filesep date filesep 'ephys'  filesep 'kilosort2' filesep site_dir filesep recording_dir];
        [filename,file_exists] = check_locations(filepattern,server_location);
            end
         if ~file_exists
             filepattern = [animal filesep date filesep 'ephys'  filesep 'kilosort2' filesep site_dir filesep 'experiment' num2str(experiment) filesep recording_dir];
        [filename,file_exists] = check_locations(filepattern,server_location);
        
         end
         if ~file_exists
             filepattern = [animal filesep date filesep 'ephys'  filesep 'kilosort2' filesep site_dir   ];
        [filename,file_exists] = check_locations(filepattern,server_location);
        
         end
         if ~file_exists
             filepattern = [animal filesep date filesep 'ephys'  filesep 'kilosort3' filesep site_dir   ];
        [filename,file_exists] = check_locations(filepattern,server_location);
        
         end
         


          if ~file_exists
             filepattern = [animal filesep date filesep 'ephys'  filesep 'pykilosort' filesep site_dir filesep 'output' ];
        [filename,file_exists] = check_locations(filepattern,server_location);
        
          end
         

        
        if file_exists
            filename = fileparts(filename{end});
        end
       
        
    case 'ephys_ks1'
        % folder with kilosort/phy outputs
        
        kilosort_version = 1; % (kilosort 2 by default)
        
        % Drop the kilosort version in the base workspace
        assignin('base','kilosort_version',kilosort_version);
        
        filepattern = [animal filesep date filesep 'ephys' filesep 'kilosort' filesep site_dir filesep recording_dir];
        [filename,file_exists] = check_locations(filepattern,server_location);
        if file_exists
            filename = fileparts(filename{1});
        end
        
    case 'probe_ccf'
        % Histology probe location from AP histology
        % (sometimes upper/lowecase "Histology" folder)
        filepattern = [animal filesep '*istology' filesep 'slices' filesep 'probe_ccf.mat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        
end
end


function [filename,file_exists] = check_locations(filepattern,server_location)

% Loop through all server locations and look for file
for curr_location = 1:length(server_location)
    curr_filename = [server_location{curr_location} filesep filepattern];
    curr_filename_dir = dir(curr_filename);
    file_exists = ~isempty(curr_filename_dir);
    
    if file_exists
        % If found, break and use this filename
        if length(curr_filename_dir) == 1
            filename = [curr_filename_dir.folder filesep curr_filename_dir.name];
        else
            filename = cellfun(@(path,fn) [path filesep fn], ...
                {curr_filename_dir.folder},{curr_filename_dir.name},'uni',false);
        end
        break
    else
        % If not found, clear filename
        filename = [];
    end
end

end





