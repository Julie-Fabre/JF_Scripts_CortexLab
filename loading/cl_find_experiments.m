function experiments = cl_find_experiments(animal, protocol, flexible_name, site)
cl_myPaths;
% experiments = AP_find_experiments(animal,protocol,flexible_name)
%
% flexible_name - if true, uses strfind(lower)
%
% Find all experiments from an animal with a given protocol name
% 
% Note: file locations changed many times so this is necessarily totally
% disorganized and ridiculous

% If no protocol specified, return all experiments
if ~exist('protocol','var') || isempty(protocol)
    protocol = [];
end

% If no flexible name specified, use exact
if ~exist('flexible_name','var') || isempty(flexible_name)
    flexible_name = false;
end

% If no flexible name specified, use exact
if ~exist('site','var') || isempty(site)
    site = [];
end
% Initialize pathname, add to it with each server location
days_combined = {};
days_pathnames_combined = {};

% (look in server 1 expInfo - old)
expInfo_path = [zserverPath filesep animal];
expInfo_dir = dir(expInfo_path);
day_paths = cellfun(@(x) ~isempty(regexp(x,'\d\d\d\d-\d\d-\d\d')),{expInfo_dir.name}) &...
    [expInfo_dir.isdir];
curr_days = {expInfo_dir(day_paths).name};
curr_days_pathname = cellfun(@(x) [expInfo_path filesep x],curr_days,'uni',false);

days_combined = [days_combined,curr_days];
days_pathnames_combined = [days_pathnames_combined,curr_days_pathname];

% (look in server 1 subjects - new)
expInfo_path = [zinuPath filesep animal];
expInfo_dir = dir(expInfo_path);
day_paths = cellfun(@(x) ~isempty(regexp(x,'\d\d\d\d-\d\d-\d\d')),{expInfo_dir.name}) &...
    [expInfo_dir.isdir];
curr_days = {expInfo_dir(day_paths).name};
curr_days_pathname = cellfun(@(x) [expInfo_path filesep x],curr_days,'uni',false);

days_combined = [days_combined,curr_days];
days_pathnames_combined = [days_pathnames_combined,curr_days_pathname];

% (look in server 2)
expInfo_path = [znasPath filesep animal];
expInfo_dir = dir(expInfo_path);
day_paths = cellfun(@(x) ~isempty(regexp(x,'\d\d\d\d-\d\d-\d\d')),{expInfo_dir.name}) &...
    [expInfo_dir.isdir];
curr_days = {expInfo_dir(day_paths).name};
curr_days_pathname = cellfun(@(x) [expInfo_path filesep x],curr_days,'uni',false);

days_combined = [days_combined,curr_days];
days_pathnames_combined = [days_pathnames_combined,curr_days_pathname];


% (look in server 3)
expInfo_path = [zaruPath filesep animal];
expInfo_dir = dir(expInfo_path);
day_paths = cellfun(@(x) ~isempty(regexp(x,'\d\d\d\d-\d\d-\d\d')),{expInfo_dir.name}) &...
    [expInfo_dir.isdir];
curr_days = {expInfo_dir(day_paths).name};
curr_days_pathname = cellfun(@(x) [expInfo_path filesep x],curr_days,'uni',false);

days_combined = [days_combined,curr_days];
days_pathnames_combined = [days_pathnames_combined,curr_days_pathname];


% If multiple days: experiment number folders should be preserved across
% all folders so just pick one
[days,unique_day_idx] = unique(days_combined);
days_pathnames = days_pathnames_combined(unique_day_idx);

% Find experiments with chosen protocol and which modalities were recorded
protocol_expts = cell(size(days));
imaging_expts = cell(size(days));
ephys_expts = cell(size(days));
ephys_paths = cell(size(days));
ephys_ks_paths = cell(size(days));
ephys_raw_paths = cell(size(days));

for curr_day = 1:length(days)  
    
    day = days{curr_day};
    % Find all experiments of that day (defined as number-only folders)
    expDay_dir = dir(days_pathnames{curr_day});
    exp_folders = cellfun(@any,regexp({expDay_dir.name},'^\d*$'));
    exp_nums = cellfun(@str2num,{expDay_dir(exp_folders).name});
    
    % If looking for specific protocol, find amongst days's experiments
    if ~isempty(protocol)
        use_exp = false(size(exp_nums));
        for curr_exp = 1:length(exp_nums)
            % Check for signals or MPEP
            [block_filename, block_exists] = cl_cortexlab_filename(animal,day,exp_nums(curr_exp),'block');
            [protocol_filename,protocol_exists] = cl_cortexlab_filename(animal,day,exp_nums(curr_exp),'protocol');
            
            if block_exists
                % If signals
                try
                    load(block_filename)
                
                    
                if isfield(block,'expType') % old name
                    [~,expDef] = fileparts(block.expType);
                    if flexible_name
                        use_exp(curr_exp) = contains(lower(expDef),lower(protocol));
                    elseif ~flexible_name
                        use_exp(curr_exp) = strcmp(expDef,protocol);
                    end
                end
                if isfield(block,'expDef') % new name
                    [~,expDef] = fileparts(block.expDef);
                    if flexible_name
                        use_exp(curr_exp) = contains(lower(expDef),lower(protocol));
                    elseif ~flexible_name
                        use_exp(curr_exp) = strcmp(expDef,protocol);
                    end
                    
                end
                catch
                end
            elseif protocol_exists
                % If MPEP
                load(protocol_filename)
                [~,expDef] = fileparts(Protocol.xfile);
                if flexible_name
                    use_exp(curr_exp) = contains(lower(expDef),lower(protocol));
                elseif ~flexible_name
                    use_exp(curr_exp) = strcmp(expDef,protocol);
                end
            else
                continue
            end
        end
    else
        use_exp = true(size(exp_nums));
    end
    
    % Find days with imaging/electrophysiology
    if any(use_exp)
        protocol_expts{curr_day} = exp_nums(use_exp);
        imaging_path = cl_cortexlab_filename(animal,day,exp_nums(use_exp),'imaging');
        imaging_expts{curr_day} = exist([imaging_path filesep 'meanImage_blue.npy'],'file') > 0;
        [ephys_paths{curr_day},ephys_expts{curr_day}] = cl_cortexlab_filename(animal,day,exp_nums(use_exp),'ephys_dir');
        [ephys_ks_paths{curr_day},~] = cl_cortexlab_filename(animal,day,exp_nums(use_exp),'ephys');
        [ephys_all_raw, ~] = cl_cortexlab_filename(animal,day,exp_nums(use_exp),'ephys_includingCompressed',site);
        if size(ephys_all_raw, 2) ==2
            ephys_raw_paths{curr_day} = ephys_all_raw{1};
        else
            ephys_raw_paths{curr_day} = ephys_all_raw;
        end

    end
    
end

% Package experiment info
use_days = ~cellfun(@isempty,protocol_expts);
experiments = struct('day',cell(sum(use_days),1),'experiment',cell(sum(use_days),1),...
    'imaging',cell(sum(use_days),1),'ephys',cell(sum(use_days),1), 'ephys_paths',cell(sum(use_days),1),...
    'ephys_ks_paths',cell(sum(use_days),1), 'ephys_raw_paths',cell(sum(use_days),1));

[experiments.thisDate] = deal(days{use_days});
[experiments.experiment] = deal(protocol_expts{use_days});
[experiments.imaging] = deal(imaging_expts{use_days});
[experiments.ephys] = deal(ephys_expts{use_days});
[experiments.ephys_paths] = deal(ephys_paths{use_days});
[experiments.ephys_ks_paths] = deal(ephys_ks_paths{use_days});
[experiments.ephys_raw_paths] = deal(ephys_raw_paths{use_days});







