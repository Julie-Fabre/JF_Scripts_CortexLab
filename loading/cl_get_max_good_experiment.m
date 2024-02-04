function goodExp_max_trials = cl_get_max_good_experiment(experiments_struct, recording, site, exp_idx, thisDate_idx, animal, protocol)

goodExp_info_all = table();

for iExperiment = 1:size(exp_idx, 2)
    experiment = exp_idx(iExperiment); %for now, just load 1 experiment
    thisDate = experiments_struct(thisDate_idx).thisDate;

    [goodExpInfo_name, goodExpInfo_exists] = cl_cortexlab_filename(animal, thisDate, experiment, 'goodExpInfo', site, recording);

    if goodExpInfo_exists
        load(goodExpInfo_name)
    else
        goodExpInfo = struct;

        load_parts.ephys = true;
        load_parts.cam = true;
        verbose = false;
        try
            cl_load_experiment;
        catch
            goodExpInfo.n_trials = 0;
            goodExpInfo.n_keep_trials = 0;
            goodExpInfo.facecam_exists = 0;
            goodExpInfo.bad_facecam = 1;
            goodExpInfo.eyecam_exists = 0;
            goodExpInfo.bad_eyecam = 1;
            goodExpInfo.error_sync = 1;
            goodExpInfo.bad_flipper = 1;
            goodExpInfo.expDef = expDef;
        end
        
        if exist('no_move_trials', 'var')% only in passive protocols
            keep_trials = no_move_trials;
        else
            keep_trials = ones(n_trials(end),1);
        end

        goodExpInfo.n_trials = length(keep_trials);
        goodExpInfo.n_keep_trials = sum(keep_trials);
        goodExpInfo.facecam_exists = facecam_exists;
        goodExpInfo.bad_facecam = bad_facecam;
        goodExpInfo.eyecam_exists = eyecam_exists;
        goodExpInfo.bad_eyecam = bad_eyecam;
        goodExpInfo.error_sync = error_sync;
        goodExpInfo.bad_flipper = bad_flipper;
        goodExpInfo.expDef = expDef;

        expInfo_name = cl_cortexlab_filename(animal, thisDate, experiment, 'expInfo', site, recording);

        save([expInfo_name filesep, num2str(experiment), ...
            filesep, ['goodExpInfo_', num2str(site), '_' num2str(recording) '.mat']], 'goodExpInfo')

    end
    
    % select bext experiment: correct protocol, no errors (or at least no
    % ephys ones) and max. n trials.
    goodExp_info_all.n_trials(iExperiment) = goodExpInfo.n_trials;
    goodExp_info_all.n_keep_trials(iExperiment)  = goodExpInfo.n_keep_trials ;
    goodExp_info_all.facecam_exists(iExperiment) = goodExpInfo.facecam_exists;
    goodExp_info_all.bad_facecam(iExperiment) = goodExpInfo.bad_facecam;
    goodExp_info_all.eyecam_exists(iExperiment) = goodExpInfo.eyecam_exists;
    goodExp_info_all.bad_eyecam(iExperiment) = goodExpInfo.bad_eyecam;
    goodExp_info_all.error_sync(iExperiment) = goodExpInfo.error_sync;
    goodExp_info_all.bad_flipper(iExperiment) = goodExpInfo.bad_flipper;
    goodExp_info_all.expDef{iExperiment} = goodExpInfo.expDef;


end



% find experiments with no pbs 
goodExps = find(goodExp_info_all.facecam_exists == 1 & goodExp_info_all.eyecam_exists == 1 ... %there are vids
    & goodExp_info_all.bad_facecam == 0 & goodExp_info_all.bad_eyecam == 0 ...%those vids have no issues
    & goodExp_info_all.error_sync == 0 & goodExp_info_all.bad_flipper == 0 ...% no problems with ephys alignement
    & contains(goodExp_info_all.expDef, protocol));% correct expdef

if isempty(goodExps)
    % if not, find experiments with no timeline/ephys pbs 
    goodExps = find(goodExp_info_all.error_sync == 0 & goodExp_info_all.bad_flipper == 0 ...% no problems with ephys alignement
        & contains(goodExp_info_all.expDef, protocol));
    if isempty(goodExps) % if there are no experiments with no timeline/ephys pbs, return no exps.
        goodExp_max_trials = [];
    else
        goodExp_max_trials = exp_idx(goodExps(find(goodExp_info_all.n_keep_trials(goodExps) == max(goodExp_info_all.n_keep_trials(goodExps)))));
    end
else
    goodExp_max_trials = goodExps(find(goodExp_info_all.n_keep_trials(goodExps) == max(goodExp_info_all.n_keep_trials(goodExps))));
end

% of those, get the experiments with the maximum number of trials
goodExp_max_trials = goodExps(find(goodExp_info_all.n_keep_trials(goodExps) == max(goodExp_info_all.n_keep_trials(goodExps))));


end