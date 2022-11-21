function [spike_times_timeline, bad_flipper] = JF_align_ephys_to_timeline(animal, day, isSpikeGLX, flipper_flip_times_timeline, ...
    ephys_sync_folder, flipper_sync_idx, experiment_idx, acqLive_sync_idx, spike_times)
%flipper_flip_times_timelinesync(flipper_sync_idx).timestamps

%% load ephys sync (flipper)
if isSpikeGLX
    load(sprintf('%s//sync.mat', ephys_sync_folder)); % ephys flipper
    if unique(sync) < 2
        error('check flipper was inputed in correct IMEC card')
    end
    sample_rate = 30000; % IMEC sample rate.
else
    load([ephys_sync_folder, filesep, 'sync.mat']);
    if length(sync) >= flipper_sync_idx
        bad_flipper = false;
        ephys_sync_timestamps = sync(flipper_sync_idx).timestamps;
        ephys_sync_values = sync(flipper_sync_idx).values;
    else
        bad_flipper = true;
    end
end

%% find experiment times
% find long gaps in flip times = between experiments/when they end
flip_diff_thresh = 1; % time between flips to define experiment gap (s)
flipper_expt_idx = [1; find(abs(diff(ephys_sync_timestamps)) > ...
    flip_diff_thresh) + 1; length(ephys_sync_timestamps) + 1];
possibilities = diff(flipper_expt_idx);
[val, idx] = min(abs(possibilities-length(flipper_flip_times_timeline)));
if length(flipper_expt_idx) < find(experiment_idx) + 1
    experiment_idx = idx;
    flipper_flip_times_ephys = ephys_sync_timestamps( ...
        flipper_expt_idx(find(experiment_idx)):flipper_expt_idx(find(experiment_idx)+1)-1);
else
    flipper_flip_times_ephys = ephys_sync_timestamps( ...
        flipper_expt_idx(find(experiment_idx)):flipper_expt_idx(find(experiment_idx)+1)-1);
end

%% match ephys and timeline
if length(flipper_flip_times_ephys) == length(flipper_flip_times_timeline)
    % If same number of flips in ephys/timeline, use all
    sync_timeline = flipper_flip_times_timeline;
    sync_ephys = flipper_flip_times_ephys;
elseif length(flipper_flip_times_ephys) ~= length(flipper_flip_times_timeline) ...
        && val == 0
    experiment_idx = idx;
    flipper_flip_times_ephys = ephys_sync_timestamps( ...
        flipper_expt_idx(experiment_idx):flipper_expt_idx(experiment_idx+1)-1);
    sync_timeline = flipper_flip_times_timeline;
    sync_ephys = flipper_flip_times_ephys;
elseif length(flipper_flip_times_ephys) ~= length(flipper_flip_times_timeline)
    % If different number of flips in ephys/timeline, best
    % contiguous set via xcorr of diff
    warning([animal, ' ', day, ':Flipper flip times different in timeline/ephys. /n' ...
        'The fix for this is probably not robust: always check'])
    [flipper_xcorr, flipper_lags] = ...
        xcorr(diff(flipper_flip_times_timeline), diff(flipper_flip_times_ephys));
    [~, flipper_lag_idx] = max(flipper_xcorr);
    flipper_lag = flipper_lags(flipper_lag_idx);
    % (at the moment, assuming only dropped from ephys)
    sync_ephys = flipper_flip_times_ephys;
    try
        sync_timeline = flipper_flip_times_timeline(flipper_lag+1: ...
            flipper_lag+1:flipper_lag+length(flipper_flip_times_ephys));
    catch
        sync_timeline = flipper_flip_times_timeline;
    end
    if length(diff(sync_ephys)) ~= length(diff(sync_timeline))
        experiment_idx = idx;
        flipper_flip_times_ephys = ephys_sync_timestamps( ...
            flipper_expt_idx(experiment_idx):flipper_expt_idx(experiment_idx+1)-1);
        flipper_flip_times_ephys = ephys_sync_timestamps( ...
            flipper_expt_idx(idx):flipper_expt_idx(idx+1)-1);
        % If different number of flips in ephys/timeline, best
        % contiguous set via xcorr of diff
        warning([animal, ' ', day, ':Flipper flip times different in timeline/ephys. /n ' ...
            'The fix for this is probably not robust: always check'])
        [flipper_xcorr, flipper_lags] = ...
            xcorr(diff(flipper_flip_times_timeline), diff(flipper_flip_times_ephys));
        [~, flipper_lag_idx] = max(flipper_xcorr);
        flipper_lag = flipper_lags(flipper_lag_idx);
        % (at the moment, assuming only dropped from ephys)
        sync_ephys = flipper_flip_times_ephys;
        try
            sync_timeline = flipper_flip_times_timeline(flipper_lag+1: ...
                flipper_lag+1:flipper_lag+length(flipper_flip_times_ephys));
        catch
            try
                sync_timeline = flipper_flip_times_timeline(1: ...
                    1:-flipper_lag+length(flipper_flip_times_ephys));
            catch
                sync_timeline = flipper_flip_times_timeline;
            end
        end
        if length(diff(sync_ephys)) ~= length(diff(sync_timeline))
            bad_flipper = true;
        end
    end
    bad_flipper = true;
end

%% if bad flipper, use acq live 
if bad_flipper
    % (if no flipper or flipper problem, use acqLive)
    if isSpikeGLX
    else
        if length(sync) >= acqLive_sync_idx
            bad_flipper = false;
            ephys_sync_timestamps = sync(acqLive_sync_idx).timestamps;
            ephys_sync_values = sync(acqLive_sync_idx).values;
            
        else
            error([animal, ' ', day, ': bad flipper and acq live signal'])
        end
    end
    % Get acqLive times for current experiment
    experiment_ephys_starts = ephys_sync_timestamps(ephys_sync_values == max(ephys_sync_values));
    experiment_ephys_stops = ephys_sync_timestamps(ephys_sync_values == min(ephys_sync_values));
    acqlive_ephys_currexpt = [experiment_ephys_starts(idx), ...
        experiment_ephys_stops(idx)];

    sync_timeline = acqLive_timeline;
    sync_ephys = acqlive_ephys_currexpt;

    % Check that the experiment time is the same within threshold
    % (it should be almost exactly the same)
    if abs(diff(acqLive_timeline)-diff(acqlive_ephys_currexpt)) > 1
        error([animal, ' ', day, ': acqLive duration different in timeline and ephys']);
    end
end

%% Get spike times in timeline time
spike_times_timeline = interp1(sync_ephys, sync_timeline, spike_times, 'linear', 'extrap');
end