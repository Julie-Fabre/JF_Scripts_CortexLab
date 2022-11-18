% Check number of stim matches photodiode
% (once I saw some weird extra flip at the end of an
% experiment? so added case to only use first n flips)

% correct for incorrect photodiode flip number 
if length(signals_events.stimOnTimes) < length(stimOn_times) % too many photodiode flips
    warning([animal, ' ', day, ': different stim number signals and photodiode']);
    % do any of them violate min ITI - if so, get rid of both
    % photodiode flips that violate this
    iti_violations = diff(photodiode_flip_times(2:2:end)) < (min_ITI * timeline_sample_rate - buffer_samples);
    if any(iti_violations)
        stimOn_times(logical([0; iti_violations])) = []; % remove culprit
        stimOn_times(logical([iti_violations])) = NaN; % NaN-out the time of closest culprit too- this time can't be trusted
    end
    if length(signals_events.stimOnTimes) < length(stimOn_times) % if no violation or correction didn't work
        % else, use x first photodiode flips
        stimOn_times = stimOn_times(1:length(signals_events.stimOnTimes));
    end
elseif length(signals_events.stimOnTimes) > length(stimOn_times) % some photodiode flips missing
    warning([animal, ' ', day, ': different stim number signals and photodiode']);
    % do any of them violate min ITI - if so, get rid of both
    % photodiode flips that violate this
    iti_violations = find(diff(photodiode_flip_times(2:2:end)) > (max_ITI * timeline_sample_rate + buffer_samples));
    if ~isempty(iti_violations)
        for iITI_violation = 1:length(iti_violations)
            this_ITI_violation = iti_violations(iITI_violation);
            if iITI_violation == 1 %first occurence
                stimOn_times = [stimOn_times(1:this_ITI_violation), NaN, stimOn_times(this_ITI_violation+1:end)];
            else
                stimOn_times = [stimOn_times(1:iti_violations(iITI_violation-1)+iITI_violation), ...
                    NaN, stimOn_times(iti_violations(iITI_violation-1)+1+iITI_violation:end)];
            end

        end
    end

    % else, use x first photodiode flips
    if length(signals_events.stimOnTimes) > length(stimOn_times) % if no violation or correction didn't work
        stimOn_times = stimOn_times(1:length(signals_events.stimOnTimes));
    end

end

% sanity checks
% - times between stim on times in signals
signals_photodiode_iti_diff = diff(signals_events.stimOnTimes(2:end)) - diff(stimOn_times) - 0.5';
if any(signals_photodiode_iti_diff > 0.1)
    warning('mismatching signals/photodiode stim ITIs')
end
% - ITI s makes sense 
timeline_sample_rate = 0.001;
buffer_samples = 200;
min_ITI = min(signals_events.stimITIsValues);
max_ITI = max(signals_events.stimITIsValues);
photodiode_iti = diff(stimOn_times);
if any(photodiode_iti < (min_ITI * timeline_sample_rate - buffer_samples)) ||...
        any(photodiode_iti > (min_ITI * timeline_sample_rate + buffer_samples))
    warning('unusual ITIs')
end
