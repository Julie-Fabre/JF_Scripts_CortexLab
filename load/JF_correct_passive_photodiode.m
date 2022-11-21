% Check number of stim matches photodiode
% (once I saw some weird extra flip at the end of an
% experiment? so added case to only use first n flips)
timeline_sample_rate = 0.001;
buffer_samples = 200;
min_ITI = min(signals_events.stimITIsValues);
max_ITI = max(signals_events.stimITIsValues);

% correct for incorrect photodiode flip number
flip_problem = 0;
if length(signals_events.stimOnTimes) < length(stimOn_times) % too many photodiode flips
    warning([animal, ' ', day, ': different stim number signals and photodiode']);
    flip_problem = 1;
    

    % do any of them violate min ITI - if so, get rid of both
    % photodiode flips that violate this
    iti_violations = diff(photodiode_flip_times(2:2:end)) < (min_ITI * timeline_sample_rate - buffer_samples);
    if any(iti_violations)
        stimOn_times(logical([0; iti_violations])) = []; % remove culprit
        stimOn_times(logical(iti_violations)) = NaN; % NaN-out the time of closest culprit too- this time can't be trusted
    end
    if length(signals_events.stimOnTimes) < length(stimOn_times) % if no violation or correction didn't work
        % else, use x first photodiode flips
        stimOn_times = stimOn_times(1:length(signals_events.stimOnTimes));
    end
elseif length(signals_events.stimOnTimes) > length(stimOn_times) % some photodiode flips missing
    warning([animal, ' ', day, ': different stim number signals and photodiode']);
    flip_problem = 2;
    
    % do any of them violate min ITI / stim on time - if so, pad with NaNs at missing spots  
    % stim On times
    stimOn_average_samples = nanmedian([photodiode_flip(3:2:end) - photodiode_flip(2:2:end-1)]);
    on_photodiode = photodiode_flip(3:2:end);
    off_photodiode = photodiode_flip(2:2:end);
    photodiode_flip_added_stim_on = photodiode_flip(1:end);
    bad_on_count = 0; 
    for iFlip = 1:size(photodiode_flip_added_stim_on(3:2:end),1)
        if on_photodiode(iFlip) - off_photodiode(iFlip) > (stimOn_average_samples + buffer_samples)
            bad_on_count = bad_on_count + 1;
            photodiode_flip_added_stim_on = [photodiode_flip_added_stim_on(1:(iFlip)*2  ); NaN; ...
                photodiode_flip_added_stim_on((iFlip )*2  + 1: end)];
            on_photodiode = photodiode_flip_added_stim_on(3:2:end);
            off_photodiode = photodiode_flip_added_stim_on(2:2:end);
        end
    end
    stimOn_violations = find(isnan(photodiode_flip_added_stim_on)); % QQ make more robust:

    
    % iti violations 
    iti_violations = find([photodiode_flip_added_stim_on(4:2:end) - photodiode_flip_added_stim_on(3:2:end-1)] ...
        > (max_ITI ./ timeline_sample_rate + buffer_samples)); % QQ make more robust:

    photodiode_flip_added_stim_off = photodiode_flip_added_stim_on;
    if ~isempty(iti_violations)
        for iITI_violation = 1:length(iti_violations)
            this_ITI_violation = iti_violations(iITI_violation);

            if iITI_violation == 1 %first occurence
                photodiode_flip_added_stim_off = [photodiode_flip_added_stim_off(1:this_ITI_violation); NaN; NaN; photodiode_flip_added_stim_off(this_ITI_violation+1:end)];
            else
                photodiode_flip_added_stim_off = [photodiode_flip_added_stim_off(1:iti_violations(iITI_violation-1)+iITI_violation); ...
                    NaN; NaN; photodiode_flip_added_stim_off(iti_violations(iITI_violation-1)+1+iITI_violation:end)];
            end

        end
    end
    stimOn_times = photodiode_flip_added_stim_off(2:2:end);
    stimOn_times(~isnan(stimOn_times))= stimScreen_on_t(stimOn_times(~isnan(stimOn_times)));
    
    % attempt to get missing stimOn_times ? 
    

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
% - ITI s makes sense % QQ implement correction for this - happens sometimes
% "fast flips" of photodiode - especially present with short stimon times
% and short itis

photodiode_iti = diff(stimOn_times);
if any(photodiode_iti < (min_ITI * timeline_sample_rate - buffer_samples)) || ...
        any(photodiode_iti > (min_ITI * timeline_sample_rate + buffer_samples))
    warning('unusual ITIs - not correcting this yet')
end

if debug
    samples_to_plot = 305000;
    figure('Color', 'white');
    clf;
    title('photodiode')
    hold on;
    % raw trace
    plot(Timeline.rawDAQData(1:samples_to_plot, photodiode_idx))
    % original detected flips
    scatter(photodiode_flip(find(photodiode_flip <= samples_to_plot)), ones(size(find(photodiode_flip <= samples_to_plot), 1), 1), 8, 'filled')
    if flip_problem == 2
        % if too few flips, detected "missing" periods
        scatter(nanmean([photodiode_flip_added_stim_on(stimOn_violations-1), photodiode_flip_added_stim_on(stimOn_violations + 1)],2), ...
            ones(size(stimOn_violations,1),1)*2.5, filled')
        scatter(nanmean([photodiode_flip_added_stim_on(iti_violations*2+1), photodiode_flip_added_stim_on(iti_violations*2+2)],2), ...
            ones(size(iti_violations,1),1)*2)
        % if too few flips, estimated missing flips times

        legend({'raw trace', 'detected flips', 'missing stim on flips periods', 'missing stim off flips periods'})
   
    elseif flip_problem == 1
        % if too many flips, detected "surplus" periods
    end
    xlabel('time (in samples)')
     %plot(photodiode_trace_medfilt(1:samples_to_plot))


end
