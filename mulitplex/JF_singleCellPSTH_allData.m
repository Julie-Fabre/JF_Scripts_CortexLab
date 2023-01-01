function allData_singleCellPSTH = JF_singleCellPSTH_allData(allDataStruct,trialGroups, alignTo, raster_window, psth_bin_size)
unitCount = 0;
for iRecording = 1:size(allDataStruct,2)
    uniqueUnits = unique(allDataStruct(iRecording).spike_templates);
    for iUnit = 1:length(uniqueUnits)
        unitCount = unitCount + 1;
        if strcmp(trialGroups, 'stim_id')
            align_group = allDataStruct(iRecording).trial_conditions(:,1) + 10*(2+allDataStruct(iRecording).trial_conditions(:,2));
        end
        if strcmp(alignTo, 'stim')
            align_times = allDataStruct(iRecording).stimOn_times;
        end

        [curr_smoothed_psth, ~, ~, ~, ~] = JF_raster_PSTH(allDataStruct(iRecording).spike_templates,...
            allDataStruct(iRecording).spike_times, ...
            uniqueUnits(iUnit), raster_window, psth_bin_size, align_times, align_group, [],[], 0, 1);
        allData_singleCellPSTH.psth(unitCount,:,:) = curr_smoothed_psth;

    end
end
end