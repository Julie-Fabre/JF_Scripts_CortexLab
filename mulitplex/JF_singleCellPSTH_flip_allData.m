function JF_singleCellPSTH_flip_allData(allDataStruct,trialGroups, alignTo, raster_window, psth_bin_size)
%F_singleCellPSTH_flip_allData(allDataStruct,'stim_id', 'stim', raster_window, psth_bin_size);
%nitCount=uni
for iRecording = 1:size(allDataStruct,2)
    uniqueUnits = unique(allDataStruct(iRecording).spike_templates);
    for iUnit = 1:length(uniqueUnits)
        %unitCount = unitCount + 1;
        iUnit = iUnit + 1;
        if strcmp(trialGroups, 'stim_id')
            align_group = allDataStruct(iRecording).trial_conditions(:,1) + 10*(2+allDataStruct(iRecording).trial_conditions(:,2));
        end
        if strcmp(alignTo, 'stim')
            align_times = allDataStruct(iRecording).stimOn_times;
        end
     align_group = isnan(allDataStruct(iRecording).stim_to_move(allDataStruct(iRecording).trial_conditions(:,1)==1,1));
     align_times = allDataStruct(iRecording).stimOn_times(allDataStruct(iRecording).trial_conditions(:,1)==1,1);
     sort_by =  allDataStruct(iRecording).stim_to_move(allDataStruct(iRecording).trial_conditions(:,1)==1,1);
        [curr_smoothed_psth, ~, raster_x, raster_y, raster_color, alignedVector] = JF_raster_PSTH(allDataStruct(iRecording).spike_templates,...
            allDataStruct(iRecording).spike_times, ...
            uniqueUnits(iUnit), raster_window, psth_bin_size, align_times, align_group,sort_by,[], 1, 1, 5);

             align_group = isnan(allDataStruct(iRecording).stim_to_move(allDataStruct(iRecording).trial_conditions(:,1)==3,1));
     align_times = allDataStruct(iRecording).stimOn_times(allDataStruct(iRecording).trial_conditions(:,1)==3,1);
     sort_by =  allDataStruct(iRecording).stim_to_move(allDataStruct(iRecording).trial_conditions(:,1)==3,1);
        [curr_smoothed_psth, ~, raster_x, raster_y, raster_color, alignedVector] = JF_raster_PSTH(allDataStruct(iRecording).spike_templates,...
            allDataStruct(iRecording).spike_times, ...
            uniqueUnits(iUnit), raster_window, psth_bin_size, align_times, align_group,sort_by,[], 1, 1, 6);

    end
end
end