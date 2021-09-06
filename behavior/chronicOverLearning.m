animalsAll = {'JF051'};
implantDays = {'2021-08-27'};
for iAnimal = 1:size(animalsAll, 2)


    animal = animalsAll{1, iAnimal}; %this animal
    implantDay = implantDays{1, iAnimal}; %this animal
    protocol1 = 'location'; %protocol name contains this name
    protocol2 = 'goNoGo';
    protocol3 = 'Hack';
    protocol0 = 'JF_choiceworldStimuli';
    flexible_name = true; %protocol name can be slightly different
    clearvars experiments
    experiments0 = AP_find_experimentsJF(animal, protocol0, flexible_name);
    experiments1 = AP_find_experimentsJF(animal, protocol1, flexible_name);
    startTraining = experiments1(1).day;
    startTraining = split(startTraining, '-');
    startTraining = datetime(str2num(startTraining{1}), str2num(startTraining{2}), str2num(startTraining{3}));
    keepPassive = zeros(size(experiments0, 1), 1);
    for iPassiveDate = 1:min(size(experiments0, 1), 10) %look max 10 first days
        expD = split(experiments0(iPassiveDate).day, '-');
        keepPassive(iPassiveDate) = datetime(str2num(expD{1}), str2num(expD{2}), str2num(expD{3})) < startTraining;
    end
    experiments0 = experiments0(logical(keepPassive));
    %only keep experiments0 befoire training sdtarts
    experiments2 = AP_find_experimentsJF(animal, protocol2, flexible_name);
    experiments3 = AP_find_experimentsJF(animal, protocol3, flexible_name);
    experiments(1:size(experiments0, 1), 1) = experiments0;
    experiments(size(experiments0, 1)+1:size(experiments0, 1)+size(experiments1, 1), 1) = experiments1;
    experiments(size(experiments0, 1)+size(experiments1, 1)+1:size(experiments0, 1)+size(experiments1, 1)+size(experiments2, 1), 1) = experiments2;
    experiments(size(experiments0, 1)+size(experiments1, 1)+size(experiments2, 1)+1:size(experiments0, 1)+size(experiments1, 1)+size(experiments2, 1)+size(experiments3, 1), 1) = experiments3;
    %experiments(end)=[];
    bhv = struct; %initialize structure
    keep_day = [];
    noGoDay = [];
    alldays = {experiments.day};
    for curr_day = 1:size(experiments, 1)
        alldaysNum(curr_day) = hours((datetime(alldays{curr_day}, 'InputFormat', 'yyyy-MM-dd') - datetime(implantDay, 'InputFormat', 'yyyy-MM-dd'))/24);
        % if kilosorted, load the spikes, spike_templates, amplitudes,
        % spike_depths
        thisDay = experiments(curr_day).day;
        thisExperiment = experiments(curr_day).experiment(end);
        [ephys_filename, ~] = AP_cortexlab_filenameJF(animal, thisDay, thisExperiment, 'ephys');
        if ~isempty(dir([ephys_filename, '/site1/spike_templates.npy'])) %is has been already kilosorted
            amplitudes = readNPY([ephys_filename, '/site1/amplitudes.npy']);
            spike_templates = readNPY([ephys_filename, '/site1/spike_templates.npy']);

            template_depths = readNPY([ephys_filename, '/site1/templates.npy']);
            channel_positions = readNPY([ephys_filename, '/site1/channel_positions.npy']);
            channel_map = readNPY([ephys_filename, '/site1/channel_map.npy']);

            unitCount(curr_day) = length(unique(spike_templates));
            deadChannels(curr_day) = length(unique(channel_map));
            figure(1);
            subplot(1, size(experiments, 1), curr_day)
            [spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap([ephys_filename, '/site1/']);
            plotDriftmap(spikeTimes, amplitudes, spikeDepths);
            makeprettyLarge;
            xlim([0, max(spikeTimes)])
            if curr_day == 1
                xticks([0, max(spikeTimes)])
                xticklabels({'0', num2str(max(spikeTimes)/60)})
                ylabel('Depth from tip (\mum)');
                xlabel('time (min)')
                makeprettyLarge;
            end

        else
            unitCount(curr_day) = NaN;
            deadChannels(curr_day) = NaN;
        end
    end
    figure(2);
    plot(alldaysNum, unitCount, 'Color', 'k')
    hold on;
    ylabel('# of units')
    makeprettyLarge;
    yyaxis right;
    plot(alldaysNum, 384-deadChannels, 'Color', 'b')
    makeprettyLarge;
    ax = gca;
    ax.YAxis(2).Color = 'b';
    xlabel('post-implant day #')
    ylabel('# of dead channels')
    xlim([alldaysNum(1), alldaysNum(end)])


end
