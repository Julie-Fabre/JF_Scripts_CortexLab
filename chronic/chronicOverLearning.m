animalsAll = {'JF067'};
implantDays = {'2022-01-24'};
implantSide = 90;
plotDriftMap = 0; 
ex = [4,6];
startTr = 3;
for iAnimal = 1:size(animalsAll, 2)

     figure(4); 
     clf;
    animal = animalsAll{1, iAnimal}; %this animal
    implantDay = implantDays{1, iAnimal}; %this animal
    protocol1 = 'JF_choiceworldStimuli_wheel'; %protocol name contains this name
    protocol2 = 'JF_choiceworldStimuli_wheel_left_center';
    protocol3 = 'stage';
    protocol0 = 'JF_choiceworldStimuli';
    flexible_name = true; %protocol name can be slightly different
    clearvars experiments
    experiments = AP_find_experimentsJF(animal, protocol0, flexible_name);
    experiments3 = AP_find_experimentsJF(animal, protocol1, flexible_name);
    startTraining = experiments3(1).day;
    startTraining = split(startTraining, '-');
    startTraining = datetime(str2num(startTraining{1}), str2num(startTraining{2}), str2num(startTraining{3}));
    keepPassive = zeros(size(experiments, 1), 1);
%     for iPassiveDate = 1:min(size(experiments0, 1), 10) %look max 10 first days
%         expD = split(experiments0(iPassiveDate).day, '-');
%         keepPassive(iPassiveDate) = datetime(str2num(expD{1}), str2num(expD{2}), str2num(expD{3})) < startTraining;
%     end
%    experiments0 = experiments0(logical(keepPassive));
    %only keep experiments0 befoire training sdtarts
%     experiments2 = AP_find_experimentsJF(animal, protocol2, flexible_name);
%     experiments1 = AP_find_experimentsJF(animal, protocol3, flexible_name);
%     experiments(1:size(experiments0, 1), 1) = experiments0;
%     experiments(size(experiments0, 1)+1:size(experiments0, 1)+size(experiments1, 1), 1) = experiments1;
%     experiments(size(experiments0, 1)+size(experiments1, 1)+1:size(experiments0, 1)+size(experiments1, 1)+size(experiments2, 1), 1) = experiments2;
%     %experiments(size(experiments0, 1)+size(experiments1, 1)+size(experiments2, 1)+1:size(experiments0, 1)+size(experiments1, 1)+size(experiments2, 1)+size(experiments3, 1), 1) = experiments3;
    %experiments(end)=[];
    bhv = struct; %initialize structure
    keep_day = [];
    noGoDay = [];
    experiments(ex)=[];
    movingFrac = keep(iAnimal).movingFrac;
    movingFrac(ex - startTr + 1,:)=[];
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
            
            if plotDriftMap
                figure(1);
                subplot(1, size(experiments, 1), curr_day)
                [spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap([ephys_filename, '/site1/']);
                plotDriftmap(spikeTimes, amplitudes, spikeDepths);
                makeprettyLarge;
                xlim([0, max(spikeTimes)])
                ylim([0, 2880])

                xticks([0, max(spikeTimes)])
                xticklabels({'0', num2str(max(spikeTimes)/60)})
                if curr_day == 1
                    ylabel('Depth from tip (\mum)');
                    xlabel('time (min)')
                    makeprettyLarge;
                else
                    set(gca, 'ytick', [])
                    set(gca, 'yticklabel', [])
                    ylabel('')
                end
            end
            
            day = experiments(curr_day).day;
            experiment = experiments(curr_day).experiment(end); % last one = passive 
            site = 1;%1,1; 2,4; 3,7
            recording = []; 
            [ephysAPfile,aa] = AP_cortexlab_filenameJF(animal,day,experiment,'ephys_ap',site,recording);
            isSpikeGlx = contains(ephysAPfile, 'g0') | contains(ephysAPfile, 'g1')  | contains(ephysAPfile, 'g2')  | contains(ephysAPfile, 'g3') ;%spike glx (2.0 probes) or open ephys (3A probes)? 
            loadClusters=0;
            AP_load_experimentJF;
            
            figure(3)
            subplot(1, size(experiments, 1), curr_day)
            
            norm_spike_n = mat2gray(log10(accumarray(spike_templates, 1)+1));
            unit_dots = plot(norm_spike_n, template_depths, '.k', 'MarkerSize', 20);
            xlim([-0.1, 1]);
            ylim([-50, 2880 + 50]);
            if curr_day == 1
                ylabel('Depth (\mum)')
                xlabel('Normalized log rate')
            else
                set(gca, 'ytick', [])
                set(gca, 'yticklabel', [])
            end
            makeprettyLarge;
            
            figure(4);
            if curr_day >= startTr
                subplot(2, size(experiments, 1), curr_day)
                plot(keep(iAnimal).binBorders(1:end-1), movingFrac(curr_day - startTr +1 , :), 'b');
                box off;
                if curr_day == startTr
                    ylabel('frac. mov.')
                    %legend({'go1, contra','go1, center','other, contra','other, center'})
                else
                    set(gca,'ytick',[])
                    set(gca,'yticklabel',[])
                    set(gca, 'XColor','white')
                    set(gca, 'YColor','white')
                end
                line([0, 0], [0 1],'Color','k')
                makepretty;
                ylim([ 0 1])
                
            end
            stimType = nan(size(trial_conditions,1),1);
            stimType(trial_conditions(:,1) == 4 & trial_conditions(:,2) == -implantSide(iAnimal))=1;% go1 stim, contra
            stimType(trial_conditions(:,1) == 4 & trial_conditions(:,2) == 0)=2;% go1 stim, center
            stimType(ismember(trial_conditions(:,1), [1,2,3,5,6]) & trial_conditions(:,2) == -implantSide(iAnimal))=3;% other stims, contra
            stimType(ismember(trial_conditions(:,1), [1,2,3,5,6]) & trial_conditions(:,2) == 0)=4;% % other stims, center
      
            
            subplot(2, size(experiments, 1), curr_day+size(experiments,1))
            colorsO = [rgb('Purple');...
                    rgb('RoyalBlue');rgb('Turquoise'); rgb('Green'); ];
            alphaV=1;
            szV=2;

            [pltTime1, pltY1, pltTime2, pltY2,pltTime3, pltY3, pltCol3, pltTime4, pltY4, pltG4] =...
                plotExCellRasterPSTH(stimType, unique(stimType(~isnan(stimType))), ...
                colorsO, stimOn_times,spike_templates,...
                spike_times_timeline, unique(spike_templates),alphaV,szV,9,0, 0);
            pltY4= (pltY4 - nanmean(pltY4(:, pltTime4<0),2))./nanmean(pltY4(:, pltTime4<0),2);
            psth_lines = plot(pltTime4,pltY4); hold on;
            arrayfun(@(align_group) set(psth_lines(align_group), ...
                'XData',pltTime4,'YData',pltY4(align_group,:), ...
                'Color',colorsO(align_group,:)),1:size(pltY4,1));
            
            xlim([pltTime4(1), 0.25])
            %ylim([-300, 1000])
            line([0, 0], [-0.37, 1.6],'Color','k')
            box off;
            
            if curr_day == 1
                %ylabel('sp/s')
                legend({'go1, contra','go1, center','other, contra','other, center'})
            else
                set(gca,'ytick',[])
                set(gca,'yticklabel',[])
                set(gca, 'XColor','white')
                set(gca, 'YColor','white')
            end
            if curr_day == min(2, size(experiments,1))
                xlabel('time from stim onset')
            end
            
             makepretty;
             ylim([-0.37, 1.6])
             
%             pos = get(gca, 'Position');
%             spacing = 0.01:1/(size(experiments, 1)+1):0.99;
%             pos(1) = spacing(curr_day);
%             pos(3) = spacing(curr_day+1);
%             set(gca, 'Position', pos)
        else
            unitCount(curr_day) = NaN;
            deadChannels(curr_day) = NaN;
        end
    end
    figure(2);
    clf;
    plot(alldaysNum, unitCount, 'Color', 'k')
    hold on;

    grid minor;
    ylabel('# of units')
    makeprettyLarge;
    yyaxis right;
    plot(alldaysNum, 384-deadChannels, 'Color', 'b')
    makeprettyLarge;
    ax = gca;
    ax.YAxis(2).Color = 'b';
    xlabel('post-implant day #')
    ylabel('# of channels with no units')
    xlim([alldaysNum(1), alldaysNum(end)])
    
    
    

end
figure(6)
subplot(2, 2, 1)
plot(keep(iAnimal).binBorders(1:end-1), movingFrac(curr_day-1 - startTr +1 , :), 'b');
box off;
if curr_day == startTr
    ylabel('frac. mov.')
    %legend({'go1, contra','go1, center','other, contra','other, center'})
else
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    set(gca, 'XColor','white')
    set(gca, 'YColor','white')
end
makepretty;

subplot(2, 2, 2)
plot(keep(iAnimal).binBorders(1:end-1), movingFrac(curr_day - startTr +1 , :), 'b');
box off;
if curr_day == startTr
    ylabel('frac. mov.')
    %legend({'go1, contra','go1, center','other, contra','other, center'})
else
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    set(gca, 'XColor','white')
    set(gca, 'YColor','white')
end
makepretty;