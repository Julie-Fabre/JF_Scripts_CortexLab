
%% SUAbugsBgone
% with some good old for loops
%clear all;
%load('E:\analysis\wf_ephys_choiceworld\paper\data\trial_activity_naive.mat')
if ~exist('removeMissingSpikesUnits', 'var')
    removeMissingSpikesUnits = 0;
end

for iAnimal = 1:size(trial_data_all.goodUnits, 1)
    for iRecording = 1:size(trial_data_all.goodUnits{iAnimal}, 1)

        %% 1. if flag remove spikes missing, add these to the goodUnits
        if removeMissingSpikesUnits
            trial_data_all.goodUnits{iAnimal}{iRecording} = trial_data_all.goodUnits{iAnimal}{iRecording} ...
                & nanmin(trial_data_all.percent_missing_ndtr{iAnimal}{iRecording}(2:end, :)) < removeMissingSpikesUnits;
        end

        %% 2. re-classify to: add noise guys to uins, remove large post spike suppression FSI and UINs, and short waveform guys from TANs. 
        if ~isempty(trial_data_all.acg{iAnimal}{iRecording})
         for iCell = 1:size(trial_data_all.acg{iAnimal}{iRecording}, 1)
            pss_allcat2temp = find(trial_data_all.acg{iAnimal}{iRecording}(iCell, 500:1000) >= nanmean(trial_data_all.acg{iAnimal}{iRecording}(iCell, 600:900)));
            pss_allcat2(iCell) = pss_allcat2temp(1);
         end
        allDepths = zeros(size(trial_data_all.allGroups{iAnimal}{iRecording},2),1);
        allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [1,4,7,10,13,16]))=1; 
        allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [2,5,8,11,14,17]))=2;
        allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [3,6,9,12,15,18]))=3;
        largePssShorties = find(trial_data_all.templateDuration{iAnimal}{iRecording} < 400 & pss_allcat2' > 40); 
        
        fsi = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 < 40 & trial_data_all.prop_long_isi{iAnimal}{iRecording} <= 0.1;
        uin = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 < 40 & trial_data_all.prop_long_isi{iAnimal}{iRecording} > 0.1;
        tan = trial_data_all.templateDuration{iAnimal}{iRecording} > 400 & pss_allcat2 >= 40;
        msn = trial_data_all.templateDuration{iAnimal}{iRecording} > 400 & pss_allcat2 < 40;
        shortDurLongPss = trial_data_all.templateDuration{iAnimal}{iRecording} < 400 & pss_allcat2 >= 40;
        
        trial_data_all.allGroups{iAnimal}{iRecording}(msn) = allDepths(msn);
        trial_data_all.allGroups{iAnimal}{iRecording}(fsi) = allDepths(fsi) + 3;
        trial_data_all.allGroups{iAnimal}{iRecording}(tan) = allDepths(tan) + 6;
        trial_data_all.allGroups{iAnimal}{iRecording}(uin) = allDepths(uin) + 9;
        trial_data_all.allGroups{iAnimal}{iRecording}(shortDurLongPss) = allDepths(shortDurLongPss) + 12;
        trial_data_all.allGroups{iAnimal}{iRecording}(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [16,17,18])) =...
            allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [16,17,18]))+12;%ones un-caught that should be in the shortLong category
        clearvars pss_allcat2
        end
    end
end



