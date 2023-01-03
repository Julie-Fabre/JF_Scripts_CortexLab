function [instHit_rate, instCR_rate ] = JF_getBehavArousalMeasures(stimIDs, response, stimTimes)

correctGos = stimIDs(1:size(response,2)) ==1 & response==-1;
correctNoGos = stimIDs(1:size(response,2)) == 3 & response==0;
instHit_rate = movmean(correctGos, 40);
instCR_rate = movmean(correctNoGos, 40);

end