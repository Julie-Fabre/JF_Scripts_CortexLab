function [instHit_rate, instCR_rate, instGo_rate ] = JF_getBehavArousalMeasures(stimIDs, response)

Gos = response==-1 | response==1;
correctGos = stimIDs(1:length(response)) ==1 & response==-1;
correctNoGos = stimIDs(1:length(response)) == 3 & response==0;
instHit_rate = correctGos;
instCR_rate = correctNoGos;
instGo_rate = Gos;

end