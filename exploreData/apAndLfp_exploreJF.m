%% look at a chunk of raw data: ap and lfp 
%% recording info
animal = 'JF017';
day = '2020-10-22';
site = 2;
apPath = AP_cortexlab_filenameJF(animal,day,[],'ephys_ap',site,[]);
chanPath = AP_cortexlab_filenameJF(animal,day,[],'ephys',site,[]);
lfpPath = 
%% load ap data
chansToPlot = find(chanDistances<170);
    
    %get spike locations 
    pull_spikeT = -40:41;
    thisC=ephysData(iCount).spike_templates(iCluster);
    theseTimesCenter=ephysData(iCount).spike_times(ephysData(iCount).spike_templates==thisC);
    theseTimesCenter=theseTimesCenter(theseTimesCenter > timeChunkStart/ephysData(iCount).ephys_sample_rate);
    theseTimesCenter=theseTimesCenter(theseTimesCenter < timeChunkStop/ephysData(iCount).ephys_sample_rate);
    if ~isempty(theseTimesCenter)
        theseTimesFull = reshape(theseTimesCenter*ephysData(iCount).ephys_sample_rate+pull_spikeT, [size(theseTimesCenter,1)*size(pull_spikeT,2),1]);
    end
    %plot
%   cCount=cumsum(repmat(abs(max(max(memMapData(chansToPlot, timeChunkStart:timeChunkStop)))),size(chansToPlot,1),1),1);
    cCount=cumsum(repmat(1000,size(chansToPlot,1),1),1);
    

    t=timeChunkStart:timeChunkStop;
    LinePlotReducer(@plot, t, double(memMapData(chansToPlot, timeChunkStart:timeChunkStop))+double(cCount),'k');
    if ~isempty(theseTimesCenter)
        hold on;
        
        LinePlotReducer(@plot, theseTimesFull,double(memMapData(chansToPlot, theseTimesFull))+double(cCount),'r');
    end
    LinePlotExplorer(gcf);
    
%% load lfp data 