

function [timeBins, posBins1_full, posBins2_full, posBins3_full, allP, normVals] = psthByPos3D(spikeTimes, spikePos1, spikePos2, spikePos3, posBinSize, timeBinSize, eventTimes, win, bslWin, varargin)
% function [timeBins, depthBins, allP] = psthByDepth(spikeTimes, ...
%   spikeDepths, depthBinSize, timeBinSize, eventTimes, win, bslWin[, bslEvents])
%
% Computes PSTH split by the depths of spikes provided
% If pass a "bslWin" will normalize to the bin counts in that period
% - if so normalized, units of allP is stdev relative to baseline mean
% - if not normalized, units of allP is spikes/sec
%
% timeBins is 1xnTimeBins
% depthBins is 1xnDepthBins
% allP is nDepthBins x nTimeBins
%


posBins1 = min(spikePos1):posBinSize:max(spikePos1);
nD1 = length(posBins1) - 1;

posBins2 = min(spikePos2):posBinSize:max(spikePos2);
nD2 = length(posBins2) - 1;

posBins3 = min(spikePos3):posBinSize:max(spikePos3);
nD3 = length(posBins3) - 1;

if ~isempty(varargin)
    bslEventTimes = varargin{1};
else
    bslEventTimes = eventTimes;
end

combs = allcomb((1:nD1,1:nD2,1:nD3);
posBins1_full = posBins1(combs(1,:));
posBins2_full = posBins2(combs(2,:));
posBins3_full = posBins3(combs(3,:));
if ~isempty(bslWin)
    normVals = zeros(size(combs,2), 2);
else
    normVals = [];
end

for iComb = 1:size(combs, 2)
    d1 = combs(1, iComb);
    d2 = combs(2, iComb);
    d3 = combs(3, iComb);
    theseSp = spikePos1 > posBins1(d1) & spikePos1 <= posBins1(d1+1) & spikePos2 > posBins2(d2) & spikePos2 <= posBins2(d2+1) & ...
        spikePos3 > posBins3(d3) & spikePos3 <= posBins3(d3+1);

    if ~isempty(bslWin)
        [psth, ~, ~, ~, ~, ~] = psthAndBA(spikeTimes(theseSp), bslEventTimes, bslWin, timeBinSize);
        normMn = mean(psth);
        normStd = std(psth);
    end

    [psth, timeBins, ~, ~, ~, ~] = psthAndBA(spikeTimes(theseSp), eventTimes, win, timeBinSize);
    if iComb == 1
            allP = zeros(iComb, length(psth));
    end
        
    if ~isempty(bslWin) && normStd > 0
        allP(iComb, :) = (psth - normMn) ./ normStd;
        normVals(iComb, :) = [normMn, normStd];
    else
        allP(iComb, :) = psth;
    end


end