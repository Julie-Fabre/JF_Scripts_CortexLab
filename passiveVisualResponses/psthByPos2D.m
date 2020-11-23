

function [timeBins, posBins1, posBins2, allP, normVals] = psthByPos2D(spikeTimes, spikePos1, spikePos2, posBinSize, timeBinSize, eventTimes, win, bslWin, varargin)
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
if ~isempty(varargin)
    bslEventTimes = varargin{1};
else
    bslEventTimes = eventTimes;
end

if ~isempty(bslWin)
    normVals = zeros(nD, 2);
else
    normVals = [];
end

for d1 = 1:nD1
    for d2 = 1:nD2
        theseSp = spikePos1 > posBins1(d1) & spikePos1 <= posBins1(d1+1) & spikePos2 > posBins2(d2) & spikePos2 <= posBins2(d2+1);

        if ~isempty(bslWin)
            [psth, ~, ~, ~, ~, ~] = psthAndBA(spikeTimes(theseSp), bslEventTimes, bslWin, timeBinSize);
            normMn = mean(psth);
            normStd = std(psth);
        end

        [psth, timeBins, ~, ~, ~, ~] = psthAndBA(spikeTimes(theseSp), eventTimes, win, timeBinSize);

        if d == 1
            allP = zeros(nD, length(psth));
        end
        if ~isempty(bslWin) && normStd > 0
            allP(d, :) = (psth - normMn) ./ normStd;
            normVals(d, :) = [normMn, normStd];
        else
            allP(d, :) = psth;
        end
    end
end