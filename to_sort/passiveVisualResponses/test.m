 for iBinX = 1:size(Xedges,2)
            for iBinY = 1:size(Yedges,2)
                theseNeurons = binX==iBinX & binY==iBinY;
                theseNeuronsInd = find(theseNeurons);
                binnedArrayTot = [];
                if ~isempty(theseNeuronsInd)
                    [h,hh] = histc(theseNeuronsInd, theseLocationsInfo);
                    uniqueRecs = unique(hh)+1;
                    for iUniqueRec = 1:size(uniqueRecs,1) %get psth per rec 
                        theseTheseNeurons = theseNeuronsInd(hh == uniqueRecs(iUniqueRec)); 
                        if ~isempty(ephysData(uniqueRecs(iUniqueRec)).spike_times_timeline) && ~isempty(ephysData(uniqueRecs(iUniqueRec)).stimOn_times)
                            theseTheseNeuronsTemplate = ismember(ephysData(uniqueRecs(iUniqueRec)).spike_templates,...
                                theseTheseNeurons-theseLocationsInfo(uniqueRecs(iUniqueRec)));
                        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = ...
                            psthAndBA(ephysData(uniqueRecs(iUniqueRec)).spike_times_timeline(theseTheseNeuronsTemplate),...
                            ephysData(uniqueRecs(iUniqueRec)).stimOn_times, [-0.2, 0.3], 0.01);
                        binnedArrayTot = [binnedArrayTot; binnedArray]; 
                            
                        end
                    end
                    
                end 
                if ~isempty(binnedArrayTot)
                    binnedArrayTotBinned(iBinX,iBinY,:) = nanmean(binnedArrayTot,1); %average 
                else
                    binnedArrayTotBinned(iBinX,iBinY,:) = nan(50,1); %average 
                end
            end
 end
 