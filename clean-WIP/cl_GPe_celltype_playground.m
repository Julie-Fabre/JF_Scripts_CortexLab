zscore_psth = (passive_data.psth - nanmean(passive_data.psth(:,1:50),2)) ./ ...
    nanstd(passive_data.psth(:,1:50),[],2);

iRegion = 2
    %2 =GPe
    % remove mostly NaN rows (to be replaced by bombcell output when that's
    % finished running)
    keep_these = sum(isnan(zscore_psth),2)<100;
    
    % get all cells
    these_units = passive_data.unit_area==iRegion  & ...
    (passive_data.unitType' ==1 | passive_data.unitType' ==2);
    
    figure();
   % subplot(131)
    stem3(passive_data.wvDur(these_units),  passive_data.pss(these_units), passive_data.fr(these_units));
    hold on;
    xlabel('waveform duration (us)')
    ylabel('post spike suppression (ms)')
    zlabel('firing rate')
    makepretty;
    xlim([150, 900])

    figure();
    stem3(passive_data.wvDur(these_units),  passive_data.propISI(these_units), passive_data.fr(these_units));
    hold on;
    xlabel('waveform duration (us)')
    ylabel('prop ISI > 2s')
    zlabel('firing rate')
    makepretty;
    xlim([150, 900])

    % different visual responses? 


