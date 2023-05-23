%% cl_population PSTH

zscore_psth = (passive_data.psth - nanmean(passive_data.psth(:,1:50),2)) ./ ...
    nanstd(passive_data.psth(:,1:50),[],2);

for iRegion = 1:size(regions,2)
    
    % remove mostly NaN rows (to be replaced by bombcell output when that's
    % finished running)
    keep_these = sum(isnan(zscore_psth),2)<100;
    
    % get all cells
    these_units = passive_data.unit_area==iRegion;


    % sort cells by activity
    this_image = zscore_psth(these_units&keep_these, :);
    [~, cell_idx] = sort(nanmean(this_image(:,55:70),2));

    % small smoothing -otherwise matlab doesn't display properly 
    smooth_filt = [10,1]; % (units x frames)
    this_image_smooth = conv2(this_image(cell_idx,:),ones(smooth_filt),'same')./ ...
    conv2(~isnan(this_image(cell_idx,:)),ones(smooth_filt),'same');

    % plot PSTH 
    figure();
    imagesc(passive_data.t, [], this_image)
    caxis([-max(max(abs(this_image_smooth))).*0.2, max(max(abs(this_image_smooth))).*0.2])
    c = colorbar;
    c.Label.String = ({'z-score (clim saturated)'});
    colormap(brewermap([],'*BrBG'));
    xlabel('time from stim onset (s)')
    ylabel('unit #')
    makepretty;
    

end