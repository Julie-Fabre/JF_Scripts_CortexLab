%% cl_population PSTH
close all;

zscore_psth = (passive_data.psth - nanmean(passive_data.psth(:,1:50),2)) ./ ...
    nanstd(passive_data.psth(:,1:50),[],2);
region_max = [0.2, 1, 1, 1, 1, 1, 1];
region_smooth = [10, 1, 1, 1, 1, 1, 1];
region_clim_string = {'z-score (clim saturated)', 'z-score', 'z-score', 'z-score', 'z-score', 'z-score', 'z-score'};
for iRegion = 1:size(regions,2)
    
    % remove mostly NaN rows (to be replaced by bombcell output when that's
    % finished running)
    keep_these = sum(isnan(zscore_psth),2)<100;
    
    % get all cells
    these_units = passive_data.unit_area==iRegion;

    if sum(these_units) > 0
    % sort cells by activity
    this_image = zscore_psth(these_units&keep_these, :);
    [~, cell_idx] = sort(nanmean(this_image(:,55:70),2));

    % small smoothing -otherwise matlab doesn't display properly 
    smooth_filt = [region_smooth(iRegion),1]; % (units x frames)
    this_image_smooth = conv2(this_image(cell_idx,:),ones(smooth_filt),'same')./ ...
    conv2(~isnan(this_image(cell_idx,:)),ones(smooth_filt),'same');

    % plot PSTH 
    figure();
    imagesc(passive_data.t, [], this_image_smooth)
    caxis([-max(max(abs(this_image_smooth))).*region_max(iRegion), max(max(abs(this_image_smooth))).*region_max(iRegion)])
    c = colorbar;
    c.Label.String = (region_clim_string{iRegion});
    colormap(brewermap([],'*BrBG'));
    xlabel('time from stim onset (s)')
    ylabel('unit #')
    title(regions{iRegion})
    makepretty;
    end
    

end



%% example cells 
CP = [3427, 3488];
cp_units = find(passive_data.unit_area==iRegion & sum(isnan(zscore_psth),2)<100);

passive_data.animal_day_site_shank(cp_units(3427),:);