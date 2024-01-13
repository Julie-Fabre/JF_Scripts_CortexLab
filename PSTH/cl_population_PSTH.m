%% cl_population PSTH
%close all;
index = 1;
average_across_images = squeeze(nanmean(passive_data.psth{index}(:,3,:,:), 3));


zscore_psth = (average_across_images - nanmean(average_across_images(:,1:200),2)) ./ ...
    (nanstd(average_across_images(:,1:200),[],2)+0.001);
region_max = [1, 1, 1, 1, 1, 1, 1];
region_smooth = [5, 1, 1, 1, 1, 1, 1];
region_clim_string = {'z-score', 'z-score', 'z-score', 'z-score', 'z-score', 'z-score', 'z-score'};
for iRegion = [1,2,3]
    
    % remove mostly NaN rows (to be replaced by bombcell output when that's
    % finished running)
    keep_these = sum(isnan(zscore_psth),2) < 100 & sum(zscore_psth==0,2) < 100 &...
        abs(nanmean(zscore_psth(:,250:450),2)) > 0.25;
    
    % get all cells
    these_units = passive_data.unit_area==iRegion  &...
        (passive_data.unitType' ==1 | passive_data.unitType' ==2);%& passive_data.pvalue_shuffled_005{index}' == 1;

    if sum(these_units) > 0
        % sort cells by activity
        this_image = zscore_psth(these_units & keep_these, :);
        [~, cell_idx] = sort(nanmean(this_image(:,250:450),2));
    
        % small smoothing -otherwise matlab doesn't display properly 
        smooth_filt = [region_smooth(iRegion),10]; % (units x frames)
        %this_image_smooth = this_image(cell_idx,:);
        this_image_smooth = conv2(this_image(cell_idx,:),ones(smooth_filt),'same')./ ...
            conv2(~isnan(this_image(cell_idx,:)),ones(smooth_filt),'same');
    
        % plot PSTH 
        figure();
        imagesc(passive_data.t, [], this_image_smooth(3:end,:))
       % caxis([-max(max(abs(this_image_smooth))).*region_max(iRegion), max(max(abs(this_image_smooth))).*region_max(iRegion)])
        c = colorbar;
        c.Label.String = (region_clim_string{iRegion});
        colormap(brewermap([],'*BrBG'));
        xlabel('time from stim onset (s)')
        ylabel('unit #')
        title(regions{iRegion})
        makepretty;
        caxis([-2.5, 2.5])
     end
    

end



%% example cells 
% CP = [3427, 3488];
% cp_units = find(passive_data.unit_area==iRegion & sum(isnan(zscore_psth),2)<100);
% 
% passive_data.animal_day_site_shank(cp_units(3427),:);