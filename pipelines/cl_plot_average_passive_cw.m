passive_data_cw_all = [passive_data_cw1, passive_data_cw2];

figure();
%clf;
cl_plottingSettings;

contra = 1;
if contra
    go1_contra = [squeeze(passive_data_cw_all(1).psth{5}(:,:,1,:));...
        squeeze(passive_data_cw_all(2).psth{4}(:,:,2,:))];
    go2_contra = [squeeze(passive_data_cw_all(1).psth{5}(:,:,5,:));...
        squeeze(passive_data_cw_all(2).psth{4}(:,:,6,:))];
    noGo_contra = [squeeze(passive_data_cw_all(1).psth{5}(:,:,3,:));...
        squeeze(passive_data_cw_all(2).psth{4}(:,:,4,:))];
else
      go1_contra = [squeeze(passive_data_cw_all(1).psth{5}(:,:,2,:));...
        squeeze(passive_data_cw_all(2).psth{4}(:,:,3,:))];
    go2_contra = [squeeze(passive_data_cw_all(1).psth{5}(:,:,6,:));...
        squeeze(passive_data_cw_all(2).psth{4}(:,:,7,:))];
    noGo_contra = [squeeze(passive_data_cw_all(1).psth{5}(:,:,4,:));...
        squeeze(passive_data_cw_all(2).psth{4}(:,:,5,:))];
end

unit_regions = [passive_data_cw_all(1).unit_area; passive_data_cw_all(2).unit_area];
pvalue_shuffled = [passive_data_cw_all(1).pvalue_shuffled_005{5}, passive_data_cw_all(2).pvalue_shuffled_005{4}];
pvalue = [passive_data_cw_all(1).pvalue{5}, passive_data_cw_all(2).pvalue{4}];

zscore_psth(:,1) = [(nanmean(squeeze(passive_data_cw_all(1).psth{5}(:,1,1,55:65)),2) - nanmean(squeeze(passive_data_cw_all(1).psth{5}(:,1,1,1:50)),2)) ./ ...
     nanstd(squeeze(passive_data_cw_all(1).psth{5}(:,1,1,1:50)),[],2);...
     (nanmean(squeeze(passive_data_cw_all(2).psth{4}(:,1,2,55:65)),2) - nanmean(squeeze(passive_data_cw_all(2).psth{4}(:,1,2,1:50)),2)) ./ ...
     nanstd(squeeze(passive_data_cw_all(2).psth{4}(:,1,2,1:50)),[],2)];

zscore_psth(:,2) = [(nanmean(squeeze(passive_data_cw_all(1).psth{5}(:,1,5,55:65)),2) - nanmean(squeeze(passive_data_cw_all(1).psth{5}(:,1,5,1:50)),2)) ./ ...
     nanstd(squeeze(passive_data_cw_all(1).psth{5}(:,1,5,1:50)),[],2);...
     (nanmean(squeeze(passive_data_cw_all(2).psth{4}(:,1,6,55:65)),2) - nanmean(squeeze(passive_data_cw_all(2).psth{4}(:,1,6,1:50)),2)) ./ ...
     nanstd(squeeze(passive_data_cw_all(2).psth{4}(:,1,2,1:50)),[],2)];

zscore_psth(:,3) = [(nanmean(squeeze(passive_data_cw_all(1).psth{5}(:,1,3,55:65)),2) - nanmean(squeeze(passive_data_cw_all(1).psth{5}(:,1,3,1:50)),2)) ./ ...
     nanstd(squeeze(passive_data_cw_all(1).psth{5}(:,1,3,1:50)),[],2);...
     (nanmean(squeeze(passive_data_cw_all(2).psth{4}(:,1,4,55:65)),2) - nanmean(squeeze(passive_data_cw_all(2).psth{4}(:,1,4,1:50)),2)) ./ ...
     nanstd(squeeze(passive_data_cw_all(2).psth{4}(:,1,4,1:50)),[],2)];
region_max = [1, 1, 1, 1, 1, 1, 1];

region_smooth = [5, 1, 1, 1, 1, 1, 1];
%region_clim_string = {'z-score (clim saturated)', 'z-score', 'z-score', 'z-score', 'z-score', 'z-score', 'z-score'};
for iRegion = 1%:size(regions,2)
    

    
    % get all cells
    these_units = unit_regions(1:size(task_data_gogogo.psth{2},1))==iRegion& any(abs(zscore_psth)>2, 2);% & passive_data_cw_all.pvalue{4}' < 0.05 &...
        %passive_data_cw_all.pvalue_shuffled_005{4}' == 1;%& passive_data_cw_all.pvalue{4}' < 0.05 &...
        %(passive_data_cw_all.unitType' ==1 | passive_data_cw_all.unitType' ==2) ;

    train_image = squeeze(go1_contra(these_units, 1, :));
    train_image = (train_image - nanmean(train_image(:,1:50),2)) ./ ...
        nanstd(train_image(:,1:50),[],2);
    [~, cell_idx] = sort(nanmean(train_image(:,55:70),2));
    
    %go 1
    test_image_go1 = squeeze(go1_contra(these_units, 2, :));
    test_image_go1 = (test_image_go1 - nanmean(test_image_go1(:,1:50),2)) ./ ...
        nanstd(test_image_go1(:,1:50),[],2);
    smooth_filt = [region_smooth(iRegion),1]; % (units x frames)
    this_image_smooth_go1 = conv2(test_image_go1(cell_idx,:),ones(smooth_filt),'same')./ ...
    conv2(~isnan(test_image_go1(cell_idx,:)),ones(smooth_filt),'same');

    %go 2
    test_image_go2 = squeeze([go2_contra(these_units, 2, :)]);
    test_image_go2 = (test_image_go2 - nanmean(test_image_go2(:,1:50),2)) ./ ...
        nanstd(test_image_go2(:,1:50),[],2);
    smooth_filt = [region_smooth(iRegion),1]; % (units x frames)
    this_image_smooth_go2 = conv2(test_image_go2(cell_idx,:),ones(smooth_filt),'same')./ ...
    conv2(~isnan(test_image_go2(cell_idx,:)),ones(smooth_filt),'same');

    %no go
    test_image_noGo = squeeze([noGo_contra(these_units, 2, :)]);
    test_image_noGo = (test_image_noGo - nanmean(test_image_noGo(:,1:50),2)) ./ ...
        nanstd(test_image_noGo(:,1:50),[],2);
    smooth_filt = [region_smooth(iRegion),1]; % (units x frames)
    this_image_smooth_noGo = conv2(test_image_noGo(cell_idx,:),ones(smooth_filt),'same')./ ...
    conv2(~isnan(test_image_noGo(cell_idx,:)),ones(smooth_filt),'same');
    
    % remove mostly NaN rows (to be replaced by bombcell output when that's
    % finished running)
    keep_these = sum(isnan(this_image_smooth_noGo),2)<100 & sum(isnan(this_image_smooth_go1),2)<100 &...
        sum(isnan(this_image_smooth_go2),2)<100 ;

    % plot PSTH 
    
    subplot(131)
    imagesc(passive_data_cw_all(1).t, [], this_image_smooth_go1(keep_these,:))
    caxis([-max(max(abs(this_image_smooth_go1(keep_these,:)))).*region_max(iRegion), max(max(abs(this_image_smooth_go1(keep_these,:)))).*region_max(iRegion)])
    c = colorbar;
    %c.Label.String = (region_clim_string{iRegion});
    colormap(brewermap([],'*BrBG'));
    xlabel('time from stim onset (s)')
    ylabel('unit # (sorted on 1/2 trials)')
    title('Go 1')
    makepretty;

    subplot(132)
    imagesc(passive_data_cw_all(1).t, [], this_image_smooth_go2(keep_these,:))
    caxis([-max(max(abs(this_image_smooth_go2(keep_these,:)))).*region_max(iRegion), max(max(abs(this_image_smooth_go2(keep_these,:)))).*region_max(iRegion)])
    c = colorbar;
    %c.Label.String = (region_clim_string{iRegion});
    colormap(brewermap([],'*BrBG'));
    xlabel('time from stim onset (s)')
    ylabel('unit #')
    title('Go 2')
    makepretty;

    subplot(133)
    imagesc(passive_data_cw_all(1).t, [], this_image_smooth_noGo(keep_these,:))
    caxis([-max(max(abs(this_image_smooth_noGo(keep_these,:)))).*region_max(iRegion), max(max(abs(this_image_smooth_noGo(keep_these,:)))).*region_max(iRegion)])
    c = colorbar;
    %c.Label.String = (region_clim_string{iRegion});
    colormap(brewermap([],'*BrBG'));
    xlabel('time from stim onset (s)')
    ylabel('unit #')
    title('No go')
    makepretty;



end

%figure(); 
%scatter 
passiveColor = [0, 0, 0]; %[214, 37, 41]./256;
figure(2); 
plot(passive_data_cw_all(1).t, nanmean(this_image_smooth_go1(keep_these,:)), 'Color', passiveColor); 
plotshaded(task_data.t, [-nanstd(this_image_smooth_go1(keep_these,:))./sqrt(sum(keep_these)) + nanmean(this_image_smooth_go1(keep_these,:));...
    nanstd(this_image_smooth_go1(keep_these,:))./sqrt(sum(keep_these)) + nanmean(this_image_smooth_go1(keep_these,:))],  passiveColor);

xlabel('time from stim onset (s)')
ylabel('zscore')
legend('pre-learning', 'post-learning')
ylim([-0.4, 3.5])
xlim([-0.5, 0.8])
makepretty;

figure(3); 
plot(passive_data_cw_all(1).t, nanmean(this_image_smooth_go1(keep_these,:)), 'Color', passiveColor); 
plotshaded(task_data.t, [-nanstd(this_image_smooth_go2(keep_these,:))./sqrt(sum(keep_these)) + nanmean(this_image_smooth_go2(keep_these,:));...
    nanstd(this_image_smooth_go2(keep_these,:))./sqrt(sum(keep_these)) + nanmean(this_image_smooth_go2(keep_these,:))],  passiveColor);
 xlabel('time from stim onset (s)')
ylabel('zscore')
legend('pre-learning', 'post-learning')
xlim([-0.5, 0.8])
makepretty;
% 
figure(4); 
plot(passive_data_cw_all(1).t, nanmean(this_image_smooth_noGo(keep_these,:)), 'Color', passiveColor); 
plotshaded(task_data.t, [-nanstd(this_image_smooth_noGo(keep_these,:))./sqrt(sum(keep_these)) + nanmean(this_image_smooth_noGo(keep_these,:));...
    nanstd(this_image_smooth_noGo(keep_these,:))./sqrt(sum(keep_these)) + nanmean(this_image_smooth_noGo(keep_these,:))],  passiveColor);
 xlabel('time from stim onset (s)')
ylabel('zscore')
legend('pre-learning', 'post-learning')
xlim([-0.5, 0.8])
makepretty;

