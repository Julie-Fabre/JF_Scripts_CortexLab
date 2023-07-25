

if gogogo
    task_data_here = task_data_gogogo;
else
    task_data_here = task_data; 
end
figure();
%clf;
cl_plottingSettings;

contra = 1;
if contra
    go1_contra = squeeze(task_data_here.psth{2}(:,:,2,:));
    go2_contra = squeeze(task_data_here.psth{2}(:,:,6,:));
    noGo_contra = squeeze(task_data_here.psth{2}(:,:,4,:));
else
    go1_contra = squeeze(task_data_here.psth{2}(:,:,3,:));
    go2_contra = squeeze(task_data_here.psth{2}(:,:,7,:));
    noGo_contra = squeeze(task_data_here.psth{2}(:,:,5,:));
end


% zscore_psth(1,:) = (task_data.psth{2}(:,1,1,:) - nanmean(task_data.psth(:,1:50),2)) ./ ...
%     nanstd(task_data.psth(:,1:50),[],2);
region_max = [1, 1, 1, 1, 1, 1, 1];
region_smooth = [5, 1, 1, 1, 1, 1, 1];
%region_clim_string = {'z-score (clim saturated)', 'z-score', 'z-score', 'z-score', 'z-score', 'z-score', 'z-score'};
for iRegion = 1:size(regions,2)
    

    
    % get all cells %1:301
    these_units = task_data_here.unit_area==iRegion &...
        task_data_here.pvalue_shuffled_005{2}(1:size(task_data_here.unit_area,1))' == 1;% & task_data.pvalue{2}' < 0.05 &...
        %task_data.pvalue_shuffled_005{2}' == 1;%& task_data.pvalue{2}' < 0.05 &...
        %(task_data.unitType' ==1 | task_data.unitType' ==2) ;

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
    test_image_go2 = squeeze([go2_contra(these_units, 2, :); go2_contra(these_units, 1, :)]);
    test_image_go2 = (test_image_go2 - nanmean(test_image_go2(:,1:50),2)) ./ ...
        nanstd(test_image_go2(:,1:50),[],2);
    smooth_filt = [region_smooth(iRegion),1]; % (units x frames)
    this_image_smooth_go2 = conv2(test_image_go2(cell_idx,:),ones(smooth_filt),'same')./ ...
    conv2(~isnan(test_image_go2(cell_idx,:)),ones(smooth_filt),'same');

    %no go
    test_image_noGo = squeeze([noGo_contra(these_units, 2, :); noGo_contra(these_units, 2, :)]);
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
    imagesc(task_data_here.t, [], this_image_smooth_go1(keep_these,:))
    caxis([-max(max(abs(this_image_smooth_go1(keep_these,:)))).*region_max(iRegion), max(max(abs(this_image_smooth_go1(keep_these,:)))).*region_max(iRegion)])
    c = colorbar;
    %c.Label.String = (region_clim_string{iRegion});
    colormap(brewermap([],'*BrBG'));
    xlabel('time from stim onset (s)')
    ylabel('unit # (sorted on 1/2 trials)')
    title('Go 1')
    makepretty;

    subplot(132)
    imagesc(task_data_here.t, [], this_image_smooth_go2(keep_these,:))
    caxis([-max(max(abs(this_image_smooth_go2(keep_these,:)))).*region_max(iRegion), max(max(abs(this_image_smooth_go2(keep_these,:)))).*region_max(iRegion)])
    c = colorbar;
    %c.Label.String = (region_clim_string{iRegion});
    colormap(brewermap([],'*BrBG'));
    xlabel('time from stim onset (s)')
    ylabel('unit #')
    title('Go 2')
    makepretty;

    subplot(133)
    imagesc(task_data_here.t, [], this_image_smooth_noGo(keep_these,:))
    caxis([-max(max(abs(this_image_smooth_noGo(keep_these,:)))).*region_max(iRegion), max(max(abs(this_image_smooth_noGo(keep_these,:)))).*region_max(iRegion)])
    c = colorbar;
    %c.Label.String = (region_clim_string{iRegion});
    colormap(brewermap([],'*BrBG'));
    xlabel('time from stim onset (s)')
    ylabel('unit #')
    title('No go')
    makepretty;



end
%plotshaded(0:0.001:0.5, [-acgstd + acgmean; acgstd + acgmean], 'g');
keep_these2 = find(keep_these);
%keep_these2 = keep_these2(80:end);
colorMtx = bc_colors(4);
if gogogo
    thisColor = [214, 37, 41]./256;
    num = 1.4;
else
    thisColor = [143, 103, 169]./256;
    num=1;
end
figure(2); 
title('Go1');
 hold on;
plot(task_data_here.t, nanmean(this_image_smooth_go1(keep_these2,:)), 'Color',thisColor); hold on;
plotshaded(task_data_here.t, [-nanstd(this_image_smooth_go1(keep_these2,:))./sqrt(sum(keep_these2)) + nanmean(this_image_smooth_go1(keep_these2,:));...
    nanstd(this_image_smooth_go1(keep_these2,:))./sqrt(sum(keep_these2)) + nanmean(this_image_smooth_go1(keep_these2,:))], thisColor);
legend('Go no go -trained', '', 'pre-learning', '', 'Go go go trained')
makepretty;
xlim([-0.5, 0.8])
ylim([-1, 6])

figure(3); 
title('Go2'); hold on;
plot(task_data_here.t, nanmean(this_image_smooth_go2(keep_these2,:))*num, 'Color',thisColor); hold on;
plotshaded(task_data_here.t, [-nanstd(this_image_smooth_go2(keep_these2,:))./sqrt(sum(keep_these2)) + nanmean(this_image_smooth_go2(keep_these2,:))*num;...
    nanstd(this_image_smooth_go2(keep_these2,:))./sqrt(sum(keep_these2)) + nanmean(this_image_smooth_go2(keep_these2,:))*num],  thisColor);
legend('Go no go -trained', '', 'pre-learning', '', 'Go go go trained')
makepretty;
xlim([-0.5, 0.8])
ylim([-1, 6])

figure(4); 
title('No Go'); hold on;
plot(task_data_here.t, nanmean(this_image_smooth_noGo(keep_these2,:))*num, 'Color',thisColor); hold on;
plotshaded(task_data_here.t, [-nanstd(this_image_smooth_noGo(keep_these2,:))./sqrt(sum(keep_these2)) + nanmean(this_image_smooth_noGo(keep_these2,:))*num;...
    nanstd(this_image_smooth_noGo(keep_these2,:))./sqrt(sum(keep_these2)) + nanmean(this_image_smooth_noGo(keep_these2,:))*num],  thisColor);
legend('Go no go -trained', '', 'pre-learning', '', 'Go go go trained')
makepretty;
xlim([-0.5, 0.8])
ylim([-1, 6])
