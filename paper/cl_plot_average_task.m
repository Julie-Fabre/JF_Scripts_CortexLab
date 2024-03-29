
function cl_plot_average_task(task_data, idx, contra, center, keepUnits, keepVis)


cl_plottingSettings;


if contra
    % 1,2,3 = all the +90s
    % 13,14,15 %go1 [4]
    % 19,20,21 % nogo [6]
    % 37,38,39 %go2 [12]
    % 40,41,42 %go likes [13]
    % 16,17,18; 31,32,33; 34,35,36 %nogo likes [5, 10, 11]
    if size(task_data.psth{idx},3)==6
        go1_contra = squeeze(task_data.psth{idx}(:, :, 4, :));
        go2_contra = squeeze(task_data.psth{idx}(:, :, 6, :));
        noGo_contra = squeeze(task_data.psth{idx}(:, :, 5, :));
    else
        go1_contra = squeeze(task_data.psth{idx}(:, :, 10, :));
        go2_contra = squeeze(task_data.psth{idx}(:, :, 34, :));
        noGo_contra = squeeze(task_data.psth{idx}(:, :, 16, :));
        noGoLike_contra = squeeze(task_data.psth{idx}(:, :, 19, :));
        goLike_contra = squeeze(task_data.psth{idx}(:, :, 37, :));
        noGoLike2_contra = squeeze(task_data.psth{idx}(:, :, 28, :));
        noGoLike3_contra = squeeze(task_data.psth{idx}(:, :, 31, :));
    end
elseif center
    if size(task_data.psth{idx},3)==6
        go1_contra = squeeze(task_data.psth{idx}(:, :, 1, :));
        go2_contra = squeeze(task_data.psth{idx}(:, :, 3, :));
        noGo_contra = squeeze(task_data.psth{idx}(:, :, 2, :));
    else
    go1_contra = squeeze(task_data.psth{idx}(:, :, 11, :));
    go2_contra = squeeze(task_data.psth{idx}(:, :, 35, :));
    noGo_contra = squeeze(task_data.psth{idx}(:, :, 17, :));
    noGoLike_contra = squeeze(task_data.psth{idx}(:, :, 20, :));
    goLike_contra = squeeze(task_data.psth{idx}(:, :, 38, :));
    noGoLike2_contra = squeeze(task_data.psth{idx}(:, :, 29, :));
    noGoLike3_contra = squeeze(task_data.psth{idx}(:, :, 32, :));
    end
else
    go1_contra = squeeze(task_data.psth{idx}(:, :, 12, :));
    go2_contra = squeeze(task_data.psth{idx}(:, :, 36, :));
    noGo_contra = squeeze(task_data.psth{idx}(:, :, 18, :));
    noGoLike_contra = squeeze(task_data.psth{idx}(:, :, 21, :));
    goLike_contra = squeeze(task_data.psth{idx}(:, :, 39, :));
    noGoLike2_contra = squeeze(task_data.psth{idx}(:, :, 30, :));
    noGoLike3_contra = squeeze(task_data.psth{idx}(:, :, 33, :));
end


% zscore_psth(1,:) = (task_data.psth{idx}(:,1,1,:) - nanmean(task_data.psth(:,val_t_1),2)) ./ ...
%     nanstd(task_data.psth(:,val_t_1),[],2);
region_max = [1, 1, 1, 1, 1, 1, 1];
region_smooth = [1, 1, 1, 1, 1, 1, 1];
region_lims = [5, 2, 2, 1, 1, 1];
plot_regions = [1, 2, 3];



%region_clim_string = {'z-score (clim saturated)', 'z-score', 'z-score', 'z-score', 'z-score', 'z-score', 'z-score'};
figure();
passive = 1; %to disable tons of plots
for thisRegion = 1:size(plot_regions, 2)
try
    iRegion = plot_regions(thisRegion);

    % get all cells %1:301
    if keepVis
     these_units = task_data.unit_area == iRegion & ...
          ismember(task_data.unitType', keepUnits) &...
    (task_data.pvalue_shuffled_005{1,idx}(1:size(task_data.unit_area, 1))' == 1);
    else
    these_units = task_data.unit_area == iRegion & ...
        ismember(task_data.unitType', keepUnits);


    end% &...
    %task_data.wvDur'<=400&task_data.pss'<10000 );% |...
    %any(task_data.pvalue{idx}(:,:)>0.975,2) | any(task_data.pvalue{idx}(:,:)<0.025,2) );% & ...
    %task_data.pvalue_shuffled_005{idx}(1:size(task_data.unit_area, 1))' == 1; % & task_data.pvalue{idx}' < 0.05 &...
    %task_data.pvalue_shuffled_005{idx}' == 1;%& task_data.pvalue{idx}' < 0.05 &...
    %(task_data.unitType' ==1 | task_data.unitType' ==2) ;

    train_image = squeeze(go1_contra(these_units, 1, :));
    if size(train_image, 2) == 800
        val_t_1 = [1:150];
        val_t_2 = [250:400];
    else
        val_t_1 = [1:20];
        val_t_2 = [55:75];

    end
    % if sum(these_units)<10
    %     region_smooth(iRegion) = 1;
    % elseif sum(these_units)<20
    %     region_smooth(iRegion) = 2;
    % elseif sum(these_units)<30
    %     region_smooth(iRegion) = 3;
    % elseif sum(these_units)<40
    %     region_smooth(iRegion) = 4;
    % else
    %     region_smooth(iRegion) = 5;
    % end

    train_image = (train_image - nanmean(train_image(:, val_t_1), 2)) ./ ...
        nanstd(train_image(:, val_t_1), [], 2);
    [~, cell_idx] = sort(nanmean(train_image(:, val_t_2), 2));

    %go 1
    test_image_go1 = squeeze(go1_contra(these_units, 2, :));
    test_image_go1 = (test_image_go1 - nanmean(test_image_go1(:, val_t_1), 2)) ./ ...
        nanstd(test_image_go1(:, val_t_1), [], 2);


    smooth_filt = [region_smooth(iRegion), 10]; % (units x frames)
    this_image_smooth_go1 = conv2(test_image_go1(cell_idx, :), ones(smooth_filt), 'same') ./ ...
        conv2(~isnan(test_image_go1(cell_idx, :)), ones(smooth_filt), 'same');

    %go 2
    test_image_go2 = squeeze(go2_contra(these_units, 2, :));
    test_image_go2 = (test_image_go2 - nanmean(test_image_go2(:, val_t_1), 2)) ./ ...
        nanstd(test_image_go2(:, val_t_1), [], 2);
    smooth_filt = [region_smooth(iRegion), 10]; % (units x frames)
    this_image_smooth_go2 = conv2(test_image_go2(cell_idx, :), ones(smooth_filt), 'same') ./ ...
        conv2(~isnan(test_image_go2(cell_idx, :)), ones(smooth_filt), 'same');

    %no go
    test_image_noGo = squeeze(noGo_contra(these_units, 2, :));
    test_image_noGo = (test_image_noGo - nanmean(test_image_noGo(:, val_t_1), 2)) ./ ...
        nanstd(test_image_noGo(:, val_t_1), [], 2);
    smooth_filt = [region_smooth(iRegion), 10]; % (units x frames)
    this_image_smooth_noGo = conv2(test_image_noGo(cell_idx, :), ones(smooth_filt), 'same') ./ ...
        conv2(~isnan(test_image_noGo(cell_idx, :)), ones(smooth_filt), 'same');

    if ~passive
        %no go like
        test_image_noGoLike = squeeze([noGoLike_contra(these_units, 2, :); noGoLike_contra(these_units, 2, :)]);
        test_image_noGoLike = (test_image_noGoLike - nanmean(test_image_noGoLike(:, val_t_1), 2)) ./ ...
            nanstd(test_image_noGoLike(:, val_t_1), [], 2);
        smooth_filt = [region_smooth(iRegion), 10]; % (units x frames)
        this_image_smooth_noGoLike = conv2(test_image_noGoLike(cell_idx, :), ones(smooth_filt), 'same') ./ ...
            conv2(~isnan(test_image_noGoLike(cell_idx, :)), ones(smooth_filt), 'same');

        %no go like 2
        test_image_noGoLike2 = squeeze([noGoLike2_contra(these_units, 2, :); noGoLike2_contra(these_units, 2, :)]);
        test_image_noGoLike2 = (test_image_noGoLike2 - nanmean(test_image_noGoLike2(:, val_t_1), 2)) ./ ...
            nanstd(test_image_noGoLike2(:, val_t_1), [], 2);
        smooth_filt = [region_smooth(iRegion), 10]; % (units x frames)
        this_image_smooth_noGoLike2 = conv2(test_image_noGoLike2(cell_idx, :), ones(smooth_filt), 'same') ./ ...
            conv2(~isnan(test_image_noGoLike2(cell_idx, :)), ones(smooth_filt), 'same');

        %no go like 3
        test_image_noGoLike3 = squeeze([noGoLike3_contra(these_units, 2, :); noGoLike3_contra(these_units, 2, :)]);
        test_image_noGoLike3 = (test_image_noGoLike3 - nanmean(test_image_noGoLike3(:, val_t_1), 2)) ./ ...
            nanstd(test_image_noGoLike3(:, val_t_1), [], 2);
        smooth_filt = [region_smooth(iRegion), 10]; % (units x frames)
        this_image_smooth_noGoLike3 = conv2(test_image_noGoLike3(cell_idx, :), ones(smooth_filt), 'same') ./ ...
            conv2(~isnan(test_image_noGoLike3(cell_idx, :)), ones(smooth_filt), 'same');

        %go like
        test_image_goLike = squeeze(goLike_contra(these_units, 2, :));
        test_image_goLike = (test_image_goLike - nanmean(test_image_goLike(:, val_t_1), 2)) ./ ...
            nanstd(test_image_goLike(:, val_t_1), [], 2);
        smooth_filt = [region_smooth(iRegion), 10]; % (units x frames)
        this_image_smooth_goLike = conv2(test_image_goLike(cell_idx, :), ones(smooth_filt), 'same') ./ ...
            conv2(~isnan(test_image_goLike(cell_idx, :)), ones(smooth_filt), 'same');
    end

    % remove mostly NaN rows (to be replaced by bombcell output when that's
    % finished running)
    keep_these = sum(isnan(this_image_smooth_noGo), 2) < 100 & sum(isnan(this_image_smooth_go1), 2) < 100 & ...
        sum(isnan(this_image_smooth_go2), 2) < 100;

    % plot PSTH

    if ~passive
        subplot(size(plot_regions, 2), 7, (7 * (thisRegion - 1))+1)
    else
        subplot(size(plot_regions, 2), 3, (3 * (thisRegion - 1))+1)
    end
    %subplot(131)
    imagesc(task_data.t, [], this_image_smooth_go1(keep_these, :))
    caxis([-max(max(abs(this_image_smooth_go1(keep_these, :)))) .* region_max(iRegion), max(max(abs(this_image_smooth_go1(keep_these, :)))) .* region_max(iRegion)])
    c = colorbar;
    %c.Label.String = (region_clim_string{iRegion});
    colormap(brewermap([], '*BrBG'));
    if iRegion == 1
        xlabel('time from stim onset (s)')
        ylabel('unit # (sorted on 1/2 trials)')
    end
    title('Go 1')
    makepretty;
    clim([-region_lims(iRegion), region_lims(iRegion)])

    if ~passive
        subplot(size(plot_regions, 2), 7, (7 * (thisRegion - 1))+2)
    else
        subplot(size(plot_regions, 2), 3, (3 * (thisRegion - 1))+2)
    end

    imagesc(task_data.t, [], this_image_smooth_go2(keep_these, :))
    caxis([-max(max(abs(this_image_smooth_go2(keep_these, :)))) .* region_max(iRegion), max(max(abs(this_image_smooth_go2(keep_these, :)))) .* region_max(iRegion)])
    c = colorbar;
    %c.Label.String = (region_clim_string{iRegion});
    colormap(brewermap([], '*BrBG'));

    title('Go 2')
    makepretty;
    clim([-region_lims(iRegion), region_lims(iRegion)])

    if ~passive
        subplot(size(plot_regions, 2), 7, (7 * (thisRegion - 1))+3)
    else
        subplot(size(plot_regions, 2), 3, (3 * (thisRegion - 1))+3)
    end
    imagesc(task_data.t, [], this_image_smooth_noGo(keep_these, :))
    caxis([-max(max(abs(this_image_smooth_noGo(keep_these, :)))) .* region_max(iRegion), max(max(abs(this_image_smooth_noGo(keep_these, :)))) .* region_max(iRegion)])
    c = colorbar;
    %c.Label.String = (region_clim_string{iRegion});
    colormap(brewermap([], '*BrBG'));

    title('Stim 3')
    makepretty;
    clim([-region_lims(iRegion), region_lims(iRegion)])

    if ~passive
        subplot(size(plot_regions, 2), 7, (7 * (thisRegion - 1))+4)

        imagesc(task_data.t, [], this_image_smooth_noGoLike(keep_these, :))
        caxis([-max(max(abs(this_image_smooth_noGoLike(keep_these, :)))) .* region_max(iRegion), max(max(abs(this_image_smooth_noGoLike(keep_these, :)))) .* region_max(iRegion)])
        c = colorbar;
        %c.Label.String = (region_clim_string{iRegion});
        colormap(brewermap([], '*BrBG'));

        title('No go -like')
        makepretty;
        %clim([-region_lims(iRegion),region_lims(iRegion)] )

        subplot(size(plot_regions, 2), 7, (7 * (thisRegion - 1))+5)

        imagesc(task_data.t, [], this_image_smooth_noGoLike2(keep_these, :))
        caxis([-max(max(abs(this_image_smooth_noGoLike2(keep_these, :)))) .* region_max(iRegion), max(max(abs(this_image_smooth_noGoLike2(keep_these, :)))) .* region_max(iRegion)])
        c = colorbar;
        %c.Label.String = (region_clim_string{iRegion});
        colormap(brewermap([], '*BrBG'));

        title('No go -like')
        makepretty;
        %clim([-region_lims(iRegion),region_lims(iRegion)] )

        subplot(size(plot_regions, 2), 7, (7 * (thisRegion - 1))+6)

        imagesc(task_data.t, [], this_image_smooth_noGoLike3(keep_these, :))
        caxis([-max(max(abs(this_image_smooth_noGoLike3(keep_these, :)))) .* region_max(iRegion), max(max(abs(this_image_smooth_noGoLike3(keep_these, :)))) .* region_max(iRegion)])
        c = colorbar;
        %c.Label.String = (region_clim_string{iRegion});
        colormap(brewermap([], '*BrBG'));

        title('No go -like')
        makepretty;
        %clim([-region_lims(iRegion),region_lims(iRegion)] )

        subplot(size(plot_regions, 2), 7, (7 * (thisRegion - 1))+7)

        imagesc(task_data.t, [], this_image_smooth_goLike(keep_these, :))
        caxis([-max(max(abs(this_image_smooth_goLike(keep_these, :)))) .* region_max(iRegion), max(max(abs(this_image_smooth_goLike(keep_these, :)))) .* region_max(iRegion)])
        c = colorbar;
        %c.Label.String = (region_clim_string{iRegion});
        colormap(brewermap([], '*BrBG'));

        title('Go -like')
        makepretty;
        %clim([-region_lims(iRegion),region_lims(iRegion)] )

    end
catch
end
end
more = 0;
if more
    %plotshaded(0:0.001:0.5, [-acgstd + acgmean; acgstd + acgmean], 'g');
    keep_these2 = find(keep_these);
    %keep_these2 = keep_these2(80:end);
    colorMtx = bc_colors(4);
    if gogogo
        thisColor = [214, 37, 41] ./ 256;
        num = 1.4;
    else
        thisColor = [143, 103, 169] ./ 256;
        num = 1;
    end
    figure(2);
    title('Go1');
    hold on;
    plot(task_data.t, nanmean(this_image_smooth_go1(keep_these2, :)), 'Color', thisColor);
    hold on;
    plotshaded(task_data.t_det, [-nanstd(this_image_smooth_go1(keep_these2, :)) ./ sqrt(sum(keep_these2)) + nanmean(this_image_smooth_go1(keep_these2, :)); ...
        nanstd(this_image_smooth_go1(keep_these2, :)) ./ sqrt(sum(keep_these2)) + nanmean(this_image_smooth_go1(keep_these2, :))], thisColor);
    legend('Go no go -trained', '', 'pre-learning', '', 'Go go go trained')
    makepretty;
    xlim([-0.5, 0.8])
    ylim([-1, 6])

    figure(3);
    title('Go2');
    hold on;
    plot(task_data.t, nanmean(this_image_smooth_go2(keep_these2, :))*num, 'Color', thisColor);
    hold on;
    plotshaded(task_data.t_det, [-nanstd(this_image_smooth_go2(keep_these2, :)) ./ sqrt(sum(keep_these2)) + nanmean(this_image_smooth_go2(keep_these2, :)) * num; ...
        nanstd(this_image_smooth_go2(keep_these2, :)) ./ sqrt(sum(keep_these2)) + nanmean(this_image_smooth_go2(keep_these2, :)) * num], thisColor);
    legend('Go no go -trained', '', 'pre-learning', '', 'Go go go trained')
    makepretty;
    xlim([-0.5, 0.8])
    ylim([-1, 6])

    figure(4);
    title('No Go');
    hold on;
    plot(task_data.t, nanmean(this_image_smooth_noGo(keep_these2, :))*num, 'Color', thisColor);
    hold on;
    plotshaded(task_data.t_det, [-nanstd(this_image_smooth_noGo(keep_these2, :)) ./ sqrt(sum(keep_these2)) + nanmean(this_image_smooth_noGo(keep_these2, :)) * num; ...
        nanstd(this_image_smooth_noGo(keep_these2, :)) ./ sqrt(sum(keep_these2)) + nanmean(this_image_smooth_noGo(keep_these2, :)) * num], thisColor);
    legend('Go no go -trained', '', 'pre-learning', '', 'Go go go trained')
    makepretty;
    xlim([-0.5, 0.8])
    ylim([-1, 6])
end
end