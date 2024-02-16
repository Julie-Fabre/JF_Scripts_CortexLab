
%% cl_population PSTH
%close all;
passive_data = load('/home/julie/Dropbox/MATLAB/naive_data1.mat');


method = 1;
regions = {'CP', 'GPe', 'SNr'};

index = 1;
average_across_images = squeeze(nanmean(passive_data.psth{index}(:, 3, :, :), 3));

% if size(for_baseline_per_neuron_per_condition, 2) == 39 || size(for_baseline_per_neuron_per_condition, 2) == 66
%                 cond_inds = [10, 34; 16, 34; 16, 10];
%             elseif size(for_baseline_per_neuron_per_condition, 2) == 26
%                 cond_inds = [7, 23; 11, 23; 11, 7];
%             elseif size(for_baseline_per_neuron_per_condition, 2) == 4
%                 cond_inds = [2, 4; 3, 4; 3, 2];
%             else
%                 disp('wtf')
%                 continue;
%             end
zscore_psth = (average_across_images - nanmean(average_across_images(:, 1:200), 2)) ./ ...
    (nanstd(average_across_images(:, 1:200), [], 2) + 0.001);
region_max = [1, 1, 1, 1, 1, 1, 1];
region_smooth = [3, 1, 1, 1, 1, 1, 1];
region_clim_string = {'z-score', 'z-score', 'z-score', 'z-score', 'z-score', 'z-score', 'z-score'};
figure();
for iRegion = 1:3

    % remove mostly NaN rows (to be replaced by bombcell output when that's
    % finished running)

    
    keep_these = sum(isnan(zscore_psth), 2) < 100 & sum(zscore_psth == 0, 2) < 100 &  abs(nanmean(zscore_psth(:, 220:250), 2)) < 0.25;
    
        
    if method == 0
         these_units = passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2);
    elseif method == 1%method 1: high zscore
         these_units = passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2) & abs(nanmean(zscore_psth(:, 260:300), 2)) > 0.25;
    elseif method ==  2
    %method 2: mean shuffle
      these_units = passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2) &  passive_data.pvalue_shuffled_005{index}' == 1;
    elseif method ==  3
    % method 3: peak shuffle 
      these_units = passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2) &...
        (passive_data.pvalue_shuffled_0025_peak{index}' == 1 | passive_data.pvalue_shuffled_0025_trough{index}' == 1);
    end
    
    frac = [];
    

    if sum(these_units) > 0
        % sort cells by activity
        this_image = zscore_psth(these_units & keep_these, :);
        [~, cell_idx] = sort(nanmean(this_image(:, 250:450), 2));

        % small smoothing -otherwise matlab doesn't display properly
        smooth_filt = [region_smooth(iRegion), 10]; % (units x frames)
        %this_image_smooth = this_image(cell_idx,:);
        this_image_smooth = conv2(this_image(cell_idx, :), ones(smooth_filt), 'same') ./ ...
            conv2(~isnan(this_image(cell_idx, :)), ones(smooth_filt), 'same');

        % plot PSTH
        subplot(1,3,iRegion)
        imagesc(passive_data.t, [], this_image_smooth(1:end, :))
        % caxis([-max(max(abs(this_image_smooth))).*region_max(iRegion), max(max(abs(this_image_smooth))).*region_max(iRegion)])
        c = colorbar;
        c.Label.String = (region_clim_string{iRegion});
        colormap(brewermap([], '*BrBG'));
        xlabel('time from stim onset (s)')
        ylabel('unit #')
        title(regions{iRegion})
        makepretty;
        caxis([-2.5, 2.5])
    end


end

%% striatum: seperate by region
index = 1;
average_across_images = squeeze(nanmean(passive_data.psth{index}(:, 3, :, :), 3));


zscore_psth = (average_across_images - nanmean(average_across_images(:, 1:200), 2)) ./ ...
    (nanstd(average_across_images(:, 1:200), [], 2) + 0.001);
region_max = [1, 1, 1, 1, 1, 1, 1];
region_smooth = [3, 1, 1, 1, 1, 1, 1];
region_clim_string = {'z-score', 'z-score', 'z-score', 'z-score', 'z-score', 'z-score', 'z-score'};
for iRegion = 1%:3

    figure();

    % get all cells
    
        
    if method == 0
         these_units = passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2);
    elseif method == 1%method 1: high zscore
         these_units = passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2) & abs(nanmean(zscore_psth(:, 260:300), 2)) > 0.25;
    elseif method ==  2
    %method 2: mean shuffle
      these_units = passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2) &  passive_data.pvalue_shuffled_005{index}' == 1;
    elseif method ==  3
    % method 3: peak shuffle 
      these_units = passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2) &...
        (passive_data.pvalue_shuffled_0025_peak{index}' == 1 | passive_data.pvalue_shuffled_0025_trough{index}' == 1);
    end

    %these_units = passive_data.unit_area == iRegion & ...
     %   (passive_data.unitType' == 1 | passive_data.unitType' == 2); %& passive_data.pvalue_shuffled_005{index}' == 1;

    [regionClassification, unitClassification] = cl_subsection_region(passive_data.unit_coords, ...
        passive_data.unit_area, passive_data.pss, passive_data.templateDuration, passive_data.propLongISI);
    unique_subReg = unique(regionClassification(passive_data.unit_area ==iRegion));
    unique_subReg(unique_subReg == "") =[];


    for iSubReg = 1:size(unique_subReg,1)
    % remove mostly NaN rows (to be replaced by bombcell output when that's
    % finished running)
    %keep_these = sum(isnan(zscore_psth), 2) < 100 & sum(zscore_psth == 0, 2) < 100 & ...
     %   regionClassification == unique_subReg(iSubReg);%abs(nanmean(zscore_psth(:, 250:450), 2)) > 0.25 & 

    keep_these = sum(isnan(zscore_psth), 2) < 100 & sum(zscore_psth == 0, 2) < 100 &  abs(nanmean(zscore_psth(:, 220:250), 2)) < 0.25 & ...
        regionClassification == unique_subReg(iSubReg);
    
   
        if sum(these_units) > 0
            % sort cells by activity
            this_image = zscore_psth(these_units & keep_these, :);
            [~, cell_idx] = sort(nanmean(this_image(:, 250:450), 2));

            % small smoothing -otherwise matlab doesn't display properly
            smooth_filt = [region_smooth(iRegion), 10]; % (units x frames)
            %this_image_smooth = this_image(cell_idx,:);
            this_image_smooth = conv2(this_image(cell_idx, :), ones(smooth_filt), 'same') ./ ...
                conv2(~isnan(this_image(cell_idx, :)), ones(smooth_filt), 'same');

            % plot PSTH
            subplot(1, size(unique_subReg,1), iSubReg)
            imagesc(passive_data.t, [], this_image_smooth(1:end, :))
            % caxis([-max(max(abs(this_image_smooth))).*region_max(iRegion), max(max(abs(this_image_smooth))).*region_max(iRegion)])
            c = colorbar;
            c.Label.String = (region_clim_string{iRegion});
            colormap(brewermap([], '*BrBG'));
            xlabel('time from stim onset (s)')
            ylabel('unit #')
            title(unique_subReg{iSubReg})
            makepretty;
            caxis([-2.5, 2.5])
        end
    end

end


%% seperate by cell type and region striatum 

index = 1;
average_across_images = squeeze(nanmean(passive_data.psth{index}(:, 3, :, :), 3));


zscore_psth = (average_across_images - nanmean(average_across_images(:, 1:200), 2)) ./ ...
    (nanstd(average_across_images(:, 1:200), [], 2) + 0.001);
region_max = [1, 1, 1, 1, 1, 1, 1];
region_smooth = [3, 1, 1, 1, 1, 1, 1];
region_clim_string = {'z-score', 'z-score', 'z-score', 'z-score', 'z-score', 'z-score', 'z-score'};
for iRegion = 1%:3

    figure();

   if method == 0
         these_units = passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2);
    elseif method == 1%method 1: high zscore
         these_units = passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2) & abs(nanmean(zscore_psth(:, 260:300), 2)) > 0.25;
    elseif method ==  2
    %method 2: mean shuffle
      these_units = passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2) &  passive_data.pvalue_shuffled_005{index}' == 1;
    elseif method ==  3
    % method 3: peak shuffle 
      these_units = passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2) &...
        (passive_data.pvalue_shuffled_0025_peak{index}' == 1 | passive_data.pvalue_shuffled_0025_trough{index}' == 1);
    end


    [regionClassification, unitClassification] = cl_subsection_region(passive_data.unit_coords, ...
        passive_data.unit_area, passive_data.pss, passive_data.templateDuration, passive_data.propLongISI);
    unique_subReg = unique(regionClassification(passive_data.unit_area ==iRegion));
    unique_subReg(unique_subReg == "") =[];

    unique_unit = unique(unitClassification(passive_data.unit_area ==iRegion));
    unique_unit(unique_unit == "") =[];


    for iSubReg = 1:size(unique_subReg,1)
        for iUnitType = 1:size(unique_unit,1)
    % remove mostly NaN rows (to be replaced by bombcell output when that's
    % finished running)
    keep_these = sum(isnan(zscore_psth), 2) < 100 & sum(zscore_psth == 0, 2) < 100 & ...
        regionClassification == unique_subReg(iSubReg) & unitClassification == unique_unit(iUnitType) &  abs(nanmean(zscore_psth(:, 220:250), 2)) < 0.25;
   
        if sum(these_units) > 0
            % sort cells by activity
            this_image = zscore_psth(these_units & keep_these, :);
            [~, cell_idx] = sort(nanmean(this_image(:, 250:450), 2));

            % small smoothing -otherwise matlab doesn't display properly
            smooth_filt = [region_smooth(iRegion), 10]; % (units x frames)
            %this_image_smooth = this_image(cell_idx,:);
            this_image_smooth = conv2(this_image(cell_idx, :), ones(smooth_filt), 'same') ./ ...
                conv2(~isnan(this_image(cell_idx, :)), ones(smooth_filt), 'same');

            % plot PSTH
            subplot(size(unique_unit,1), size(unique_subReg,1), (iSubReg-1)*size(unique_unit,1)+iUnitType)
            imagesc(passive_data.t, [], this_image_smooth(1:end, :))
            % caxis([-max(max(abs(this_image_smooth))).*region_max(iRegion), max(max(abs(this_image_smooth))).*region_max(iRegion)])
            c = colorbar;
            c.Label.String = (region_clim_string{iRegion});
            colormap(brewermap([], '*BrBG'));
            xlabel('time from stim onset (s)')
            ylabel('unit #')
            title(unique_subReg{iSubReg})
            makepretty;
            caxis([-1, 1])
        end
        end
    end

end

%% seperate all by "cell type" 
method = 0 ; 
currF = figure();
index = 1;
average_across_images = squeeze(nanmean(passive_data.psth{index}(:, 3, :, :), 3));

regions = {'CP', 'GPe', 'SNr'};
zscore_psth = (average_across_images - nanmean(average_across_images(:, 1:200), 2)) ./ ...
    (nanstd(average_across_images(:, 1:200), [], 2) + 0.001);
region_max = [1, 1, 1, 1, 1, 1, 1];
region_smooth = [2, 2, 2, 1, 1, 1, 1];
region_clim_string = {'z-score', 'z-score', 'z-score', 'z-score', 'z-score', 'z-score', 'z-score'};
for iRegion = 1:3

    

    % get all cells
   if method == 0
         these_units = passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2);
    elseif method == 1%method 1: high zscore
         these_units = passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2) & abs(nanmean(zscore_psth(:, 260:300), 2)) > 0.25;
    elseif method ==  2
    %method 2: mean shuffle
      these_units = passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2) &  passive_data.pvalue_shuffled_005{index}' == 1;
    elseif method ==  3
    % method 3: peak shuffle 
      these_units = passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2) &...
        (passive_data.pvalue_shuffled_0025_peak{index}' == 1 | passive_data.pvalue_shuffled_0025_trough{index}' == 1);
    end



    [regionClassification, unitClassification] = cl_subsection_region(passive_data.unit_coords, ...
        passive_data.unit_area, passive_data.pss, passive_data.templateDuration, passive_data.propLongISI);
    unique_subReg = unique(regionClassification(passive_data.unit_area ==iRegion));
    unique_subReg(unique_subReg == "") =[];

    unique_unit = unique(unitClassification(passive_data.unit_area ==iRegion));
    unique_unit(unique_unit == "") =[];

      keep_these = sum(isnan(zscore_psth), 2) < 100 & sum(zscore_psth == 0, 2) < 100 & ...
          abs(nanmean(zscore_psth(:, 220:250), 2)) < 0.25;

    figure(); 
subplot(1,3,1)
scatter(passive_data.templateDuration(these_units&keep_these), passive_data.pss(these_units&keep_these),3, 'filled' )
xlabel('peak to trough duration (us)')
ylabel('post spike suppression (ms)')
subplot(1,3,2)
scatter(passive_data.templateDuration(these_units&keep_these), passive_data.propLongISI(these_units&keep_these),3, 'filled' )
xlabel('peak to trough duration (us)')
ylabel('ISI > 2s (frac)')
subplot(1,3,3)
scatter(passive_data.pss(these_units&keep_these), passive_data.propLongISI(these_units&keep_these),3, 'filled' )
ylabel('ISI > 2s (frac)')
xlabel('post spike suppression (ms)')
prettify_plot;


    %for iSubReg = 1:size(unique_subReg,1)
        for iUnitType = 1:size(unique_unit,1)
    % remove mostly NaN rows (to be replaced by bombcell output when that's
    % finished running)
    keep_these = sum(isnan(zscore_psth), 2) < 100 & sum(zscore_psth == 0, 2) < 100 & ...
         unitClassification == unique_unit(iUnitType)  &  abs(nanmean(zscore_psth(:, 220:250), 2)) < 0.25;%abs(nanmean(zscore_psth(:, 250:450), 2)) > 0.25 & 

   
        if sum(these_units) > 0
            % sort cells by activity
            this_image = zscore_psth(these_units & keep_these, :);
            [~, cell_idx] = sort(nanmean(this_image(:, 250:450), 2));

            % small smoothing -otherwise matlab doesn't display properly
            smooth_filt = [region_smooth(iRegion), 10]; % (units x frames)
            %this_image_smooth = this_image(cell_idx,:);
            this_image_smooth = conv2(this_image(cell_idx, :), ones(smooth_filt), 'same') ./ ...
                conv2(~isnan(this_image(cell_idx, :)), ones(smooth_filt), 'same');
figure(currF)
            % plot PSTH
            subplot(3, size(unique_unit,1), (iRegion-1)*size(unique_unit,1)+iUnitType)
            imagesc(passive_data.t, [], this_image_smooth(1:end, :))
            % caxis([-max(max(abs(this_image_smooth))).*region_max(iRegion), max(max(abs(this_image_smooth))).*region_max(iRegion)])
            c = colorbar;
            c.Label.String = (region_clim_string{iRegion});
            colormap(brewermap([], '*BrBG'));
            xlabel('time from stim onset (s)')
            ylabel('unit #')
            title([regions{iRegion}, ' ', unique_unit(iUnitType)])
            makepretty;
            caxis([-2.5,2.5])
        end
        end
    %end

end

%% percentage cells 
method = 1;
regions = {'CP', 'GPe', 'SNr'};
cl_plottingSettings;
index = 1;

figure();
for iRegion = 1:3

    % remove mostly NaN rows (to be replaced by bombcell output when that's
    % finished running)

    
    keep_these = sum(isnan(zscore_psth), 2) < 100 & sum(zscore_psth == 0, 2) < 100 &  abs(nanmean(zscore_psth(:, 220:250), 2)) < 0.25;
    
        
    if method == 0
         these_units = passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2);
    elseif method == 1%method 1: high zscore
         these_units = passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2) & (abs(nanmean(zscore_psth(:, 260:300), 2)) > 0.25 | abs(nanmean(zscore_psth(:, 300:360), 2)) > 0.25);
    elseif method ==  2
    %method 2: mean shuffle
      these_units = passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2) &  passive_data.pvalue_shuffled_005{index}' == 1;
    elseif method ==  3
    % method 3: peak shuffle 
      these_units = passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2) &...
        (passive_data.pvalue_shuffled_0025_peak{index}' == 1 | passive_data.pvalue_shuffled_0025_trough{index}' == 1);
    end
    
    frac_down = sum(these_units & keep_these & nanmean(zscore_psth(:, 260:360), 2) < 0)/sum(passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2) & keep_these );

    frac_up = sum(these_units & keep_these & nanmean(zscore_psth(:, 260:360), 2) > 0)/sum(passive_data.unit_area == iRegion & ...
        (passive_data.unitType' == 1 | passive_data.unitType' == 2) & keep_these);
    

subplot(3, 1, iRegion)

barh([1,2],[frac_down.*100, frac_up.*100], 'FaceColor', regionColors{iRegion})


end

prettify_plot('XLimits', 'all')

%% example cells
% CP = [3427, 3488];
% cp_units = find(passive_data.unit_area==iRegion & sum(isnan(zscore_psth),2)<100);
%
% passive_data.animal_day_site_shank(cp_units(3427),:);