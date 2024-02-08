function [regionClassification, unitClassification] = cl_subsection_region(unit_coords, unit_area, postSpikeSup, waveformDur, propISI)

bregma = [540, 0, 570];
bregma_ml_point = bregma(1) / 2.5; %2.5 is difference in scaling between 
    % brainreg (25 um resolution) and allen (10um resolution, where this bregma value comes from)

%% ~~ striatum  ~~
%% classify striatal subregions
% striatum, seperate into dorsomedial, ventromedial, dorsolateral and posterior
% AP, DV, ML format
mediodorsal_limit = [ ...
    152, 155, 177; ...
    160, 155, 178; ...
    167, 148, 171; ...
    175, 148, 165; ...
    190, 142, 159; ...
    198, 135, 152; ...
    206, 128, 146; ...
    214, 128, 140; ...
    221, 128, 134; ...
    229, 128, 127; ...
    237, 128, 121; ...
    245, 128, 115; ...
    ];
mediodorsal_limit = mediodorsal_limit + [0, 15, 0];

posterior_limit = [245, 0, 0]; % AP only rule

striatum_units = find(unit_area == 1);

regionClassification = strings(length(unit_area), 1);

for iUnit = 1:length(striatum_units)
    unit = unit_coords(striatum_units(iUnit), :);
    
    if unit(1) > posterior_limit(1)
        regionClassification(striatum_units(iUnit)) = 'posterior_striatum';
        continue;
    end
    % Find closest AP value
    [~, closest_ap_idx] = min(abs(mediodorsal_limit(:, 1)-unit(1)));

    % Get the ML and DV limits for the closest AP
    ml_limit = mediodorsal_limit(closest_ap_idx, 2);
    dv_limit = mediodorsal_limit(closest_ap_idx, 3);

    % Classification based on ML and DV
    if unit(2) <= dv_limit && unit(3) >= ml_limit
        regionClassification(striatum_units(iUnit)) = 'dorsomedial_striatum';
    elseif unit(2) > dv_limit && bregma_ml_point - abs(unit(3) - bregma_ml_point) <= ml_limit
        regionClassification(striatum_units(iUnit)) = 'ventromedial_striatum';
    else
        regionClassification(striatum_units(iUnit)) = 'dorsolateral_striatum';
    end
end

% dorsomedial_count = sum(classification == 'dorsomedial_striatum');
% ventromedial_count = sum(classification == 'ventromedial_striatum');
% dorsolateral_count = sum(classification == 'dorsolateral_striatum');
% posterior_count = sum(classification == 'posterior_striatum');
%% classify striatal cell types 
paramEP.propISI_CP_threshold = 0.2;
paramEP.templateDuration_CP_threshold = 400;
paramEP.postSpikeSup_CP_threshold = 40;

unitClassification = strings(length(unit_area), 1);

unitClassification(waveformDur > paramEP.templateDuration_CP_threshold &...
    postSpikeSup < paramEP.postSpikeSup_CP_threshold) = 'MSN';

unitClassification(waveformDur <= paramEP.templateDuration_CP_threshold &...
    propISI <= paramEP.propISI_CP_threshold) = 'FSI';

unitClassification(waveformDur > paramEP.templateDuration_CP_threshold &...
    postSpikeSup >= paramEP.postSpikeSup_CP_threshold) = 'TAN';

unitClassification(waveformDur <= paramEP.templateDuration_CP_threshold &...
    propISI > paramEP.propISI_CP_threshold) = 'UIN';

%% GPe and SNr - use similar classif as striatum for now
  

end