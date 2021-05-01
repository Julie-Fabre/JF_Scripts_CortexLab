
%% plot histology

animal = {'AP080'};
probeColors = {rgb('DeepPink'); rgb('Yellow'); rgb('Green'); rgb('Red'); rgb('OrangeRed'); ...
    rgb('Blue'); rgb('Lime'); rgb('Aqua'); rgb('LightSeaGreen'); rgb('Purple'); ...
    rgb('DarkKhaki'); rgb('OliveDrab'); rgb('DarkRed'); rgb('Orchid'); rgb('RosyBrown')}; %15 colors

%% plot slices overlaid with probe
%side by side: histology with probe one side, allen fit other side, ephys
%(average resp in 50um depth bins)
load(['\\znas.cortexlab.net\Subjects\', animal{:}, '\Histology\processed\probe2ephys.mat'])
load('\\znas.cortexlab.net\Subjects\AP080\Histology\processed\slices\histology_ccf.mat')
load('\\znas.cortexlab.net\Subjects\AP080\Histology\processed\slices\atlas2histology_tform.mat')
load('\\znas.cortexlab.net\Subjects\AP080\Histology\processed\slices\probe_ccf.mat')


imgs = dir(['\\znas.cortexlab.net\Subjects\', animal{:}, '\Histology\processed\slices\*.tif']);

% load all slices, and get info about which probe is where

for iImg = 1:length(imgs)
    %load image
    fullHistology(iImg, :, :, :) = imread_rgb([imgs(1).folder, filesep, imgs(iImg).name]);
    %draw probe(s) on image
    %probes = size(slice_slide_locations{1,iImg}{:},2);
    curr_av_slice = histology_ccf(iImg).av_slices;
    curr_av_slice(isnan(curr_av_slice)) = 1;
    iImg_im = fullHistology(iImg, :, :, :);

    tform = affine2d;
    h = imwarp(curr_av_slice, tform);

    for iProbe = 1:size(probe_ccf, 1)
        [idap, locap] = ismemberf(round(histology_ccf(1).plane_ap), round(probe_ccf(iProbe).trajectory_coords(:, 3)));
        [idml, locml] = ismemberf(histology_ccf(1).plane_ml, round(probe_ccf(iProbe).trajectory_coords(:, 2)));
        [iddv, locdv] = ismemberf(histology_ccf(1).plane_dv, round(probe_ccf(iProbe).trajectory_coords(:, 1)));

        [apFind1, apFind2] = find(locap);
        [mlFind1, mlFind2] = find(locml);
        [dvFind1, dvFind2] = find(locdv);
       
        
        if ~isempty(apFind1) && ~isempty(mlFind1) && ~isempty(dvFind1)
            %find common bits 
             [apMLfind, apMLLoc] = ismemberf(apFind1,mlFind1);
        [apdvfind, apdvLoc]= ismemberf(apFind1,dvFind1);
        [mldvfind, mldvLoc]= ismemberf(mlFind1,dvFind1);
        theseIndicesX = intersect(intersect(locap(find(apMLLoc)),locap(find(apdvLoc))),locml(find(mldvLoc)));
        [apMLfindy, apMLLocy] = ismemberf(apFind2,mlFind2);
        [apdvfindy, apdvLocy]= ismemberf(apFind2,dvFind2);
        [mldvfindy, mldvLocy]= ismemberf(mlFind2,dvFind2);
        theseIndicesY = intersect(intersect(locap(find(apMLLocy)),locap(find(apdvLocy))),locml(find(mldvLocy)));
        
        end
    end


end

% for each probe, one figure: side by side of slice, and allen, overlaid
% with probe, and probe trajectory + ephys (units/depth, aligned to visual stim for each experiment (overlaid)) on the side

  figure();
    subplot(321)
    set(gcf, 'color', 'k');
    imagesc(squeeze(fullHistology(iImg, :, :, :)))
    axis square;

    subplot(322)
    set(gcf, 'color', 'k');
    imagesc(squeeze(fullHistology(iImg, :, :, :)))
    hold on;


    allenBorders1 ...
        = [diff(curr_av_slice, [], 1) >  0; squeeze(zeros(1, size(curr_av_slice, 2)))];
    allenBorders1(allenBorders1 == 0) = NaN;

    allenBorders2 = [diff(curr_av_slice, [], 2) >  0, squeeze(zeros(size(curr_av_slice, 1), 1))];
    allenBorders2 = double(allenBorders2);
    allenBorders2(allenBorders2 == 0) = NaN;

    allenBorders = cat(3, allenBorders1, allenBorders2);
    allenBorders = nansum(allenBorders, 3);
    allenBordersAligned = imwarp(allenBorders, tform);
    allenBordersAligned(allenBordersAligned > 0) = 1;
    allenBordersAligned(allenBordersAligned == 0) = NaN;

    imagesc(allenBordersAligned, 'AlphaData', ~isnan(allenBordersAligned))

    subplot(323)
    imagesc(h)
    hold on;
    for iProbe = 1:size(probe_ccf, 1)
        for iPoint = 1:size(probe_ccf(iProbe).points,1)
            thisPoint = probe_ccf(iProbe).points(1,:);
            
        end
        [idap, locap] = ismemberf([round(histology_ccf(1).plane_ap), histology_ccf(1).plane_ml, histology_ccf(1).plane_dv], ...
            [round(probe_ccf(iProbe).points(:, 3)),  round(probe_ccf(iProbe).points(:, 2)),round(probe_ccf(iProbe).points(:, 1))]);
        [idml, locml] = ismemberf(histology_ccf(1).plane_ml, round(probe_ccf(iProbe).points(:, 2)));
        [iddv, locdv] = ismemberf(histology_ccf(1).plane_dv, round(probe_ccf(iProbe).points(:, 1)));

        [apFind1, apFind2] = find(locap);
        [mlFind1, mlFind2] = find(locml);
        [dvFind1, dvFind2] = find(locdv);
        if ~isempty(apFind1) && ~isempty(mlFind1) && ~isempty(dvFind1)
            hold on;

            line([mlFind2(1), mlFind2(end)], [dvFind2(1), dvFind2(end)], 'Color', 'k', 'LineWidth', 5);
        end
    end
figure();
    plot(theseIndicesX(1:size(theseIndicesY,2)),theseIndicesY)