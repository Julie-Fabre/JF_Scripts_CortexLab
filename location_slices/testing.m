theseNeurons = binX == iBinX & binY == iBinY & ...
    theseLocationsBregmaAbs(:, 2) >= round(chunks_region(iRegion, iChunk)) & ...
    theseLocationsBregmaAbs(:, 2) < round(chunks_region(iRegion, iChunk+1)) & ...
    passive_data.unit_area == iRegion & ...
    (passive_data.unitType' == 1 | passive_data.unitType' == 2) & passive_data.animal_thisDate_site_shank(:, 1) == iMouse;
for iChunk = 1:nChunks
for iBinX = 1:size(Xedges, 2)
                    for iBinY = 1:size(Yedges, 2)
theseNeurons = find(binX == iBinX & binY == iBinY & ...
    passive_data.unit_area == iRegion & ...
    (passive_data.unitType' == 1 | passive_data.unitType' == 2) &...
    passive_data.animal_thisDate_site_shank(:, 1) == iMouse)
if ~isempty(theseNeurons)
    keyboard;
end
                    end
end
end


theseNeurons = find(...
    passive_data.unit_area == iRegion & ...
    passive_data.animal_thisDate_site_shank(:, 1) == iMouse)


for curr_probe = 1:length(probe_ccf)
    % Plot points and line of best fit
    r0 = mean(probe_ccf(curr_probe).points,1);
    xyz = bsxfun(@minus,probe_ccf(curr_probe).points,r0);
    [~,~,V] = svd(xyz,0);
    histology_probe_direction = V(:,1);
    % (make sure the direction goes down in DV - flip if it's going up)
    if histology_probe_direction(2) < 0
        histology_probe_direction = -histology_probe_direction;
    end

    line_eval = [-1000,1000];
    probe_fit_line = bsxfun(@plus,bsxfun(@times,line_eval',histology_probe_direction'),r0);
    plot3(probe_ccf(curr_probe).points(:,1), ...
        probe_ccf(curr_probe).points(:,3), ...
        probe_ccf(curr_probe).points(:,2), ...
        '.','color',rgb('HotPink'),'MarkerSize',20);
    line(probe_fit_line(:,1),probe_fit_line(:,3),probe_fit_line(:,2), ...
        'color',gui_data.probe_color(curr_probe,:),'linewidth',2)
end