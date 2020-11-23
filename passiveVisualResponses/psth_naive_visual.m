
%find == 'CP'
striatumData = find(arrayfun(@(x) contains(ephysData(x).location, 'CP'), 1:numel(ephysData)));

%find striatum limits - allen atlas
allen_atlas_path = 'C:\Users\Julie\Dropbox\Atlas\allenCCF';
tv = readNPY([allen_atlas_path, filesep, 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allen_atlas_path, filesep, 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allen_atlas_path, filesep, 'structure_tree_safe_2017.csv']); % a table of what all the labels mean
curr_plot_structure = find(contains(st.name, 'audoputamen'));
slice_spacing = 10;
plot_structure_color = hex2dec(reshape(st.color_hex_triplet{curr_plot_structure}, 2, [])') ./ 255;

%get limits and bins
binSize = 100;
structure_3d_lims = permute(av(1:slice_spacing:end, ...
    1:slice_spacing:end, 1:slice_spacing:end) == curr_plot_structure, [3, 1, 2]);
[r, c, v] = ind2sub(size(structure_3d_lims), find(structure_3d_lims));
x_lim_um = [min(r), max(r)] .* 100;
y_lim_um = [min(c), max(c)] .* 100;
z_lim_um = [min(v), max(v)] .* 100;
dataBinsX = x_lim_um(1):binSize:x_lim_um(2);
dataBinsY = y_lim_um(1):binSize:y_lim_um(2);
dataBinsZ = z_lim_um(1):binSize:z_lim_um(2);

%find which bins striatal templates belong to
uniqueRecordings = arrayfun(@(x) [ephysData(x).animal, ephysData(x).date, num2str(ephysData(x).site)], ...
    1:numel(ephysData), 'UniformOutput', false); %unique date + site
[uniqueId, uniqueIdx] = unique(uniqueRecordings, 'rows');
%psth - combining all protocols
posBinSize = 50;
timeBinSize = 0.01;
win = [0, 0.5];
bslWin = [-0.2, 0];

allP_allprotocols_allrecs = [];
posBinsX_allprotocols_allrecs = [];
posBinsY_allprotocols_allrecs = [];

for iUniqueRec = 1:size(uniqueId, 2)
    theseRecs = find(contains(uniqueRecordings, uniqueId(iUniqueRec)));
    allP_allprotocols = [];
    posBinsX_allprotocols = [];
    posBinsY_allprotocols = [];
    for iProtocol = 1:size(theseRecs, 2)

        spikeTimes = ephysData(theseRecs(iProtocol)).spike_times_timeline;
        uniqueTemp = unique(ephysData(theseRecs(iProtocol)).spike_templates);
        spikePos = nan(size(spikeTimes, 1), 3);
        for iTemp = 1:size(uniqueTemp, 1)
            theseSp = find(ephysData(theseRecs(iProtocol)).spike_templates == uniqueTemp(iTemp));
            spikePos(theseSp, 1) = ephysData(theseRecs(iProtocol)).template_location(iTemp, 1);
            spikePos(theseSp, 2) = ephysData(theseRecs(iProtocol)).template_location(iTemp, 2);
            spikePos(theseSp, 3) = ephysData(theseRecs(iProtocol)).template_location(iTemp, 3);
        end
        eventTimes = ephysData(theseRecs(iProtocol)).stimOn_times;

        if (max(spikePos(:, 1)) - min(spikePos(:, 1))) > posBinSize && (max(spikePos(:, 2)) - min(spikePos(:, 2))) > posBinSize
            [timeBins, posBinsX, posBinsY, allP, normVals] = psthByPos2D(spikeTimes, spikePos(:, 1), squeeze(spikePos(:, 2)), ...
                posBinSize, timeBinSize, eventTimes, win, bslWin);

        elseif (max(spikePos(:, 1)) - min(spikePos(:, 1))) > posBinSize
            [timeBins, posBinsX, allP, normVals] = psthByPos1D(spikeTimes, spikePos(:, 1), ...
                posBinSize, timeBinSize, eventTimes, win, bslWin);
            posBinsY = nan(size(posBinsX, 2), 1);
            posBinsY(:) = nanmean(spikePos(:, 2));

        elseif (max(spikePos(:, 2)) - min(spikePos(:, 2))) > posBinSize
            [timeBins, posBinsY, allP, normVals] = psthByPos1D(spikeTimes, squeeze(spikePos(:, 2)), ...
                posBinSize, timeBinSize, eventTimes, win, bslWin);
            posBinsX = nan(size(posBinsY, 2), 1);
            posBinsX(:) = nanmean(spikePos(:, 1));
        else
            [psth, ~, ~, ~, ~, ~] = psthAndBA(spikeTimes, eventTimes, bslWin, timeBinSize);
            normMn = mean(psth);
            normStd = std(psth);

            [psth, ~, ~, ~, ~, ~] = psthandBA(spikeTimes, ...
                eventTimes, win, timeBinSize);
            allP = (psth - normMn) ./ normStd;
            posBinsX = nanmean(spikePos(:, 1));
            posBinsY = nanmean(spikePos(:, 2));
        end
        allP_allprotocols = [allP_allprotocols, allP];
        posBinsX_allprotocols = [posBinsX_allprotocols, posBinsX];
        posBinsY_allprotocols = [posBinsY_allprotocols, posBinsY];
    end
    allP_allprotocols_allrecs = [allP_allprotocols_allrecs, allP_allprotocols];
    posBinsX_allprotocols_allrecs = [posBinsX_allprotocols_allrecs, posBinsX_allprotocols];
    posBinsY_allprotocols_allrecs = [posBinsY_allprotocols_allrecs, posBinsY_allprotocols];

end

%psth per grating, location, orientation, sp. frequency, nat. image

%cell type
