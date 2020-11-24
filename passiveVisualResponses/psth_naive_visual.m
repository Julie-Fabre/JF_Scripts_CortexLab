
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
win = [-0.2, 0.5];
bslWin = [];% [-0.2, 0];

allP_allprotocols_allrecs = [];
posBinsX_allprotocols_allrecs = [];
posBinsY_allprotocols_allrecs = [];
posBinsZ_allprotocols_allrecs = [];
for iUniqueRec = 1:size(uniqueId, 2)
    theseRecs = find(contains(uniqueRecordings, uniqueId(iUniqueRec)));
    allP_allprotocols = [];
    posBinsX_allprotocols = [];
    posBinsY_allprotocols = [];
    posBinsZ_allprotocols = [];
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
        %disp(nanmean(nanmean(spikePos)))
        if (max(spikePos(:, 1)) - min(spikePos(:, 1))) > posBinSize && (max(spikePos(:, 3)) - min(spikePos(:, 3))) > posBinSize
            [timeBins, posBinsX, posBinsY, posBinsZ, allP, normVals] = psthByPos3D(spikeTimes, spikePos(:, 1), squeeze(spikePos(:, 3)), ...
                squeeze(spikePos(:, 2)), ...
                posBinSize, timeBinSize, eventTimes, win, bslWin);
            posBinsX(1:end-1) = nanmean([posBinsX(1:end-1); posBinsX(2:end)]);
            posBinsX(end) = [];
            posBinsY(1:end-1) = nanmean([posBinsY(1:end-1); posBinsY(2:end)]);
            posBinsY(end) = [];
            posBinsZ(1:end-1) = nanmean([posBinsZ(1:end-1); posBinsZ(2:end)]);
            posBinsZ(end) = [];
        elseif (max(spikePos(:, 1)) - min(spikePos(:, 1))) > posBinSize
            [timeBins, posBinsX, posBinsZ, allP, normVals] = psthByPos2D(spikeTimes, spikePos(:, 1), ...
                squeeze(spikePos(:, 2)), posBinSize, timeBinSize, eventTimes, win, bslWin);
            posBinsY = nan(1, size(posBinsX, 2)-1);
            posBinsY(:) = nanmean(spikePos(:, 3));
            posBinsX(1:end-1) = nanmean([posBinsX(1:end-1); posBinsX(2:end)]);
            posBinsX(end) = [];
            posBinsZ(1:end-1) = nanmean([posBinsZ(1:end-1); posBinsZ(2:end)]);
            posBinsZ(end) = [];
        elseif (max(spikePos(:, 3)) - min(spikePos(:, 3))) > posBinSize
            [timeBins, posBinsY, posBinsZ, allP, normVals] = psthByPos2D(spikeTimes, squeeze(spikePos(:, 2)), ...
                squeeze(spikePos(:, 3)), posBinSize, timeBinSize, eventTimes, win, bslWin);
            posBinsX = nan(1, size(posBinsY, 2)-1);
            posBinsX(:) = nanmean(spikePos(:, 1));
            posBinsY(1:end-1) = nanmean([posBinsY(1:end-1); posBinsY(2:end)]);
            posBinsY(end) = [];
            posBinsZ(1:end-1) = nanmean([posBinsZ(1:end-1); posBinsZ(2:end)]);
            posBinsZ(end) = [];
        else
            [timeBins, posBinsZ, allP, normVals] = psthByPos1D(spikeTimes, ...
                squeeze(spikePos(:, 2)), posBinSize, timeBinSize, eventTimes, win, bslWin);

            posBinsZ(1:end-1) = nanmean([posBinsZ(1:end-1); posBinsZ(2:end)]);
            posBinsZ(end) = [];
            posBinsX = nan(1, size(posBinsZ, 2));
            posBinsX(:) = nanmean(spikePos(:, 1));
            posBinsY = nan(1, size(posBinsZ, 2));
            posBinsY(:) = nanmean(spikePos(:, 3));
        end

        allP_allprotocols = [allP_allprotocols; allP];
        posBinsX_allprotocols = [posBinsX_allprotocols, round(posBinsX/posBinSize) * posBinSize];
        posBinsY_allprotocols = [posBinsY_allprotocols, round(posBinsY/posBinSize) * posBinSize];
        posBinsZ_allprotocols = [posBinsZ_allprotocols, round(posBinsZ/posBinSize) * posBinSize];
if iProtocol == size(theseRecs, 2)
figure(); imagesc(allP_allprotocols)
end
    end
    allP_allprotocols_allrecs = [allP_allprotocols_allrecs; allP_allprotocols];
    posBinsX_allprotocols_allrecs = [posBinsX_allprotocols_allrecs, posBinsX_allprotocols];
    posBinsY_allprotocols_allrecs = [posBinsY_allprotocols_allrecs, posBinsY_allprotocols];
    posBinsZ_allprotocols_allrecs = [posBinsZ_allprotocols_allrecs, posBinsZ_allprotocols];
end

figure(); %position of recordings (binned by posBinSize)
scatter3(posBinsX_allprotocols_allrecs, posBinsY_allprotocols_allrecs, posBinsZ_allprotocols_allrecs)

figure();
imagesc(allP_allprotocols_allrecs);


%psth per grating, location, orientation, sp. frequency, nat. image
%average same

[uniqueXY, uniqueXY_id] = unique([posBinsX_allprotocols_allrecs', posBinsY_allprotocols_allrecs'], 'rows');
uniqueX = unique([posBinsX_allprotocols_allrecs]);
uniqueXrange = min(uniqueX):posBinSize:max(uniqueX);
uniqueY = unique([posBinsY_allprotocols_allrecs]);
uniqueYrange = min(uniqueY):posBinSize:max(uniqueY);
for iUniqueX = 1:length(uniqueXrange)
    for iUniqueY = 1:length(uniqueYrange)
        thisVal = nanmean(nanmax(allP_allprotocols_allrecs(posBinsX_allprotocols_allrecs == uniqueXrange(iUniqueX) & ...
            posBinsY_allprotocols_allrecs == uniqueYrange(iUniqueY), 1:20)));
        if ~isempty(thisVal)
            allP_meanUniK_XY(iUniqueX, iUniqueY) = nanmean(nanmax(thisVal));
        else
            allP_meanUniK_XY(iUniqueX, iUniqueY) = NaN;
        end

    end
end
figure();
b1 = imagesc(uniqueXrange, uniqueYrange, allP_meanUniK_XY);
set(b1, 'AlphaData', ~isnan(allP_meanUniK_XY) & allP_meanUniK_XY ~= 0)
colorbar;

uniqueZ = unique([posBinsZ_allprotocols_allrecs]);
uniqueZrange = min(uniqueZ):posBinSize:max(uniqueZ);
for iUniqueX = 1:length(uniqueXrange)
    for iUniqueZ = 1:length(uniqueZrange)
       thisVal = nanmean(nanmax(allP_allprotocols_allrecs(posBinsX_allprotocols_allrecs == uniqueXrange(iUniqueX) & ...
            posBinsZ_allprotocols_allrecs == uniqueZrange(iUniqueZ), 1:20)));
        if ~isempty(thisVal)
            allP_meanUniK_XZ(iUniqueX, iUniqueZ) = nanmean(nanmax(thisVal));
        else
            allP_meanUniK_XZ(iUniqueX, iUniqueZ) = NaN;
        end
    end
end
figure();
b2 = imagesc(uniqueXrange, uniqueZrange, allP_meanUniK_XZ);
set(b2, 'AlphaData', ~isnan(allP_meanUniK_XZ) & allP_meanUniK_XZ ~= 0)
colorbar;

for iUniqueY = 1:length(uniqueYrange)
    for iUniqueZ = 1:length(uniqueZrange)
        thisVal = nanmean(nanmax(allP_allprotocols_allrecs(posBinsY_allprotocols_allrecs == uniqueYrange(iUniqueY) & ...
            posBinsZ_allprotocols_allrecs == uniqueZrange(iUniqueZ), 1:20)));
        if ~isempty(thisVal)
            allP_meanUniK_YZ(iUniqueY, iUniqueZ) = nanmean(nanmax(thisVal));
        else
            allP_meanUniK_YZ(iUniqueY, iUniqueZ) = NaN;
        end
    end
end
figure();
b = imagesc(uniqueYrange, uniqueZrange, allP_meanUniK_YZ);
set(b, 'AlphaData', ~isnan(allP_meanUniK_YZ) & allP_meanUniK_YZ ~= 0)
colorbar;
%cell type
% tod o replace zweros by nans, plot all *grey*