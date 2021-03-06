function AP_process_histologyJF(im_path, im_type)
% AP_process_histology(im_path)
%
% Resize and white balance histology images and extract images of each slice
% Andy Peters (peters.andrew.j@gmail.com)
% modified JF to include other formats. need matlab toolbox bioformats:
% https://www.openmicroscopy.org/bio-formats/downloads/
%
if strcmp(im_type, 't')
    % Get and sort image files
    im_path_dir = dir([im_path, filesep, '*.tif*']);
    im_fn = natsortfiles(cellfun(@(path, fn) [path, filesep, fn], ...
        {im_path_dir.folder}, {im_path_dir.name}, 'uni', false));

    % Get microns/pixel from metadata (if ome.tiff)
    im_info = imfinfo(im_fn{1});
    im_description = im_info(1).ImageDescription;
    im_um = regexp(im_description, 'PhysicalSizeX="(\S*)".*PhysicalSizeY="(\S*)"', 'tokens');
    im_um_x = str2num(im_um{1}{1});
    im_um_y = str2num(im_um{1}{2});
elseif strcmp(im_type, 'brainSaw') == 1
    im_path_dir = dir([im_path, filesep, '*Composite-Substack.tif*']);
    im_fn = natsortfiles(cellfun(@(path, fn) [path, filesep, fn], ...
        {im_path_dir.folder}, {im_path_dir.name}, 'uni', false));
    im_fn = im_fn(end:-1:1);
    im_info = imfinfo(im_fn{1});
    im_um_x = 25; %/10000 for cm -> um
    im_um_y = 25;
    
elseif strcmp(im_type, 'ot') == 1
    im_path_dir = dir([im_path, filesep, '*.tif*']);
    im_fn = natsortfiles(cellfun(@(path, fn) [path, filesep, fn], ...
        {im_path_dir.folder}, {im_path_dir.name}, 'uni', false));

    % Get microns/pixel from metadata (if ome.tiff)
    im_info = imfinfo(im_fn{1});
    % im_description = im_info(1).UnknownTags.Value;
    %im_um = regexp(im_description,'PhysicalSizeX="(\S*)".*PhysicalSizeY="(\S*)"','tokens');
    %     im_um_x = im_info(1).XResolution/10000; %/10000 for cm -> um
    %     im_um_y = im_info(1).YResolution/10000;
    %im_um = im_info(1).ResolutionUnit;
    im_um_x = 3.45;
    im_um_y = 3.45; %QQ hardcoded BAD
elseif strcmp(im_type, 'tiffUnmerged')
    im_path_dir = dir([im_path, filesep, '*.tif*']);
    im_fn = natsortfiles(cellfun(@(path, fn) [path, filesep, fn], ...
        {im_path_dir.folder}, {im_path_dir.name}, 'uni', false));

    % Get microns/pixel from metadata (if ome.tiff)
    im_info = imfinfo(im_fn{1});
    % im_description = im_info(1).UnknownTags.Value;
    %im_um = regexp(im_description,'PhysicalSizeX="(\S*)".*PhysicalSizeY="(\S*)"','tokens');
    %     im_um_x = im_info(1).XResolution/10000; %/10000 for cm -> um
    %     im_um_y = im_info(1).YResolution/10000;
    %im_um = im_info(1).ResolutionUnit;
    im_um_x = 3.45;
    im_um_y = 3.45; %QQ hardcoded BAD
elseif strcmp(im_type, 'tiffUnmergedNoDef')
    im_path_dir = dir([im_path, filesep, '*.tif*']);
    im_fn = natsortfiles(cellfun(@(path, fn) [path, filesep, fn], ...
        {im_path_dir.folder}, {im_path_dir.name}, 'uni', false));

    % Get microns/pixel from metadata (if ome.tiff)
    im_info = imfinfo(im_fn{1});
    im_um_x = 2.3; %/10000 for cm -> um
    im_um_y = 2.3;
elseif strcmp(im_type, 'n')
    im_path_dir = dir([im_path, filesep, '*.nd2*']);
    im_fn = natsortfiles(cellfun(@(path, fn) [path, filesep, fn], ...
        {im_path_dir.folder}, {im_path_dir.name}, 'uni', false));
    %dd=dir([im_path filesep '*nd2']);
    reader = bfGetReader([im_fn{1}]);
    omeMeta = reader.getMetadataStore();
    im_um_x = double(omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER));
    im_um_y = double(omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROMETER));
    im_info = struct;
    for iImg = 1:numel(im_fn)
        reader = bfGetReader([im_fn{iImg}]);
        omeMeta = reader.getMetadataStore();
        im_info(iImg).Height = double(omeMeta.getPixelsSizeY(0).getValue());
        im_info(iImg).Width = double(omeMeta.getPixelsSizeX(0).getValue());
    end
    im_path_dir = dir([im_path, filesep, '*.tif*']);
    im_fn = natsortfiles(cellfun(@(path, fn) [path, filesep, fn], ...
        {im_path_dir.folder}, {im_path_dir.name}, 'uni', false));
end
if im_um_x ~= im_um_y
    error('Pixel X/Y values different')
end

% Get number of imaged channels (ome.tiff has extra page with nothing?)
%if strcmp(im_type, 'brainSaw') == 1
%    n_channels = length(im_fn);
%else
    n_channels = sum(any([im_info.Height; im_info.Width], 1));
%end
%n_channels =2; %QQ hard coded here, change
% Set resize factor to match to Allen CCF
allen_um2px = 10; % Allen CCF: 10 um/voxel
im_rescale_factor = im_um_x / allen_um2px;

% Load and resize images
n_im = length(im_fn);
im_resized = cell(n_im, n_channels);

%h = waitbar(0,'Loading and resizing images...');
% imC = imread(im_fn{2},'tif',2);
% size(imC)
% figure();
% imshow(imC)
% %
% imo = imread(im_fn{1},'tif', curr_channel);
% size(imo)
% figure();
% imshow(imo)
%
% im2 = imread(im_fn{3},'tif', 2);
% size(im2)
% figure();
% imshow(im2)
if strcmp(im_type, 'tiffUnmerged')
    n_channels = 3; %QQ hard-coded
    for curr_im = 1:n_im / 3
        for curr_channel = 1:3
            im_resized{curr_im, curr_channel, 1} = imresize(imread(im_fn{(curr_im - 1)*3+curr_channel}, 'tif'), im_rescale_factor);
            imwrite(im_resized{curr_im, curr_channel, 1}, [im_path, '/resized/', strcat(num2str(curr_im), '_', num2str(curr_channel), '.tiff')], 'tiff');
            nameS{:, curr_channel} = [im_path, '/resized/', strcat(num2str(curr_im), '_', num2str(curr_channel), '.tiff')];
        end
        %waitbar(curr_im/n_im,h,['Loading and resizing images (' num2str(curr_im) '/' num2str(n_im) ')...']);
    end
else
    if strcmp(im_type, 'brainSaw')
        n_im = 2;
        %n_channels=2;
    end
    %if isempty(dir([im_path, '/resized/*tiff']))
    for curr_im = 1:n_im
        for curr_channel = 1:n_channels
            
            
            if strcmp(im_type, 'brainSaw') == 1
                %if size(dir([im_path, '/resized/']),1) <= 2
                im_resized{curr_channel, curr_im, 1} = imresize(imread(im_fn{curr_im}, 'tif', curr_channel), im_rescale_factor);
                imwrite(im_resized{curr_channel, curr_im, 1}, [im_path, '/resized/', strcat(num2str(curr_channel), '_', num2str(curr_im), '.tiff')], 'tiff');
                
                %end
                nameS{:, curr_channel} = [im_path, '/resized/', strcat(num2str(curr_channel), '_', num2str(curr_im), '.tiff')];
            else
                im_resized{curr_im, curr_channel, 1} = imresize(imread(im_fn{curr_im}, 'tif', curr_channel), im_rescale_factor);
                imwrite(im_resized{curr_im, curr_channel, 1}, [im_path, '/resized/', strcat(num2str(curr_im), '_', num2str(curr_channel), '.tiff')], 'tiff');
                nameS{:, curr_channel} = [im_path, '/resized/', strcat(num2str(curr_im), '_', num2str(curr_channel), '.tiff')];
            end

        end
        %waitbar(curr_im/n_im,h,['Loading and resizing images (' num2str(curr_im) '/' num2str(n_im) ')...']);
   % end
    end
end
%close(h);


% Estimate white balance within each channel
% (dirty: assume one peak for background, one for signal)
h = figure;
if ~strcmp(im_type, 'brainSaw') == 1
   
    
im_montage = cell(n_channels, 1);
channel_caxis = nan(n_channels, 2);
channel_color = cell(n_channels, 1);
end
if n_channels == 2
    chC = {'blue', 'red'};
elseif n_channels == 3
    if strcmp(im_type, 'brainSaw') == 1
        chC = {'red', 'green', 'blue'};
        n_im = n_channels;
        n_channels = 3;
        
    else
        chC = {'blue', 'red', 'green'};
    end
end
if strcmp(im_type, 'brainSaw') == 1
        chC = {'red', 'green', 'blue'};
        n_im = n_channels;
        n_channels = 2;
        im_montage = cell(n_channels, 1);
channel_caxis = nan(n_channels, 2);
channel_color = cell(n_channels, 1);
else
end
for curr_channel = 1:n_channels
    %image(im_resized{1,curr_channel})
    hh = figure();
    curr_montage = montage(nameS(1, curr_channel));

    im_montage{curr_channel} = curr_montage.CData;

    im_hist = histcounts(im_montage{curr_channel}(im_montage{curr_channel} > 0), 0:max(im_montage{curr_channel}(:)));
    im_hist_deriv = diff(smooth(im_hist, 10));

    [~, bg_median] = max(im_hist);
    bg_signal_min = find(im_hist_deriv(1:end) > 0, 1); %-1 b/c of diff
    [~, bg_median_rel] = max(im_hist(1:end));
    signal_median = bg_median_rel + bg_median;

    cmin = min(min(curr_montage.CData));
    cmax = max(max(curr_montage.CData));
    caxis([cmin, cmax]);

    channel_caxis(curr_channel, :) = [cmin, cmax];

    channel_color{curr_channel} = chC{curr_channel};

end
close(h)

% Get order of colors
color_order_gun = { 'green'; 'blue'; 'red'};
[~, color_order_slide] = ismember(channel_color, color_order_gun);
%color_order_slide(color_order_slide==0)=[];
% Display montage of final balanced image, sort color channels by RGB
im_montage_rgb = zeros(size(im_montage{1}, 1), size(im_montage{1}, 2), 3);
im_montage_rgb(:, :, color_order_slide) = ...
    cell2mat(arrayfun(@(ch) rescale(im_montage{ch}, ...
    'InputMin', channel_caxis(ch, 1), 'InputMax', channel_caxis(ch, 2)), ...
    permute(1:n_channels, [1, 3, 2]), 'uni', false));
figure;
imshow(im_montage_rgb);
title('Overview of all images');

% Store RGB for each slide
im_rgb = cellfun(@(x) zeros(size(x, 1), size(x, 2), 3), im_resized(:, 1), 'uni', false);
for curr_im = 1:n_im
    im_rgb{curr_im}(:, :, color_order_slide) = ...
        cell2mat(arrayfun(@(ch) rescale(im_resized{curr_im, ch}, ...
        'InputMin', channel_caxis(ch, 1), 'InputMax', channel_caxis(ch, 2)), ...
        permute(1:n_channels, [1, 3, 2]), 'uni', false));
end

% Set up GUI to pick slices on slide to extract
slice_fig = figure('KeyPressFcn', @slice_keypress);

% Initialize data
slice_data = struct;
slice_data.im_path = im_path;
slice_data.im_fn = im_fn;
slice_data.im_rescale_factor = im_rescale_factor;
slice_data.im_rgb = im_rgb;
slice_data.curr_slide = 0;
slice_data.slice_mask = cell(0, 0);
slice_data.slice_rgb = cell(0, 0);

% Update gui data
guidata(slice_fig, slice_data);

% Update slide
update_slide(slice_fig);

end

function slice_click(slice_fig, eventdata)
% On slice click, mark to extract

slice_data = guidata(slice_fig);

selected_slice_bw = bwselect(slice_data.mask, eventdata.IntersectionPoint(1), eventdata.IntersectionPoint(2));

if eventdata.Button == 1
    % If left button pressed, create new slice ROI
    roi_num = size(slice_data.user_masks, 3) + 1;

    % Make new mask with object
    slice_data.user_masks(:, :, roi_num) = selected_slice_bw;

    % Draw bounding box around object
    box_x = find(any(slice_data.user_masks(:, :, roi_num), 1), 1);
    box_y = find(any(slice_data.user_masks(:, :, roi_num), 2), 1);
    box_w = find(any(slice_data.user_masks(:, :, roi_num), 1), 1, 'last') - box_x;
    box_h = find(any(slice_data.user_masks(:, :, roi_num), 2), 1, 'last') - box_y;
    slice_data.user_rectangles(roi_num) = ...
        rectangle('Position', [box_x, box_y, box_w, box_h], 'EdgeColor', 'w');

elseif eventdata.Button == 3
    % if right button pressed, join to last slice ROI
    roi_num = size(slice_data.user_masks, 3);

    % Join old and new objects in mask
    slice_data.user_masks(:, :, roi_num) = ...
        slice_data.user_masks(:, :, roi_num) | selected_slice_bw;

    % Draw bounding box around object
    box_x = find(any(slice_data.user_masks(:, :, roi_num), 1), 1);
    box_y = find(any(slice_data.user_masks(:, :, roi_num), 2), 1);
    box_w = find(any(slice_data.user_masks(:, :, roi_num), 1), 1, 'last') - box_x;
    box_h = find(any(slice_data.user_masks(:, :, roi_num), 2), 1, 'last') - box_y;
    set(slice_data.user_rectangles(roi_num), 'Position', ...
        [box_x, box_y, box_w, box_h]);

end

% Update gui data
guidata(slice_fig, slice_data);

end

function slice_keypress(slice_fig, eventdata)
% Move to next slide with spacebar

if strcmp(eventdata.Key, 'space')
    update_slide(slice_fig)
end

end

function update_slide(slice_fig)
% Find slices on slide by over-threshold objects of a large enough size

slice_data = guidata(slice_fig);

% Pull the images from selected slices (not during initialization)
if slice_data.curr_slide > 0
    extract_slice_rgb(slice_fig);
    slice_data = guidata(slice_fig);
end

% After the last slice, save the images and close out
%if strcmp(im_type, 'tiffUnmergedNoDef')
    if slice_data.curr_slide == floor(length(slice_data.im_rgb))
        save_slice_rgb(slice_fig);
        close(slice_fig);
        return
    end
%else
%      if slice_data.curr_slide == floor(length(slice_data.im_rgb)/3)
%         save_slice_rgb(slice_fig);
%         close(slice_fig);
%         return
%     end
% end
% else
%     if slice_data.curr_slide == length(slice_data.im_rgb) %%QQ UNCOMMENYT IF
%     NOT UNMERGED
%         save_slice_rgb(slice_fig);
%         close(slice_fig);
%         return
%     end
% end

slice_data.curr_slide = slice_data.curr_slide + 1;

min_slice = (1000 / 10)^2; % (um/10(CCF units))^2
slice_mask = imfill(bwareaopen(mean( ...
    slice_data.im_rgb{slice_data.curr_slide}, 3) > 0.01, min_slice), 'holes');
slice_conncomp = bwconncomp(slice_mask);

im_handle = imshow(slice_data.im_rgb{slice_data.curr_slide});
set(im_handle, 'ButtonDownFcn', @slice_click);
title('Finding slice boundaries...');
drawnow;

slice_boundaries = bwboundaries(slice_mask);
slice_lines = gobjects(length(slice_boundaries), 1);
for curr_slice = 1:length(slice_boundaries)
    slice_lines(curr_slice) = line(slice_boundaries{curr_slice}(:, 2), ...
        slice_boundaries{curr_slice}(:, 1), 'color', 'w', 'linewidth', 2, 'LineSmoothing', 'on', 'linestyle', '--');
end
title('Click slices to extract (left = new, right = add to last), spacebar to finish slide');

slice_data.im_h = im_handle;
slice_data.mask = slice_mask;
slice_data.lines = slice_lines;
slice_data.user_masks = zeros(size(slice_mask, 1), size(slice_mask, 2), 0, 'logical');
slice_data.user_rectangles = gobjects(0);

% Update gui data
guidata(slice_fig, slice_data);

end


function extract_slice_rgb(slice_fig)
% When changing slide, extract the selected slice images

slice_data = guidata(slice_fig);

n_slices = size(slice_data.user_masks, 3);
curr_slice_mask = cell(n_slices, 1);
curr_slice_rgb = cell(n_slices, 1);
for curr_slice = 1:n_slices
    % Pull a rectangular area, exclude spaces (e.g. between torn piece)
    dilate_size = 30;
    curr_mask = imdilate(logical(any(slice_data.user_masks(:, :, curr_slice), 2).* ...
        any(slice_data.user_masks(:, :, curr_slice), 1)), ones(dilate_size));

    curr_rgb = reshape(slice_data.im_rgb{slice_data.curr_slide}( ...
        repmat(curr_mask, 1, 3)), sum(any(curr_mask, 2)), sum(any(curr_mask, 1)), 3);

    curr_slice_mask{curr_slice} = curr_mask;
    curr_slice_rgb{curr_slice} = curr_rgb;

end

% Store the image and mask for each slice
slice_data.slice_mask{slice_data.curr_slide} = curr_slice_mask;
slice_data.slice_rgb{slice_data.curr_slide} = curr_slice_rgb;

% Update gui data
guidata(slice_fig, slice_data);

end


function save_slice_rgb(slice_fig)
% After the last slide, save the slice images

slice_data = guidata(slice_fig);

% Set save directory as subdirectory within original
save_dir = [slice_data.im_path, filesep, 'slices'];
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end

% Concatenate all slice images
slice_rgb_cat = vertcat(slice_data.slice_rgb{:});

% Write all slice images to separate files
for curr_im = 1:length(slice_rgb_cat)
    curr_fn = [save_dir, filesep, num2str(curr_im), '.tif'];
    imwrite(slice_rgb_cat{curr_im}, curr_fn, 'tif');
end

% Get rows and columns for each slice corresponding to full size image
slice_slide_locations = cell(size(slice_data.slice_mask));
for curr_slide = 1:length(slice_data.slice_mask)
    for curr_slice = 1:length(slice_data.slice_mask{curr_slide})

        curr_mask = slice_data.slice_mask{curr_slide}{curr_slice};

        mask_x = find(interp1(1:size(curr_mask, 2), +any(curr_mask, 1), ...
            linspace(1, size(curr_mask, 2), ...
            round(size(curr_mask, 2)/slice_data.im_rescale_factor)), 'nearest'));
        mask_y = find(interp1(1:size(curr_mask, 1), +any(curr_mask, 2), ...
            linspace(1, size(curr_mask, 1), ...
            round(size(curr_mask, 1)/slice_data.im_rescale_factor)), 'nearest'));

        slice_slide_locations{curr_slide}{curr_slice} = ...
            {mask_y, mask_x};

    end
end

slice_slide_locations_fn = [save_dir, filesep, 'slice_slide_locations.mat'];
save(slice_slide_locations_fn, 'slice_slide_locations');

disp(['Slices saved in ', save_dir]);

end
