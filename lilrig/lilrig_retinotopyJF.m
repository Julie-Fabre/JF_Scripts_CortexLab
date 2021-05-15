% lilrig_retinotopy
%
% Generates visual field sign map from retinotopy on lilrig
% stim_program - 'mpep' or 'signals'

if exist('Protocol','var')
    if strcmp(Protocol.xfile,'stimSparseNoiseUncorrAsync.x')
        stim_program = 'mpep';
    else
        error('Unknown MPEP retinotopy protocol');
    end
elseif exist('expDef','var')
    if strcmp(expDef,'AP_sparseNoise') || strcmp(expDef,'sparseNoiseAsync_NS2')
        stim_program = 'signals';
    else
        error('Unknown signals retinotopy expDef');
    end
end

%% Get photodiode flip times

% Threshold the photodiode trace, find flips
photodiode_thresh = 3;
photodiode_trace = Timeline.rawDAQData(stimScreen_on,photodiode_idx) > photodiode_thresh;
% (medfilt because photodiode can be intermediate value when backlight
% coming on)
photodiode_trace_medfilt = medfilt1(Timeline.rawDAQData(stimScreen_on, ...
    photodiode_idx),3) > photodiode_thresh;
photodiode_flip = find((~photodiode_trace_medfilt(1:end-1) & photodiode_trace_medfilt(2:end)) | ...
    (photodiode_trace_medfilt(1:end-1) & ~photodiode_trace_medfilt(2:end)))+1;
photodiode_flip_times = stimScreen_on_t(photodiode_flip)';


%% Get stimulus squares and times (protocol-dependent)

switch stim_program
    
    %% MPEP sparse noise retinotopy
    case 'mpep'
                
        % Generate the sparse noise stimuli from the protocol
        myScreenInfo.windowPtr = NaN; % so we can call the stimulus generation and it won't try to display anything
        stimNum = 1;
        ss = eval([Protocol.xfile(1:end-2) '(myScreenInfo, Protocol.pars(:,stimNum));']);
        stim_screen = cat(3,ss.ImageTextures{:});
        ny = size(stim_screen,1);
        nx = size(stim_screen,2);
        
        switch lower(photodiode_type)
            case 'flicker'
                
                if size(stim_screen,3) == length(photodiode_flip_times)
                    % If stim matches photodiode, use directly
                    stim_times = photodiode_flip_times;
                    
                elseif mod(size(stim_screen,3),2) == 1 && ...
                        length(photodiode_flip_times) == size(stim_screen,3) + 1
                    % Check for case of mismatch between photodiode and stimuli:
                    % odd number of stimuli, but one extra photodiode flip to come back down
                    photodiode_flip_times(end) = [];
                    stim_times = photodiode_flip_times;
                    warning('Odd number of stimuli, removed last photodiode');
                    
                elseif size(stim_screen,3) ~= length(photodiode_flip_times)
                    % If there's a different kind of mismatch, guess stim times
                    % by interpolation
                    photodiode_flip_times = photodiode_flip_times([1,end]);
                    stim_duration = diff(photodiode_flip_times)/size(stim_screen,3);
                    stim_times = linspace(photodiode_flip_times(1), ...
                        photodiode_flip_times(2)-stim_duration,size(stim_screen,3))';
                    warning('Mismatching stim and photodiode, interpolating start/end')
                end
                
            case 'steady'
                % If the photodiode is on steady: extrapolate the stim times
                if length(photodiode_flip_times) ~= 2
                    error('Steady photodiode, but not 2 flips')
                end
                stim_duration = diff(photodiode_flip_times)/size(stim_screen,3);
                stim_times = linspace(photodiode_flip_times(1), ...
                    photodiode_flip_times(2)-stim_duration,size(stim_screen,3))';
                
        end
        
    case 'signals'
                
        ny = size(block.events.stimuliOnValues,1);
        nx = size(block.events.stimuliOnValues,2)/ ...
            size(block.events.stimuliOnTimes,2);
        
        % Get stim (not firsrt: initializes to black on startup)
        stim_screen = reshape(block.events.stimuliOnValues(:,nx+1:end),ny,nx,[]);      
        
        % Each photodiode flip is a screen update
        stim_times = photodiode_flip_times;
        
        % (if more stim than times, just try matching the last n stim)
        if size(stim_screen,3) > length(stim_times)
            warning('More stims than photodiode flips - truncating beginning')
            stim_screen(:,:,1:(size(stim_screen,3)-length(stim_times))) = [];
        end
        
        % (if times than stim, just try using the last n times)
        if size(stim_screen,3) < length(stim_times)
            warning('More photodiode flips than stim - truncating beginning')
            stim_times(1:(length(stim_times)-size(stim_screen,3))) = [];
        end
        
end

% Check that photodiode times match stim number
if size(stim_screen,3) ~= length(stim_times)
   error('Mismatching stim number and photodiode times'); 
end


%% Get average response to each stimulus (bootstrap mean)

surround_window = [0.3,0.5]; % 6s = [0.3,0.5]
%frame_t=block.events.stimuliOnTimes;
if ~exist('frame_t', 'var')
    th = cam_time(2:2:end);
    frame_t = th;
end

framerate = 1./nanmedian(diff(frame_t));
surround_samplerate = 1/(framerate*1);
surround_time = surround_window(1):surround_samplerate:surround_window(2);
response_n = nan(ny,nx);
response_grid = cell(ny,nx);
for px_y = 1:ny
    for px_x = 1:nx
        
        switch stim_program
            case 'mpep'
                % Gray to either black or white: use either
                align_stims = (stim_screen(px_y,px_x,2:end)~= 0) & ...
                    (diff(stim_screen(px_y,px_x,:),[],3) ~= 0);
                align_times = stim_times(find(align_stims)+1);
            case 'signals'
                % Black to white
                align_stims = stim_screen(px_y,px_x,2:end) == 1 & ...
                    stim_screen(px_y,px_x,1:end-1) == -1;
                align_times = stim_times(find(align_stims)+1);
        end
        
        response_n(px_y,px_x) = length(align_times);
        
        % Don't use times that fall outside of imaging
        align_times(align_times + surround_time(1) < frame_t(2) | ...
            align_times + surround_time(2) > frame_t(end)) = [];
        
        % Get stim-aligned responses, 2 choices:
        
        % 1) Interpolate times (slow - but supersamples so better)
        %         align_surround_times = bsxfun(@plus, align_times, surround_time);
        %         peri_stim_v = permute(mean(interp1(frame_t,fV',align_surround_times),1),[3,2,1]);
        
        % 2) Use closest frames to times (much faster - not different)
        align_surround_times = align_times + surround_time;
        frame_edges = [frame_t,frame_t(end)+1/framerate];
        align_frames = discretize(align_surround_times,frame_edges);
                
        % Get stim-aligned baseline (at stim onset)
        align_baseline_times = align_times;
        align_frames_baseline = discretize(align_baseline_times,frame_edges);
        
        % Don't use NaN frames (delete, dirty)
        nan_stim = any(isnan(align_frames),2) | isnan(align_frames_baseline);     
        align_frames(nan_stim,:) = [];
        align_frames_baseline(nan_stim,:) = [];
     
        % Define the peri-stim V's as subtracting first frame (baseline)
        peri_stim_v = ...
            reshape(fV(:,align_frames)',size(align_frames,1),size(align_frames,2),[]) - ...
            nanmean(reshape(fV(:,align_frames_baseline)',size(align_frames_baseline,1),size(align_frames_baseline,2),[]),2);       
        
        mean_peri_stim_v = permute(mean(peri_stim_v,2),[3,1,2]);
        
        % Save V's
        response_grid{px_y,px_x} = mean_peri_stim_v;
        
    end
end

% Get position preference for every pixel
U_downsample_factor = 1; %2 if max method
screen_resize_scale = 1; %3 if max method
filter_sigma = (screen_resize_scale*2);

% Downsample U
[Uy,Ux,nSV] = size(U);
use_u_y = 1:Uy;
Ud = imresize(U(use_u_y,:,:),1/U_downsample_factor,'bilinear');

% Convert V responses to pixel responses
use_svs = 1:500; % de-noises, otherwise size(U,3)
n_boot = 10;

response_mean_bootstrap = cellfun(@(x) bootstrp(n_boot,@mean,x')',response_grid,'uni',false);

%% Get retinotopy (for each bootstrap)

use_method = 'com'; % max or com
vfs_boot = nan(size(Ud,1),size(Ud,2),n_boot);
for curr_boot = 1:n_boot
    
    response_mean = cell2mat(cellfun(@(x) x(:,curr_boot),response_mean_bootstrap(:),'uni',false)');
    stim_im_px = reshape(permute(svdFrameReconstruct(Ud(:,:,use_svs),response_mean(use_svs,:)),[3,1,2]),ny,nx,[]);
    gauss_filt = fspecial('gaussian',[ny,nx],filter_sigma);
    stim_im_smoothed = imfilter(imresize(stim_im_px,screen_resize_scale,'bilinear'),gauss_filt);
    
    switch use_method
        case 'max'
            % Upsample each pixel's response map and find maximum
            [~,mi] = max(reshape(stim_im_smoothed,[],size(stim_im_px,3)),[],1);
            [m_y,m_x] = ind2sub(size(stim_im_smoothed),mi);
            m_yr = reshape(m_y,size(Ud,1),size(Ud,2));
            m_xr = reshape(m_x,size(Ud,1),size(Ud,2));
            
        case 'com'
            % Conversely, do COM on original^2
            [xx,yy] = meshgrid(1:size(stim_im_smoothed,2),1:size(stim_im_smoothed,1));
            m_xr = reshape(sum(sum(bsxfun(@times,stim_im_smoothed.^2,xx),1),2)./sum(sum(stim_im_smoothed.^2,1),2),size(Ud,1),size(Ud,2));
            m_yr = reshape(sum(sum(bsxfun(@times,stim_im_smoothed.^2,yy),1),2)./sum(sum(stim_im_smoothed.^2,1),2),size(Ud,1),size(Ud,2));
    end
    
    % Calculate and plot sign map (dot product between horz & vert gradient)
    
    % 1) get gradient direction
    [~,Vdir] = imgradient(imgaussfilt(m_yr,1));
    [~,Hdir] = imgradient(imgaussfilt(m_xr,1));
    
    % 3) get sin(difference in direction) if retinotopic, H/V should be
    % orthogonal, so the closer the orthogonal the better (and get sign)
    angle_diff = sind(Vdir-Hdir);
    angle_diff(isnan(angle_diff)) = 0;
    
    vfs_boot(:,:,curr_boot) = angle_diff;
end


%% Plot retinotopy (median across bootstraps)

vfs_median = imgaussfilt(nanmean(vfs_boot,3),2);

figure('Name',[animal ' ' day]);
ax1 = axes;
subplot(1,2,1,ax1);
imagesc(vfs_median);
caxis([-1,1]);
axes(ax1); axis image off;
colormap(colormap_BlueWhiteRed)

ax2 = axes;
ax3 = axes;
subplot(1,2,2,ax2);
subplot(1,2,2,ax3);
h1 = imagesc(ax2,avg_im(use_u_y,:));
colormap(ax2,gray);
caxis(ax2,[0 prctile(avg_im(:),95)]);
h2 = imagesc(ax3,vfs_median);
colormap(ax3,colormap_BlueWhiteRed);
caxis([-1,1]);
set(ax2,'Visible','off');
axes(ax2); axis image off;
set(ax3,'Visible','off');
axes(ax3); axis image off;
set(h2,'AlphaData',mat2gray(abs(vfs_median))*0.5);
colormap(ax2,gray);

drawnow;


%% Align retinotopy to master retinotopy
% Plot CCF areas and coordinates aligned to master retinotopy

% Load master VFS
master_vfs_fn = ['C:\Users\Julie\Dropbox\MATLAB\Lilrig\lilrig_master_vfs.mat'];
load(master_vfs_fn);
master_align = master_vfs;

% Align animal image to master image
ref_size = size(master_align);

[optimizer, metric] = imregconfig('monomodal');
optimizer = registration.optimizer.OnePlusOneEvolutionary();
optimizer.MaximumIterations = 200;
optimizer.GrowthFactor = 1+1e-6;
optimizer.InitialRadius = 1e-4;

tformEstimate_affine = imregtform(vfs_median,master_align,'affine',optimizer,metric);
vfs_aligned = imwarp(vfs_median,tformEstimate_affine,'Outputview',imref2d(ref_size));
avg_im_aligned = imwarp(avg_im,tformEstimate_affine,'Outputview',imref2d(ref_size));

% Plot alignment
figure;
subplot(1,2,1);
imshowpair(master_align,vfs_median);
title('Unaligned VFS');
subplot(1,2,2);
imshowpair(master_align,vfs_aligned);
title('Aligned VFS');

% Plot average image with CCF overlay
figure;
overlay_color = 'b';

ax1 = axes;
ax2 = axes;
subplot(1,2,1,ax1);
subplot(1,2,1,ax2);
h1 = imagesc(ax1,avg_im_aligned);
colormap(ax1,gray);
caxis(ax1,[0 prctile(avg_im(:),98)]);
h2 = imagesc(ax2,vfs_aligned);
colormap(ax2,colormap_BlueWhiteRed);
caxis([-1,1]);
set(ax1,'Visible','off');
axes(ax1); axis image off;
set(ax2,'Visible','off');
axes(ax2); axis image off;
set(h2,'AlphaData',mat2gray(abs(vfs_aligned))*0.4);
colormap(ax1,gray);
title(ax1,'VFS/CCF overlay');
% (load and plot aligned CCF boundaries)
hold on;
load(['C:\Users\Julie\Dropbox\MATLAB\Lilrig\lilrig_cortical_area_boundaries_aligned.mat']);
h = cellfun(@(areas) cellfun(@(outline) plot(outline(:,2),outline(:,1),'color',overlay_color),areas,'uni',false), ...
    cortical_area_boundaries_aligned,'uni',false);

% Plot average image with grid overlay
subplot(1,2,2);
imagesc(avg_im_aligned);
caxis([0 prctile(avg_im(:),98)]);
colormap(gca,'gray');
axis image off
title('0.5mm grid overlay');
% (plot grid)
hold on
um2pixel = 20.6;
bregma = [540,0,570] + 0.5;
ccf_tform_fn = ['C:\Users\Julie\Dropbox\MATLAB\Lilrig\lilrig_ccf_tform.mat'];
load(ccf_tform_fn);
bregma_resize = bregma*(10/um2pixel);
bregma_align = [bregma_resize([3,1]),1]*ccf_tform.T;
bregma_offset_x = bregma_align(1);
bregma_offset_y = bregma_align(2);

spacing_um = 500;
spacing_pixels = spacing_um/um2pixel;

xlines_pos = bregma_offset_y + spacing_pixels*(ceil((min(ylim)-bregma_offset_y)./spacing_pixels):floor((max(ylim)-bregma_offset_y)./spacing_pixels));
ylines_pos = bregma_offset_x + spacing_pixels*(ceil((min(xlim)-bregma_offset_x)./spacing_pixels):floor((max(xlim)-bregma_offset_x)./spacing_pixels));

h = struct;

for curr_xline = 1:length(xlines_pos)
    h.xlines(curr_xline) = line(xlim,repmat(xlines_pos(curr_xline),1,2),'color',overlay_color,'linestyle','-');
end

for curr_yline = 1:length(ylines_pos)
    h.ylines(curr_yline) = line(repmat(ylines_pos(curr_yline),1,2),ylim,'color',overlay_color,'linestyle','-');
end

h.bregma = plot(bregma_offset_x,bregma_offset_y,'.r','MarkerSize',30);











