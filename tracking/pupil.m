% eyecam_t 

obj = VideoReader([eyecam_dir, '/eye.mj2']);
NumberOfFrames = obj.NumFrames;
% 
img=read(obj,10);
files=size(img,1);
columns=size(img,2);
% Center
centerFile=round(files/2);
centerCol=round(columns/2);
% crop image
 figure();
[Icropped, rectout] = imcrop(read(obj,1)); % rectout: [xmin ymin width height].
imagesc(read(obj,1))
 figure();
startFrame = find(~isnan(eyecam_t), 1, 'first');
stopFrame = size(eyecam_t,1) - find(~isnan(eyecam_t(end:-1:1)), 1, 'first') + 1;
nFramesToUse = stopFrame - startFrame;
for iFrame = startFrame:stopFrame

    
    img=read(obj,iFrame);
    img = img(round(rectout(2)):round(rectout(2)+rectout(4)),round(rectout(1)):round(rectout(1)+rectout(3)));
    if size(img,3)==3
        img=rgb2gray(img);
    end
    %subplot(212)
    %imagesc(img);
    skin=~im2bw(img);
    %imagesc(skin)
    %     --
    skin=bwmorph(skin,'close');
    %imagesc(skin)
    skin=bwmorph(skin,'open');
    %imagesc(skin)
    skin=bwareaopen(skin,200);
     imagesc(skin)
    skin=imfill(skin,'holes');
    imagesc(skin);
    % Tagged objects in BW image
    L=bwlabel(skin);
    % Get areas and tracking rectangle
    out_a=regionprops(L);
    % Count the number of objects
    N=size(out_a,1);
    if N < 1 || isempty(out_a) % Returns if no object in the image
        solo_cara=[ ];
        continue
    end
    % ---
    % Select larger area
    areas=[out_a.Area];
    [area_max pam]=max(areas);
    %subplot(211)
    %imagesc(img);
    %colormap gray
    %hold on
    %rectangle('Position',out_a(pam).BoundingBox,'EdgeColor',[1 0 0],...
    %    'Curvature', [1,1],'LineWidth',2)
    %center=round(out_a(pam).Centroid);
    %X=center(1);
    %Y=center(2);
    %plot(X,Y,'g+')
    %     
    %text(X+10,Y,['(',num2str(X),',',num2str(Y),')'],'Color',[1 1 1])
%     if X<centerCol && Y<centerFile
%         title('Top left')
%     elseif X>centerCol && Y<centerFile
%         title('Top right')
%     elseif X<centerCol && Y>centerFile
%         title('Bottom left')
%     else
%         title('Bottom right')
%     end
%     hold off
%     % --
%     drawnow;

    pupilArea_pixels(iFrame,:) = pi * out_a(pam).BoundingBox(3)/2 * out_a(pam).BoundingBox(4)/2;
    if floor(iFrame/(startFrame+1000)) == ceil(iFrame/(startFrame+1000))
        disp(['Got pupil area for ', num2str(floor(iFrame/(startFrame+1000))*1000),  '/', num2str(nFramesToUse), ' frames'])
    end
end



%
thisUnit = 92;
thisUnitSpikeTimes = spike_times_timeline(spike_templates==thisUnit);
binSize = 1;
timeBins = min(spike_times_timeline):binSize:max(spike_times_timeline);
[n, x] = hist(thisUnitSpikeTimes, timeBins);
n = n ./ binSize;


figure();
subplot(411)
title('unit 92'); hold on;
plot(x, smoothdata(n, 'movmean', [1, 5]))
xlabel('time (s)')
ylabel('inst. FR')
makepretty;
xlim([eyecam_t(startFrame), eyecam_t(stopFrame)])


thisUnit = 86;
thisUnitSpikeTimes = spike_times_timeline(spike_templates==thisUnit);
binSize = 1;
timeBins = min(spike_times_timeline):binSize:max(spike_times_timeline);
[n, x] = hist(thisUnitSpikeTimes, timeBins);
n = n ./ binSize;
subplot(412)
title('unit 86'); hold on;
plot(x, smoothdata(n, 'movmean', [1, 5]))
xlabel('time (s)')
ylabel('inst. FR')
makepretty;
xlim([eyecam_t(startFrame), eyecam_t(stopFrame)])


subplot(413)
plot(Timeline.rawDAQTimestamps, smoothdata(abs(wheel_velocity), 'movmean', [1 10000]))
xlabel('time (s)')
ylabel('wheel velocity')
makepretty;
xlim([eyecam_t(startFrame), eyecam_t(stopFrame)])


subplot(414)
plot(eyecam_t(startFrame:stopFrame),smoothdata(pupilArea_pixels(startFrame:end), 'movmedian', [1,100]))
xlabel('time (s)')
ylabel('pupil area (a.u.)')
makepretty;

figure();
iC = 0;
for iTime = [880, 987,1040,1072,1096,1119 1618]
thisFrame = find(eyecam_t>iTime, 1, 'first');
iC=iC+1;
subplot(3,3,iC)
imagesc(read(obj,thisFrame)); colormap(gray)
title(['time = ' num2str(iTime)]); 
end

objFace = VideoReader([eyecam_dir, '/face.mj2']);
NumberOfFramesFace  = objFace .NumFrames;
% 
imgFace =read(objFace ,10);
filesFace =size(imgFace ,1);
columnsFace =size(imgFace ,2);

figure();
iC = 0;
for iTime = [880, 987,1040,1072,1096,1119 1618]
thisFrame = find(facecam_t>iTime, 1, 'first');
iC=iC+1;
subplot(3,3,iC)
imagesc(read(objFace,thisFrame)); colormap(gray)
title(['time = ' num2str(iTime)]); 
end

% 
thisUnit = 92;
thisUnitSpikeTimes = spike_times_timeline(spike_templates==thisUnit);
%[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times_timeline(spike_templates==thisUnit), eventTimes, thisWindow, psthBinSize);



binSize = 1;
timeBins = min(spike_times_timeline):binSize:max(spike_times_timeline);
[n, x] = hist(thisUnitSpikeTimes, timeBins);
n = n ./ binSize;


figure();
subplot(511)
plot(x, smoothdata(n, 'movmean', [1, 10]))
xlabel('time (s)')
ylabel('inst. FR')
makepretty;
xlim([eyecam_t(startFrame), eyecam_t(stopFrame)])

thisUnit = 82;
thisUnitSpikeTimes = spike_times_timeline(spike_templates==thisUnit);
binSize = 1;
timeBins = min(spike_times_timeline):binSize:max(spike_times_timeline);
[n, x] = hist(thisUnitSpikeTimes, timeBins);
n = n ./ binSize;
subplot(512)
plot(x, smoothdata(n, 'movmean', [1, 10]))
xlabel('time (s)')
ylabel('inst. FR')
makepretty;
xlim([eyecam_t(startFrame), eyecam_t(stopFrame)])

thisUnit = 87;
thisUnitSpikeTimes = spike_times_timeline(spike_templates==thisUnit);
binSize = 1;
timeBins = min(spike_times_timeline):binSize:max(spike_times_timeline);
[n, x] = hist(thisUnitSpikeTimes, timeBins);
n = n ./ binSize;
subplot(513)
plot(x, smoothdata(n, 'movmean', [1, 10]))
xlabel('time (s)')
ylabel('inst. FR')
makepretty;
xlim([eyecam_t(startFrame), eyecam_t(stopFrame)])

thisUnit = 79;
thisUnitSpikeTimes = spike_times_timeline(spike_templates==thisUnit);
binSize = 1;
timeBins = min(spike_times_timeline):binSize:max(spike_times_timeline);
[n, x] = hist(thisUnitSpikeTimes, timeBins);
n = n ./ binSize;
subplot(513)
plot(x, smoothdata(n, 'movmean', [1, 10]))
xlabel('time (s)')
ylabel('inst. FR')
makepretty;
xlim([eyecam_t(startFrame), eyecam_t(stopFrame)])

thisUnit = 86;
thisUnitSpikeTimes = spike_times_timeline(spike_templates==thisUnit);
binSize = 1;
timeBins = min(spike_times_timeline):binSize:max(spike_times_timeline);
[n, x] = hist(thisUnitSpikeTimes, timeBins);
n = n ./ binSize;
subplot(514)
plot(x, smoothdata(n, 'movmean', [1, 10]))
xlabel('time (s)')
ylabel('inst. FR')
makepretty;
xlim([eyecam_t(startFrame), eyecam_t(stopFrame)])

% facemap 
face_proc = load('/home/netshare/zinu/JF070/2022-06-10/1/face_proc.mat');
eye_proc = load('/home/netshare/zinu/JF070/2022-06-10/1/eye_proc.mat');
figure();
imagesc(face_proc.avgmotion_0)
colormap(brewermap([],'*RdBu'))

figure(); 
plot(facecam_t, smoothdata(nanmean(face_proc.motSVD_0,2), 'movmean', [1, 1000]))
makepretty;
xlim([eyecam_t(startFrame), eyecam_t(stopFrame)])