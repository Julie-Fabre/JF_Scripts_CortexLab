function motionEnergy = JF_getMotionEnergy(facecam_dir, facecam_t)

vidObj = VideoReader(facecam_dir);
startFrame = find(~isnan(facecam_t), 1, 'first');
stopFrame = size(facecam_t,1) - find(~isnan(facecam_t(end:-1:1)), 1, 'first') + 1;

for iFrame = startFrame:stopFrame-1

    figure(); image(img)
    img=read(vidObj,iFrame);
    img2=read(vidObj,iFrame+1);
    motionEnergy(iFrame) = abs(sum(sum(img))) - abs(sum(sum(img2)));
% whisker movement was defined as the absolute value of the difference between consecutive frames, summed across pixels
end