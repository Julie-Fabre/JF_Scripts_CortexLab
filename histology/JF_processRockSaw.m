%% rocksaw histology process 
animal = 'JF025'; 
%% load images 
imPath = ['\\znas.cortexlab.net\Brainsaw\', animal, '\downsampled_stacks\025_micron\']; 
imgs = dir([imPath, '/*.tif']);
for iColor = 1:size(imgs,1)
    numimgs = size(imfinfo([imPath, imgs(iColor).name]),1);
    for iImage = 1:numimgs
        RGBimg(:, :, iImage, iColor) = imread([imPath, imgs(iColor).name],iImage); 
    end
end
%rgb to grayscale 

%% allen atlas 10um sagital 
allen_atlas_path = 'C:\Users\Julie\Dropbox\Atlas\allenCCF';
tv = readNPY([allen_atlas_path, filesep, 'template_volume_10um.npy']);
av = readNPY([allen_atlas_path, filesep, 'annotation_volume_10um_by_index.npy']);
st = loadStructureTreeJF([allen_atlas_path, filesep, 'structure_tree_safe_2017.csv']);

%sagital -> coronal 
tv_coronal = permute(tv, [2,3,1]);
figure();imagesc(tv_coronal(:,:,600));colormap(gray)
tv_coronal_20um = tv_coronal(1:2:end,1:2:end,1:10:end);

%% allen tlas 25um sagital 
% from http://help.brain-map.org/display/mousebrain/API
% 25 micron volume size
thisSize = [528 320 456];
% VOL = 3-D matrix of atlas Nissl volume
fid = fopen('C:\Users\Julie\Dropbox\Atlas\allenCCF\25um\P56_atlasVolume\atlasVolume\atlasVolume.raw', 'r', 'l' );
VOL = fread( fid, prod(thisSize), 'uint8' );
fclose( fid );
VOL = reshape(VOL,thisSize);

%sagital -> coronal 
VOL_coronal = permute(VOL, [2,3,1]);
figure();imagesc(VOL_coronal(:,:,200));colormap(gray)
%tv_coronal_20um = tv_coronal(:,1:2:end,1:10:end);%error if too much of a difference in size between data and ref

%% elastix parameters 
p.Transform='AffineTransform';
p.MaximumNumberOfIterations=1000;
p.NumberOfSpatialSamples=600;
p.FixedImageDimension = 3;
p.MovingImageDimension = 3;
p.SP_alpha=0.2;
p.SP_a=4000;

% p.Transform='BSplineTransform';
% p.MaximumNumberOfIterations=1E3;
% p.NumberOfSpatialSamples=1E3;
% p.SP_a=4000;
%% automatic registering using elastix 
params = {'parameters_Affine.txt',...
    'parameters_BSpline.txt'};
[reg,elastixLog] = elastix(squeeze(RGBimg(:,:,:,2)),tv_coronal_20um ,[],params);

p.Transform='BSplineTransform';
p.MaximumNumberOfIterations=1000;
p.NumberOfSpatialSamples=600;
p.FixedImageDimension = 3;
p.MovingImageDimension = 3;
p.SP_alpha=0.2;
p.SP_a=4000;

% p.Transform='BSplineTransform';
% p.MaximumNumberOfIterations=1E3;
% p.NumberOfSpatialSamples=1E3;
% p.SP_a=4000;
%% automatic registering using elastix 
params = {'parameters_Affine.txt',...
    'parameters_BSpline.txt'};
[reg,elastixLog] = elastix(squeeze(RGBimg(:,:,:,2)),tv_coronal_20um ,[],params);

%[reg,log] = transformix(squeeze(RGBimg(:,:,:,2)),elastixLog);

%save as multipage tiff 
saveastiff(reg, ['D:/JF025/reg2g2',p.Transform, '-',num2str(p.SP_alpha), '-',num2str(p.NumberOfSpatialSamples),'-',num2str(p.MaximumNumberOfIterations),'.tiff']); 

figure();
subplot(131)
imagesc(squeeze(squeeze(RGBimg(:,:,200,1))))
subplot(132)
imagesc(squeeze(reg(:,:,200)))




%% user refinement stage 