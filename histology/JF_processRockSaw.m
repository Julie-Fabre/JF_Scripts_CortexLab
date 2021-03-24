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

%% allen tlas 25umcoronal 
% from http://help.brain-map.org/display/mousebrain/API
% 25 micron volume size
thisSize = [528 320 456];
% VOL = 3-D matrix of atlas Nissl volume
fid = fopen('C:\Users\Julie\Dropbox\Atlas\allenCCF\25um\P56_atlasVolume\atlasVolume\atlasVolume.raw', 'r', 'l' );
VOL = fread( fid, prod(thisSize), 'uint8' );
fclose( fid );
VOL = reshape(VOL,thisSize);

%% elastix parameters 
p.Transform='AffineTransform';
p.MaximumNumberOfIterations=400;
p.NumberOfSpatialSamples=300;
p.FixedImageDimension = 3;
p.MovingImageDimension = 3;
p.SP_alpha=0.2;

%% automatic registering using elastix 
reg=elastix(squeeze(RGBimg(:,:,:,1)),VOL,[],'elastix_default.yml','paramstruct',p);
%save as multipage tiff 

saveastiff(reg, ['D:/JF025/reg',num2str(p.SP_alpha),'.tiff']); 

figure();
subplot(131)
imagesc(squeeze(squeeze(RGBimg(:,:,200,1))))
subplot(132)
imagesc(squeeze(reg(:,:,200)))
%% user refinement stage 