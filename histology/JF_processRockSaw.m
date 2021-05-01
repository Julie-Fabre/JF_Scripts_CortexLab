%% rocksaw histology process 
animal = 'JF030'; 
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
figure();imagesc(tv_coronal(:,:,end-10));colormap(gray)
tv_coronal_20um = tv_coronal(1:2:end,1:2:end,1:2:end);

%% fiji step: crop images: remove cerebellum form template and image, and olf. bulbs from image 

%% automatic registering using elastix 
dd=dir('D:\JF030\Substack*');
thisTiff = loadtiff(['D:\JF030\' dd.name]);
params = {'01_ARA_affine.txt',...
    '02_ARA_bspline.txt'};
[reg,elastixLog] = elastix(squeeze(RGBimg(:,:,:,1)),thisTiff ,[],params);%green channel - don't use of teto mouse 

%get red channel 
[regRed,elastixLog] = transformix(squeeze(RGBimg(:,:,:,1)),elastixLog);%red channel 
saveastiff(reg, ['D:/' animal '/regRedChan.tiff']); 

%% fiji step: save template + image overlaid, keep only relevant bits and downsample in z 

%% histology 

AP_get_probe_histologyJF(tv, av, st, slice_path);

%[reg,log] = transformix(squeeze(RGBimg(:,:,:,2)),elastixLog);

%save as multipage tiff 
saveastiff(reg, ['D:/' animal '/reg2.tiff']); 

figure();
subplot(131)
imagesc(squeeze(squeeze(RGBimg(:,:,200,1))))
subplot(132)
imagesc(squeeze(reg(:,:,200)))




%% user refinement stage : alignement

%% draw probes 

%% locate probe depth woth ephys

%% save 