%% REGISTER ROCKSAW PROCESSED BRAINS TO ALLEN ATLAS 
% dependancies: 
% - loadtiff 
% - matlabelastix 
% JF 2021/06/08
animal = 'JF025';
%% Open dialog box, select images to register + template and load 
[filename,filepath]=uigetfile('\\znas.cortexlab.net\Brainsaw\*.tif','Select autofluorescence channel imaged brain');
[filenameref,filepathref]=uigetfile({'C:\Users\Julie\Dropbox\Atlas\*.tif'},'Select reference');
greenChannel = loadtiff([filepath, filename]);
redChannel = loadtiff([filepath, filename(1:end-11), '1_red.tif']);
allenAtlas = loadtiff([filepathref, filenameref]);
allenAtlas10um = readNPY('C:\Users\Julie\Dropbox\Atlas\allenCCF\template_volume_10um.npy');
allenAtlas10um = flipud(rot90(permute(allenAtlas10um, [3,2,1])));
%% crop template to fit image to register 
StackSlider(greenChannel);
StackSlider(allenAtlas10um);
cropAllenLimits = [209, 1245]; 

allenAtlas10um = allenAtlas10um(:,:,cropAllenLimits(1):cropAllenLimits(2));

%% register, get transform, apply
dd=dir(['D:\' animal '\Substack*']);
cd('C:\Users\Julie\Dropbox\MATLAB\JF_scripts_CortexLab\histology')
params = {'01_ARA_affine.txt',...
    '02_ARA_bspline.txt'};

[reg,elastixLog] = elastix(squeeze(greenChannel),allenAtlas10um(1:2:end, 1:2:end, 1:2:end),[],params);%green channel - don't use of teto mouse 
saveastiff(reg, ['D:/' animal '/regGreenChan.tiff']); 

%get red channel 
[regRed,elastixLog] = transformix(redChannel,elastixLog);%red channel 
saveastiff(regRed, ['D:/' animal '/regRedChan.tiff']); 

% get in 'AP' format: histology_ccf.mat with tv_slices, av_slices,
% plane_ap, plane_ml, plane_dv 
allen_atlas_path = 'C:\Users\Julie\Dropbox\Atlas\allenCCF';
tv = readNPY([allen_atlas_path, filesep, 'template_volume_10um.npy']);
tv = permute(tv, [2,1,3]);
av = readNPY([allen_atlas_path, filesep, 'annotation_volume_10um_by_index.npy']);
av = permute(av, [2,1,3]);
st = loadStructureTreeJF([allen_atlas_path, filesep, 'structure_tree_safe_2017.csv']);

histology_ccf=struct;
atlas2histology_tform = cell(size(reg,3),1);
for iSlice = 1:size(reg,3)
    histology_ccf(iSlice).plane_ap = repmat(iSlice + cropAllenLimits(1), [size(reg,1), size(reg,2)]);
    histology_ccf(iSlice).plane_ml = repmat(1:2:1140, [size(reg,1),1]);
    histology_ccf(iSlice).plane_dv = repmat(1:2:800, [ size(reg,2),1])';
    histology_ccf(iSlice).tv_slices = squeeze(tv(iSlice + cropAllenLimits(1),:,:));
    histology_ccf(iSlice).av_slices = squeeze(av(iSlice + cropAllenLimits(1),:,:));
    atlas2histology_tform{iSlice} = [1 0 0; 0 1 0; 0 0 1]; % no scaling 
end

save(['D:/' animal '/slices/histology_ccf.mat'], 'histology_ccf','-v7.3')
save(['D:/' animal '/slices/atlas2histology_tform.mat'], 'atlas2histology_tform')

%% draw probes 
AP_get_probe_histologyJF(tv, av, st, ['D:/',animal,'/slices'],'rocksaw',regRed);
%% align ephys depth