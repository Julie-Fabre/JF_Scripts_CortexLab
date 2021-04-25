%% plot histology 

animal = {'AP080'};
probeColors = {rgb('DeepPink');rgb('Yellow');rgb('Green');rgb('Red');rgb('OrangeRed');...
    rgb('Blue'); rgb('Lime'); rgb('Aqua'); rgb('LightSeaGreen'); rgb('Purple');...
    rgb('DarkKhaki');rgb('OliveDrab');rgb('DarkRed');rgb('Orchid');rgb('RosyBrown')};%15 colors


%% plot slices overlaid with probe 
%side by side: histology with probe one side, allen fit other side, ephys
%(average resp in 50um depth bins) 
load(['\\znas.cortexlab.net\Subjects\',animal{:},'\Histology\processed\slices\'])
load(['\\znas.cortexlab.net\Subjects\',animal{:},'\Histology\processed\probe2ephys.mat'])
load('\\znas.cortexlab.net\Subjects\AP080\Histology\processed\slices\histology_ccf.mat')
load('\\znas.cortexlab.net\Subjects\AP080\Histology\processed\slices\atlas2histology_tform.mat')
load('\\znas.cortexlab.net\Subjects\AP080\Histology\processed\slices\probe_ccf.mat')


imgs = dir(['\\znas.cortexlab.net\Subjects\',animal{:},'\Histology\processed\slices\*.tif']);

for iImg = 1:length(imgs)
    %load image
    fullHistology(iImg, :,:,:)=imread_rgb([imgs(1).folder, filesep, imgs(iImg).name]);
    %draw probe(s) on image 
    %probes = size(slice_slide_locations{1,iImg}{:},2);
    curr_av_slice = histology_ccf(iImg).av_slices;
    curr_av_slice(isnan(curr_av_slice)) = 1;
    iImg_im = fullHistology(iImg, :,:,:);
    
    tform = affine2d;
    h=imwarp(curr_av_slice,tform);   
    
    figure();
    subplot(321)
    set(gcf,'color','k');
    imagesc(squeeze(fullHistology(iImg, :,:,:)))
    axis square;
    
    subplot(322)
    set(gcf,'color','k');
    imagesc(squeeze(fullHistology(iImg, :,:,:)))
    hold on;
  
    
    allenBorders1...
        = [diff(curr_av_slice,[], 1) > 0; squeeze(zeros(1,size(curr_av_slice,2)))];
    allenBorders1(allenBorders1==0)=NaN;
    
    allenBorders2 = [diff(curr_av_slice,[], 2) > 0, squeeze(zeros(size(curr_av_slice,1),1))];
    allenBorders2 = double(allenBorders2);
    allenBorders2(allenBorders2==0)=NaN;
    
    allenBorders = cat(3,allenBorders1,allenBorders2); allenBorders = nansum(allenBorders,3);
    allenBordersAligned = imwarp(allenBorders,tform);
    allenBordersAligned(allenBordersAligned>0)=1;
    allenBordersAligned(allenBordersAligned==0)=NaN;
    
    imagesc(allenBordersAligned,'AlphaData',~isnan(allenBordersAligned))
   
    subplot(323)
    imagesc(h)
    
    subplot(313)
   
    
end

