%make stims 

stim = zeros(375,1125);

img = stim;
img(1:125,:)=255;
imagesc(img)
colormap(gray)
save('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\locations\img1.mat', 'img')

img = stim;
img(125:250,:)=1;
imagesc(img)
save('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\locations\img2.mat', 'img')

img = stim;
img(250:end,:)=1;
imagesc(img)
save('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\locations\img3.mat', 'img')

img = stim;
img(:,1:375)=1;
imagesc(img)
save('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\locations\img4.mat', 'img')

img = stim;
img(:,375:375+375)=1;
imagesc(img)
save('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\locations\img5.mat', 'img')

img = stim;
img(:,1125-375:end)=1;
imagesc(img)
save('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\locations\img6.mat', 'img')