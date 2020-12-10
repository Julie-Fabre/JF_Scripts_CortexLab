%make stims 

stim = zeros(375,1125)+128; %grey

wDeg = 1;  %size of image (in degrees)
nPix = 1125;  %resolution of image (pixels);

[x,y] = meshgrid(linspace(-wDeg/2,wDeg/2,nPix+1));
x = x(1:end-1,1:end-1);
y = y(1:end-1,1:end-1);

orientation = 90;  %deg (counter-clockwise from horizontal)
sf = 12; %spatial frequency (cycles/deg)

ramp = sin(orientation*pi/180)*x-cos(orientation*pi/180)*y;

grating = uint8(255 * mat2gray(sin(2*pi*sf*ramp)))+1;
imagesc(grating);

img = stim;
img(1:125,:)=grating(1:125,:);
imagesc(img)
colormap(gray)
save('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\locationsGratings\img1.mat', 'img')

img = stim;
img(125:250,:)=grating(125:250,:);
imagesc(img)
save('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\locationsGratings\img2.mat', 'img')

img = stim;
img(250:end,:)=grating(250:375,:);
imagesc(img)
save('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\locationsGratings\img3.mat', 'img')

img = stim;
img(:,1:375)=grating(1:375,1:375);
imagesc(img)
colormap(gray)
save('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\locationsGratings\img4.mat', 'img')

img = stim;
img(:,375:750)=grating(1:375,375:750);
imagesc(img)
save('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\locationsGratings\img5.mat', 'img')

img = stim;
img(:,750:end)=grating(1:375,750:end);
imagesc(img)
save('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\locationsGratings\img6.mat', 'img')


