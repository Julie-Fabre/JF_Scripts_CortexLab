%% make grating yourself 

%make stims 

stim = single(zeros(375,500)); %grey

wDeg = 1;  %size of image (in degrees)
nPix = 500;  %resolution of image (pixels);

[x,y] = meshgrid(linspace(-wDeg/2,wDeg/2,nPix+1));
x = x(1:end-1,1:end-1);
y = y(1:end-1,1:end-1);

orientation = 90;  %deg (counter-clockwise from horizontal)
sf = 5; %spatial frequency (cycles/deg)

ramp = sin(orientation*pi/180)*x-cos(orientation*pi/180)*y;

grating = single(4*mat2gray(sin(2*pi*sf*ramp)))+1;
imagesc(grating);

img = stim;
img(:,:)=grating(1:375,1:500)-2;
imagesc(img)
colormap(gray)
%img = single(img); 
save('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\choiceWorld_7Stims\img2.mat', 'img')

imwrite(img, '\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\choiceWorld_7Stims\img2.jpeg', 'JPEG')
img2=img;
img = stim;
img(1:125,:)=grating(1:125,:);
imagesc(img)
colormap(gray)
save('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\choiceWorld_7Stims\img3.mat', 'img')

img = stim;
img(125:250,:)=grating(125:250,:);
imagesc(img)
save('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\choiceWorld_7Stims\img4.mat', 'img')

img = stim;
img(250:end,:)=grating(250:375,:);
imagesc(img)
save('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\choiceWorld_7Stims\img5.mat', 'img')

%imread(\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\choiceWorld_7Stims\img6.jpeg')
sf = 15; %spatial frequency (cycles/deg)

ramp = sin(orientation*pi/180)*x-cos(orientation*pi/180)*y;

grating = single(4 * mat2gray(sin(2*pi*sf*ramp)))+1;
imagesc(grating);


img = stim;
img(:,:)=grating(1:375,1:500)-2;
imagesc(img)
colormap(gray)
save('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\choiceWorld_7Stims\img1.mat', 'img')
imwrite(img, '\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\choiceWorld_7Stims\img1.jpeg', 'JPEG')



img1= imread('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\choiceWorld_7Stims\img1.jpeg');
size(img1)
img = img1(:,:);
img(268:801,:) = 127; %gray background
imwrite(img, '\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\choiceWorld_7Stims\img3.jpeg', 'JPEG')

img = stim;
img(1:125,:)=grating(1:125,:);
imagesc(img)
colormap(gray)
save('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\choiceWorld_7Stims\img3.mat', 'img')

img = stim;
img(125:250,:)=grating(125:250,:);
imagesc(img)
save('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\choiceWorld_7Stims\img4.mat', 'img')

img = stim;
img(250:end,:)=grating(250:375,:);
imagesc(img)
save('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\choiceWorld_7Stims\img5.mat', 'img')



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

%% from screen shot (the lazy way) 
img1= imread('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\choiceWorld_7Stims\img1.jpeg'); %fullscreen grating 
size(img1)
img = img1(:,:);
img(1:267,:) = img1(268:534,:);
img(268:801,:) = 127; %gray background on middle and bottom
imwrite(img, '\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\choiceWorld_7Stims\img2.jpeg', 'JPEG')

img = img1(:,:);
img(534:800,:) = img1(268:534,:);
img(1:534,:) = 127; %gray background on top and middle
imwrite(img, '\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\choiceWorld_7Stims\img3.jpeg', 'JPEG')

img = img1(:,:);
img([1:268 534:801],:) = 127; %gray background on top and bottom
imwrite(img, '\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\choiceWorld_7Stims\img4.jpeg', 'JPEG')

img1= imread('\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\choiceWorld_7Stims\img1.jpeg'); %fullscreen grating 
size(img1)
img = img1(:,:);
img(1:267,:) = img1(268:534,:);
img(268:801,:) = 127; %gray background on middle and bottom
imwrite(img, '\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\choiceWorld_7Stims\img2.jpeg', 'JPEG')
