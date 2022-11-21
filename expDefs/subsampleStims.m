%% subsample, downsample

for iImg = 1:30
    try
img1 = imread(['/home/netshare/zserver/pregenerated_textures/JulieF/naturalImages/img' num2str(iImg) '.jpeg']);
% figure(1);
% imagesc(img1)
% colormap(gray)
% 
% figure(2);
% 
% imagesc(imresize(img1, [140, 160]))
% colormap(gray)
try
im = imresize(rgb2gray(img1), [140, 160]);
imwrite(im,gray,['/home/julie/Documents/new/img' num2str(iImg) '.jpeg'])

catch
    img1 = img1(10:end-10, 10:end-10);
    im = imresize(img1, [140, 160]);
    imwrite(im,gray,['/home/julie/Documents/new/img' num2str(iImg) '.jpeg'])

end
    catch
    end
end

%% .mats 
for iImg = 1:6
    
img1 = load(['/home/netshare/zserver-data/pregenerated_textures/JulieF/locationsGratings/img' num2str(iImg) '.mat']);
% figure(1);
% imagesc(img1.img)
% colormap(gray)
% 
% figure(2);
% 
% imagesc(imresize(img1.img, [70, 80]))
% colormap(gray)

img1.img = imresize((img1.img1), [70, 80]);
save(['/home/julie/Documents/new/img' num2str(iImg) '.mat'],'img1')


end

%% .mats 
for iImg = 1:6
   % iImg= iImg +1
load(['/home/netshare/zserver-data/pregenerated_textures/JulieF/locationsGratings/img' num2str(iImg) '.mat']);
% figure(1);
% imagesc(img1.img)
% colormap(gray)
% 
% figure(2);
% 
% imagesc(imresize(img1.img, [70, 80]))
% colormap(gray)
%img= img1.img1.img;
%img1= rmfield(img1, 'infomg')
img = imresize((img), [70, 80]);
save(['/home/julie/Documents/new/img' num2str(iImg) '.mat'],'img')


end

%% convert 2 to center only 
iImg = 26
img1 = imread(['/home/netshare/zserver/pregenerated_textures/JulieF/noGoWorld_Passive/img' num2str(iImg) '.jpeg']);
imagesc(img1)
figure();
im = img1;
im([1:53, 94:end],:,:) = 128;

imagesc(im)
imwrite(rgb2gray(im),gray,['/home/julie/Documents/new/img' num2str(iImg) '.jpeg'])

%% add black square 

for iImg = 1:30
  % iImg= iImg +1
load(['/home/netshare/zserver-data/pregenerated_textures/JulieF/backup/naturalImages/img' num2str(iImg) '.mat']);
% figure(1);
% imagesc(img)
% sz = size(img);
% xg = 1:sz(1);
% yg = 1:sz(2);
% F = griddedInterpolant({xg,yg},double(img));
% 
% xq = (0:640/500:sz(1))';
% yq = (0:512/375:sz(2))';
% vq = uint8(F({xq,yq}));
% imagesc(vq)
% colormap(gray)
img = img(124:375-124, 11:500-10);

%img(:,501: 500+13) = img(:,500-12:500);
sz = size(img);
sc(1) = sz(1)/1024;
sc(2) = sz(2)/3840;
pxSt = sc * 100;
%img(sz(1)- ceil(pxSt(1)*2): sz(1)- floor(pxSt(1)), :) = -1.85;
img(sz(1)- ceil(pxSt(1)*3): sz(1)- floor(pxSt(1))*2,sz(2)- ceil(pxSt(2)*3) : sz(2)- floor(pxSt(2)*2)) = -1.85;
%img(:, sz(2)- ceil(pxSt(2)*3) : sz(2)- floor(pxSt(2)*2)) = -1.85;
figure(2)
imagesc(img)
%disp(nanmean(nanmean(img(sz(1)- ceil(pxSt(1)*2): sz(1)- floor(pxSt(1)), sz(2)- ceil(pxSt(2)*2) : sz(2)- floor(pxSt(2))))))

save(['/home/julie/Documents/new/natIm/img' num2str(iImg) '.mat'],'img')


end


for iImg = 1:6
  % iImg= iImg +1
load(['/home/netshare/zserver-data/pregenerated_textures/JulieF/backup/locationsGratings/img' num2str(iImg) '.mat']);
% figure(1);
% imagesc(img)
% colormap(gray)
% 
% figure(2);
% 
% imagesc(imresize(img1.img, [70, 80]))
% colormap(gray)
%img= img1.img1.img;
%img1= rmfield(img1, 'infomg')
%img = imresize((img), [70, 80]);
sz = size(img);
sc(1) = sz(1)/1024;
sc(2) = sz(2)/3840;
pxSt = sc * 100;
img(sz(1)- ceil(pxSt(1)*3): sz(1)- floor(pxSt(1))*2,sz(2)- ceil(pxSt(2)*3) : sz(2)- floor(pxSt(2)*2)) = -1.85;
imagesc(img)
save(['/home/julie/Documents/new/locG/img' num2str(iImg) '.mat'],'img')
end

for iImg = 1:3

  % iImg= iImg +1
  img1 = imread(['/home/netshare/zserver/pregenerated_textures/JulieF/backup/noGoWorld_Passive_onlyTask/img' num2str(iImg) '.jpeg']);
% figure(1);
% imagesc(img)
% colormap(gray)

if size(img1,3) ==3 
    img1_gray = 127 * ones(size(img1,1), size(img1,2), size(img1,3));
    img1 = [img1_gray, img1, img1_gray];
    sz = size(img1);
sc(1) = sz(1)/1024;
sc(2) = sz(2)/3840;
pxSt = sc * 100;
    width = sz(1)- ceil(pxSt(1)*2): sz(1)- floor(pxSt(1));
    height = sz(2)- ceil(pxSt(2)*2) : sz(2)- floor(pxSt(2));
    img1(sz(1)- ceil(pxSt(1)*2): sz(1)- floor(pxSt(1)),sz(2)- ceil(pxSt(2)*2) : sz(2)- floor(pxSt(2)),:) = repmat(0,...
        [max(width) - min(width) + 1, ...
        max(height) - min(height) + 1, 3]);
    img1= rgb2gray(img1);

else
    img1_gray = 127 * ones(size(img1,1), size(img1,2), size(img1,3));
    img1 = [img1_gray, img1, img1_gray];
    sz = size(img1);
sc(1) = sz(1)/1024;
sc(2) = sz(2)/3840;
pxSt = sc * 100;

img1(sz(1)- ceil(pxSt(1)*3): sz(1)- floor(pxSt(1))*2,sz(2)- ceil(pxSt(2)*3) : sz(2)- floor(pxSt(2)*2)) = 0;
end
%disp(nanmean(img1(sz(1)- ceil(pxSt(1)*2): sz(1)- floor(pxSt(1)), sz(2)- ceil(pxSt(2)*2) : sz(2)- floor(pxSt(2)))))
imagesc(img1)
imwrite(img1,gray,['/home/julie/Documents/new/cw/img' num2str(iImg) '.jpeg'])

end

for iImg = 4:6

  % iImg= iImg +1
  img1 = imread(['/home/netshare/zserver/pregenerated_textures/JulieF/backup/noGoWorld_Passive_onlyTask/img' num2str(iImg) '.jpeg']);
% figure(1);
% imagesc(img)
% colormap(gray)

if size(img1,3) ==3 
    img1_gray = 127 * ones(size(img1,1), size(img1,2), size(img1,3));
    img1 = [img1, img1_gray,img1_gray];
    sz = size(img1);
sc(1) = sz(1)/1024;
sc(2) = sz(2)/3840;
pxSt = sc * 100;
    width = sz(1)- ceil(pxSt(1)*2): sz(1)- floor(pxSt(1));
    height = sz(2)- ceil(pxSt(2)*2) : sz(2)- floor(pxSt(2));
    img1(sz(1)- ceil(pxSt(1)*2): sz(1)- floor(pxSt(1))*1,sz(2)- ceil(pxSt(2)*2) : sz(2)- floor(pxSt(2)*1),:) = repmat(0,...
        [max(width) - min(width) + 1, ...
        max(height) - min(height) + 1, 3]);
    img1= rgb2gray(img1);

else
    img1_gray = 127 * ones(size(img1,1), size(img1,2), size(img1,3));
    img1 = [img1_gray, img1, img1_gray];
    sz = size(img1);
sc(1) = sz(1)/1024;
sc(2) = sz(2)/3840;
pxSt = sc * 100;

img1(sz(1)- ceil(pxSt(1)*3): sz(1)- floor(pxSt(1))*2,sz(2)- ceil(pxSt(2)*3) : sz(2)- floor(pxSt(2)*2)) = 0;
end
%disp(nanmean(img1(sz(1)- ceil(pxSt(1)*2): sz(1)- floor(pxSt(1)), sz(2)- ceil(pxSt(2)*2) : sz(2)- floor(pxSt(2)))))
imagesc(img1)
imwrite(img1,gray,['/home/julie/Documents/new/cw/img' num2str(iImg) '.jpeg'])

end



for iImg = 1:13


  % iImg= iImg +1
  img1 = imread(['/home/netshare/zserver/pregenerated_textures/JulieF/noGoWorld_Passive/img' num2str(iImg) '.jpeg']);
% figure(1);
% imagesc(img)
% colormap(gray)

if size(img1,3) ==3 
    img1_gray = 127 * ones(size(img1,1), size(img1,2), size(img1,3));
    img1 = [img1_gray, img1, img1_gray];
    sz = size(img1);
sc(1) = sz(1)/1024;
sc(2) = sz(2)/3840;
pxSt = sc * 100;
    width = sz(1)- ceil(pxSt(1)*3): sz(1)- floor(pxSt(1))*2;
    height = sz(2)- ceil(pxSt(2)*3) : sz(2)- floor(pxSt(2)*2);
    img1(sz(1)- ceil(pxSt(1)*3): sz(1)- floor(pxSt(1))*2,sz(2)- ceil(pxSt(2)*3) : sz(2)- floor(pxSt(2)*2),:) = repmat(0,...
        [max(width) - min(width) + 1, ...
        max(height) - min(height) + 1, 3]);
    img1= rgb2gray(img1);

else
    img1_gray = 127 * ones(size(img1,1), size(img1,2), size(img1,3));
    img1 = [img1_gray, img1, img1_gray];
    sz = size(img1);
sc(1) = sz(1)/1024;
sc(2) = sz(2)/3840;
pxSt = sc * 100;

img1(sz(1)- ceil(pxSt(1)*3): sz(1)- floor(pxSt(1))*2,sz(2)- ceil(pxSt(2)*3) : sz(2)- floor(pxSt(2)*2)) = 0;
end
%disp(nanmean(img1(sz(1)- ceil(pxSt(1)*2): sz(1)- floor(pxSt(1)), sz(2)- ceil(pxSt(2)*2) : sz(2)- floor(pxSt(2)))))
imagesc(img1)
imwrite(img1,gray,['/home/julie/Documents/new/cw_all/img' num2str(iImg) '.jpeg'])

end



for iImg = 14:26

  % iImg= iImg +1
  img1 = imread(['/home/netshare/zserver/pregenerated_textures/JulieF/noGoWorld_Passive/img' num2str(iImg) '.jpeg']);
% figure(1);
% imagesc(img)
% colormap(gray)

if size(img1,3) ==3 
    img1_gray = 127 * ones(size(img1,1), size(img1,2), size(img1,3));
    img1 = [ img1, img1_gray, img1_gray];
    sz = size(img1);
sc(1) = sz(1)/1024;
sc(2) = sz(2)/3840;
pxSt = sc * 100;
    width = sz(1)- ceil(pxSt(1)*3): sz(1)- floor(pxSt(1))*2;
    height = sz(2)- ceil(pxSt(2)*3) : sz(2)- floor(pxSt(2)*2);
    img1(sz(1)- ceil(pxSt(1)*3): sz(1)- floor(pxSt(1))*2,sz(2)- ceil(pxSt(2)*3) : sz(2)- floor(pxSt(2)*2),:) = repmat(0,...
        [max(width) - min(width) + 1, ...
        max(height) - min(height) + 1, 3]);
    img1= rgb2gray(img1);

else
    img1_gray = 127 * ones(size(img1,1), size(img1,2), size(img1,3));
    img1 = [img1_gray, img1, img1_gray];
    sz = size(img1);
sc(1) = sz(1)/1024;
sc(2) = sz(2)/3840;
pxSt = sc * 100;

img1(sz(1)- ceil(pxSt(1)*3): sz(1)- floor(pxSt(1))*2,sz(2)- ceil(pxSt(2)*3) : sz(2)- floor(pxSt(2)*2)) = 0;
end
%disp(nanmean(img1(sz(1)- ceil(pxSt(1)*2): sz(1)- floor(pxSt(1)), sz(2)- ceil(pxSt(2)*2) : sz(2)- floor(pxSt(2)))))
imagesc(img1)
imwrite(img1,gray,['/home/julie/Documents/new/cw_all/img' num2str(iImg) '.jpeg'])

end
