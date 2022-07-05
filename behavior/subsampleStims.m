%% subsample

for iImg = 1:42
    try
img1 = imread(['/home/netshare/zserver/pregenerated_textures/JulieF/choiceWorld_7Stims_Passive/img' num2str(iImg) '.jpeg']);
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

%% convert 2 to center only 
iImg = 26
img1 = imread(['/home/netshare/zserver/pregenerated_textures/JulieF/noGoWorld_Passive/img' num2str(iImg) '.jpeg']);
imagesc(img1)
figure();
im = img1;
im([1:53, 94:end],:,:) = 128;

imagesc(im)
imwrite(rgb2gray(im),gray,['/home/julie/Documents/new/img' num2str(iImg) '.jpeg'])