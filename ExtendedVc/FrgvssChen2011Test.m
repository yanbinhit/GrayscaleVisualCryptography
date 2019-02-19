function FrgvssChen2011Test()
close all; clear all; 
nRowsCover = 512;
nColsCover = 512;

% cover images
cover = zeros(nRowsCover, nColsCover, 2); % 2 cover images
img1 = double(rgb2gray(imread('../Images/lena.tiff','tiff')));  % cover image 1
cover(:,:,1) = imresize(img1,[nRowsCover nColsCover]);
img2 = double(rgb2gray(imread('../Images/peppers.tiff','tiff'))); % cover image 2
cover(:,:,2) = imresize(img2,[nRowsCover nColsCover]);

% secret image
secretImgIn = imread('../Images/geometry.bmp','bmp');
secretImgIn = imresize(secretImgIn,[nRowsCover nColsCover],'nearest');
secretImgIn = secretImgIn > 0;

% FRRGVSS 
[shares, stacked] = ...
    FrgvssChen2011(cover, secretImgIn, nRowsCover, nColsCover);

figure;imshow(shares(:,:,1),[]);
figure;imshow(shares(:,:,2),[]);
figure; imshow(stacked,[]); 
%title('displayed on screen');

contrast = SecretContrast(secretImgIn, stacked)
[mssim, ssim_map(:,:)] = ssim(secretImgIn*255, stacked*255);
mssim
%----------------------------------------
% Compare quality 
%----------------------------------------
% Compare the PSNR of two filtered images
sigma = 2;
hSize = 11;
h = fspecial('gaussian', hSize, sigma);
for k = 1:2
    sharesFilt(:,:,k) = imfilter(shares(:,:,k), h, 'replicate');
end

psnr = zeros(1,2);
for k = 1:2
    psnr(k) = calPSNR(cover(:,:,k), sharesFilt(:,:,k));
end
psnr
psnrAvg = sum(psnr)/numel(psnr)

figure;subplot(1,2,1); imshow(sharesFilt(:,:,1),[]);
subplot(1,2,2); imshow(sharesFilt(:,:,2),[]);

% calculate SSIM 
for k = 1:2
 [mssim(k), ssim_map(:,:,k)] = ssim(cover(:,:,k), shares(:,:,k));
end
mssim


for k = 1:2
 [mssim(k), ssim_map(:,:,k)] = ssim(cover(:,:,k), sharesFilt(:,:,k));
end
mssim

 