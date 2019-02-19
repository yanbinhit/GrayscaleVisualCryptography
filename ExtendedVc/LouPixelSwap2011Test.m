function LouPixelSwap2011Test()
close all; clear all; 

nRowsCover = 768;
nColsCover = 768;

% cover images
cover = zeros(nRowsCover, nColsCover, 2); % 2 cover images
img1 = double(rgb2gray(imread('../Images/lena.tiff','tiff')));  % cover image 1
%img1 = double((imread('baboon512.bmp','bmp')));  % cover image 1
cover(:,:,1) = imresize(img1,[nRowsCover nColsCover]);
img2 = double(rgb2gray(imread('../Images/peppers.tiff','tiff'))); % cover image 2
%img2 = double((imread('../Images/peppers.tiff','tiff'))); % cover image 2
cover(:,:,2) = imresize(img2,[nRowsCover nColsCover]);

figure; imshow(cover(:,:,1),[]); 
figure; imshow(cover(:,:,2),[]);

% cover(:,:,:) = cover(:,:,:)-50; % 实验增加、降低平均亮度是否会改善叠加质量
% for i = 1:nRowsCover
%     for j = 1:nColsCover
%         for k = 1:2
%             if cover(i,j,k)<0
%                 cover(i,j,k) = 0;
%             elseif cover(i,j,k)>=255
%                 cover(i,j,k) = 255;
%             end
%         end
%     end
% end

% secret image
secretImgIn = imread('../Images/geometry.bmp','bmp');
secretImgIn = imresize(secretImgIn,[nRowsCover nColsCover],'nearest');
secretImgIn = secretImgIn > 0;



% generate shares and simulate the stacking operation
[shares, stacked] = ...
    LouPixelSwap2011(cover, secretImgIn, nRowsCover, nColsCover);

figure;imshow(shares(:,:,1),[]);
figure;imshow(shares(:,:,2),[]);
figure; imshow(stacked,[]); 
%title('displayed on screen');