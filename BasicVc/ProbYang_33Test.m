function ProbYang_33Test()
close all; clear;
%imageIn = imread('../images/lena.tiff','tiff');
% Use binary secret iamge
imageIn= imread('../Images/geometry.bmp','bmp'); % Use binary image as secret
if size(imageIn, 3)>1
    imageIn = rgb2gray(imageIn);
end
[shares, stacked, imageHT] = ProbYang_33(imageIn);
figure; 
for i=1:3
    figure; imshow(shares(:,:,i),[]);
end

% stack 1,2 
figure;
imshow(~(~shares(:,:,1)|~shares(:,:,2)), []);
title('Stacking 1 and 2');

% stack 1,3
figure;
imshow(~(~shares(:,:,1)|~shares(:,:,3)), []);
title('Stacking 1 and 3');

figure;
imshow(stacked,[]);