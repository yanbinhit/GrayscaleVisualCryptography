function ItoFunTest()
close all; clear;

% Use haltone test image
imageIn = imread('../Images/lena.tiff','tiff');
if size(imageIn, 3)>1
    imageIn = rgb2gray(imageIn);
end
imageHT = HalftoningED(imageIn);

% Use binary secret iamge
imageIn= imread('../Images/geometry.bmp','bmp'); % Use binary image as secret
imageHT = imageIn;

% Ito
type = 1;
[shares, stacked] = ItoFun(imageHT);

figure; imshow(shares(:,:,1),[]);
figure; imshow(shares(:,:,2),[]);
figure; imshow(stacked,[]);



% 
% % Compare with "Random grid, algorithm 1"
% type = 1;
% [shares, stacked] = RandomGridKafri(imageHT, type);
% figure;
% subplot(1,3,1); imshow(shares(:,:,1),[]);
% subplot(1,3,2); imshow(shares(:,:,2),[]);
% subplot(1,3,3); imshow(stacked,[]);
% 
% figure; imshow(stacked,[]); title('Random grid algorithm 1');