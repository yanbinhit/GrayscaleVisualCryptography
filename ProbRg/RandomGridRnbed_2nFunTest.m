function RandomGridRnbed_2nFunTest()
close all; clear;
prevDir = pwd;
[dir, dummy, dummy2] = fileparts(mfilename('fullpath'));
cd(dir);
n=2;

imageIn = imread('../Images/lena.tiff','tif');
[shares, stacked] = RandomGridRnbed_2nFun(imageIn, n);
sImg = HalftoningED(rgb2gray(imageIn));

figure; imshow(sImg,[]);
figure; imshow(shares(:,:,1),[]); title('share 1');
figure; imshow(shares(:,:,2),[]); title('share 2');
figure; imshow(stacked,[]);