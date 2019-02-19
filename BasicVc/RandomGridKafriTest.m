function RandomGridKafriTest()
close all; clear;

imageIn= imread('../Images/geometry.bmp','bmp'); % for binary secret image
imageIn= imread('../Images/lena.bmp','bmp'); % using grayscale image

if size(imageIn, 3)>1
    imageIn = rgb2gray(imageIn);
end
imageHT = HalftoningED(imageIn);

type = 1;
[shares, stacked] = RandomGridKafri(imageHT, type);
figure; imshow(shares(:,:,1),[]);
figure; imshow(shares(:,:,2),[]);
figure; imshow(stacked,[]);




