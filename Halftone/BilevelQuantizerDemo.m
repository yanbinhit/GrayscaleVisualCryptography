function BilevelQuantizer()
close all; clear all;
imageIn = imread('../Images/lena.bmp','bmp');

T = 127;
imageOut = 255*(imageIn>127);
figure;
imshow(imageIn,[]);
figure;
imshow(imageOut,[]);