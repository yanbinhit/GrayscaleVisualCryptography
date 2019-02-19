function Naor22Test()
% Use binary secret iamge
imageIn= imread('../Images/geometry.bmp','bmp'); % Use binary image as secret

[shares, stacked] = Naor22(imageIn);

figure;
imshow(imageIn,[]);
figure;
imshow(shares(:,:,1),[]);
figure;
imshow(shares(:,:,1),[]);
figure;
imshow(stacked,[]);
