function RandomGridChenTsao2009nnTest()
close all; clear;
n=3;
type = 1; % type 1 algorithm from Kafri

prevDir = pwd;
[dir, dummy, dummy2] = fileparts(mfilename('fullpath'));
cd(dir);
imageIn= imread('../Images/geometry.bmp','bmp');
[shares, stacked]= RandomGridChenTsao2009nn(imageIn, n, type);
for i= 1:n
    figure;
    imshow(shares(:,:,i));
    title(['share ' num2str(i)]);
end

% stacking of share 1 and 2
target12 = ~(~shares(:,:,1) | ~shares(:,:,2));
figure; imshow(target12,[]);
title('target 12');
% stacking of share 1 and 3
target13 = ~(~shares(:,:,1) | ~shares(:,:,3));
figure; imshow(target13,[]);
title('target 13');
% target
figure; 
imshow(stacked,[]); title('stacked');

