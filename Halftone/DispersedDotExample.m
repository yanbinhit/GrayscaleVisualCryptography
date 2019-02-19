function DispersedDotExample()
close all; clear all; 

I2 = [1,2;3,0];
I4 = [4*I2+1, 4*I2+2; 4*I2+3, 4*I2];
I8 = [4*I4+1, 4*I4+2; 4*I4+3, 4*I4];
I = I8;
% Display dithering matrix
figure;
imagesc(I);
colormap gray;
axis off;
axis square;

% Halftoing reuslt
L= 8;
T = round(255* (I+0.5)/(L^2))

img = zeros(128, 256);
for i = 1:256
    img(:,i) = 256-i;
end

figure; 
imshow(img,[]);

[M, N] = size(img);
imgHt = zeros(M,N);

for i=0:L:M-L
    for j = 0:L:N-L
        imgHt(i+1:i+8,j+1:j+8) = img(i+1:i+8,j+1:j+8) > T;
    end
end
figure;
imshow(imgHt,[]);