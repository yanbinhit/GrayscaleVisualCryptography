function  ErrorDiffErrorImage(imageIn)
%=======================================================
% Display the error image, scatter plot between error image and modified
% input image.
%=======================================================
close all; clear all;
imageIn = imread('../Images/lena.bmp','bmp');
[M, N, S] = size(imageIn); 
imageIn = double(imageIn);
imageHT = zeros(M+2, N+2, S);
imageError = zeros(M+2, N+2, S);
for k = 1:S
    imageColor(:,:,k) = padarray(imageIn(:,:,k), [1 1], 'replicate', 'both');
end

for k = 1: S
    for i = 2:M+1
        for j = 2:N+1
            if imageColor(i,j,k) > 127
                imageHT(i,j,k) = 255;
            else
                imageHT(i,j,k) = 0;
            end
            error = imageColor(i,j,k) - imageHT(i,j,k);
            imageError(i,j,k) = -error;
                        
            % error diffusion
            imageColor(i,j+1,k) = imageColor(i,j+1,k) + error * (7/16);
            imageColor(i+1,j,k) = imageColor(i+1,j,k) + error * (5/16);
            imageColor(i+1,j-1,k) = imageColor(i+1,j-1,k) + error * (3/16);
            imageColor(i+1,j+1,k) = imageColor(i+1,j+1,k) + error * (1/16);
        end
    end
end
imageHT = uint8(imageHT(2:M+1, 2:N+1, :));
imageError = imageError(2:M+1, 2:N+1, :);
imageColor = imageColor(2:M+1, 2:N+1, :);

figure;
imshow(uint8(imageError+127));
figure;
plot(imageError(:));
figure;
hist(imageError(:)+127,40);
