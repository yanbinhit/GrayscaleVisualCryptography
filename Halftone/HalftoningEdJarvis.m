function imageHT = HalftoningEdJarvis(imageIn)
%=======================================================
% Error diffusion based halftoning for color images using Jarvis kernel
%=======================================================
[M, N, S] = size(imageIn); 
imageIn = double(imageIn);
imageHT = zeros(M+4, N+4, S);
for k = 1:S
    imageColor(:,:,k) = padarray(imageIn(:,:,k), [2 2], 'replicate', 'both');
end

for k = 1: S
    for i = 3:M+2
        for j = 3:N+2
            if imageColor(i,j,k) > 127
                imageHT(i,j,k) = 255;
            else
                imageHT(i,j,k) = 0;
            end
            error = imageColor(i,j,k) - imageHT(i,j,k);
            
            % error diffusion
            imageColor(i,j+1,k) = imageColor(i,j+1,k) + error * (7/48);
            imageColor(i,j+2,k) = imageColor(i,j+2,k) + error * (5/48);
            imageColor(i+1,j-2, k) = imageColor(i+1, j-2,k) + error * (3/48);
            imageColor(i+1,j-1, k) = imageColor(i+1, j-1, k) + error * (5/48);
            imageColor(i+1,j, k) = imageColor(i+1, j, k)+ error * (7/48);
            imageColor(i+1, j+1, k) = imageColor(i+1, j+1, k) + error*(5/48);
            imageColor(i+1, j+2, k) = imageColor(i+1, j+2, k) + error*(3/48);
            imageColor(i+2,j-2, k) = imageColor(i+2, j-2,k) + error * (1/48);
            imageColor(i+2,j-1, k) = imageColor(i+2, j-1, k) + error * (3/48);
            imageColor(i+2,j, k) = imageColor(i+2, j, k)+ error * (5/48);
            imageColor(i+2, j+1, k) = imageColor(i+2, j+1, k) + error*(3/48);
            imageColor(i+2, j+2, k) = imageColor(i+2, j+2, k) + error*(1/48);
          end
    end
end
imageHT = uint8(imageHT(3:M+2, 3:N+2, :));
