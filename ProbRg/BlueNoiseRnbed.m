function imageHT = BlueNoiseRnbed(M, N, nR)
% generate a blude noise pattern using Noise Balanced Error Diffusion
% (NBED)
img = 1/2 * ones(M,N);
imageHT = zeros(M+2, N+2);
imagePad(:,:) = padarray(img(:,:), [1 1], 'replicate', 'both');
for i = 2:M+1
    for j = 2:N+1
        sai = (2*(randi(2)-1)-1);
        imagePad(i,j) = imagePad(i,j) + nR * sai;
        if imagePad(i,j) > 1/2
            imageHT(i,j) = 1;
        else
            imageHT(i,j) = 0;
        end
        error = imagePad(i,j) - imageHT(i,j) - nR*sai;
        
        % error diffusion
        imagePad(i,j+1) = imagePad(i,j+1) + error * (7/16);
        imagePad(i+1,j) = imagePad(i+1,j) + error * (5/16);
        imagePad(i+1,j-1) = imagePad(i+1,j-1) + error * (3/16);
        imagePad(i+1,j+1) = imagePad(i+1,j+1) + error * (1/16);
    end
end
imageHT = imageHT(2:M+1, 2:N+1);
