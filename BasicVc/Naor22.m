function [shares, stacked] = Naor22(imageIn)
% Implement a (2,2) Naor scheme with size expansion of 4 (so that the
% aspect ratio of the image doesn't change. 

M0 = [1 1 0 0; 1 1 0 0];
M1 = [1 1 0 0; 0 0 1 1];

[nR, nC] = size(imageIn);
imageIn = imageIn>0;
shares = zeros(2*nR, 2*nC, 2);
stacked = zeros(2*nR, 2*nC);

for i = 1:nR
    for j = 1:nC
        p = randperm(4); % a random permuation of columns
        if imageIn(i,j) == 1   % the secret pixel is a white pixel
            M = M0(:,p);
        else
            M = M1(:,p);
        end
        ii = ((i-1)*2+1): (2*i);
        jj = ((j-1)*2+1): (2*j);
        shares(ii,jj,1) = reshape(M(1,:),2,2)*255;
        shares(ii,jj,2) = reshape(M(2,:),2,2)*255;
    end
end

shares2 = shares>0;
stacked = (~(~shares2(:,:,1) | ~shares2(:,:,2)))*255;