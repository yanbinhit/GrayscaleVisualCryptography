function [shares, stacked]= RandomGridChenTsao2009nn(imageIn, n, type)
% Implement the Chen and Tsao 2009 paper on (n,n) scheme for random grid. 
% 2009 Visual secret sharing by random grids revisited
[nR, nC] = size(imageIn);
shares = zeros(nR, nC, n);
img = HalftoningED(imageIn);
for i=1:n-1
    [shares2, stacked2] = RandomGridKafri(img, type);
    shares(:,:,i) = shares2(:,:,1);
    shares(:,:,i+1) = shares2(:,:,2);
    img = shares(:,:,i+1);
end

stacked = ones(nR, nC);
for i = 1:n
    stacked = ~(~shares(:,:,i) | ~stacked);
end

