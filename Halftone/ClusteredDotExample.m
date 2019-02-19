function ClusteredDotExample()
close all; clear all;

% Index matrix for clustered dot dithering
I = [62  57  48  36  37  49  58  63; ...
   56  47  35  21  22  38  50  59;  ...
   46  34  20  10  11  23  39  51;  ...
   33  19  9  3  0  4  12  24;  ...
   32  18  8  2  1  5  13  25;  ...
   45  31  17  7  6  14  26  40;  ...
   45  44  30  16  15  27  41  52;  ...
   61  54  43  29  28  42  53  60]
L = 8;

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

