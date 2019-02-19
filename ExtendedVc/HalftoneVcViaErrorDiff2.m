function HalftoneVcViaErrorDiff2()
% Implement the scheme 2 in the paper "Halftone Visual Cryptography via Error
% Diffusioin", a 3-out-of-3 scheme.
% 注意：不要对误差扩散滤波器归一化。没有阈值调制．
close all; clear all; 

Q = 3;  % patch size is 3-by-3
gamma = 4;
noSipsPerShare = 2;
noShares = 3;

%diffusionMask = [0 0 0 8 4; 2 4 8 4 2; 1 2 4 2 1]/42;  % Stucki
diffusionMask = [0 0 0 7 0; 0 3 5 1 0; 0 0 0 0 0]/16;  % Floyd-Steinberg
%diffusionMask = [0 0 0 7 5; 3 5 7 5 3; 1 3 5 3 1]/48;  % Jarvis

% Generate locations for SIPs, and ABPs for the 3 shares, using random method.
% (not optimal)
[sips, abps] = AbpSipLocatoins3outof3Random;
save 'temp.mat' sips abps;
load 'temp.mat';


nRowsCover = 513;
nColsCover = 513;
shares = zeros(nRowsCover, nColsCover, 3); % each slice is a share

% Cover images
cover = zeros(nRowsCover, nColsCover, 3); % 3 cover images
cover1 = double(rgb2gray(imread('lena.tiff','tiff')));  % cover image 1
cover(:,:,1) = padarray(cover1, [1 1], 'replicate', 'post');
cover2 = double(rgb2gray(imread('peppers.tiff','tiff'))); % cover image 2
cover(:,:,2) = padarray(cover2, [1 1], 'replicate', 'post');
cover3 = double(rgb2gray(imread('baboon.tiff','tiff'))); % cover image 3
cover(:,:,3) = padarray(cover3, [1 1], 'replicate', 'post');
coverOriginal = cover;
clear cover1; 
clear cover2; 
clear cover3;

secretImgIn = rgb2gray(imread('abc.bmp','bmp'));
secretImgIn = imresize(secretImgIn, [171 171] );
secretImg = secretImgIn > 0; % The secret image to be shared
figure; imshow(secretImg,[]);
% Fill in the secret pixels in the secret images


% set of matices in visual cryptography
M0 = zeros(2,2);          % set of permutated elementary matrix for WHITE pixel
M1 = zeros(2,2);          % set of permutated elementary matrix for BLACK pixel
M0 = [1 1 0 0; 1 1 0 0; 1 1 0 0];
M1 = [1 1 0 0; 0 1 1 0; 0 0 1 1];

% Fill in the SIPs in halftone share images
for i = 1:ceil(nRowsCover/Q)
    for j = 1:ceil(nColsCover/Q)
        cellIdxI = ((i-1)*Q+1):(i*Q);
        cellIdxJ = ((j-1)*Q+1):(j*Q);
        [iIdxSip,jIdxSip] = find(sips(cellIdxI, cellIdxJ)==1);
        
        p = randperm(4); % a random permuation of columns
        if secretImg(i,j) == 1   % the secret pixel is a white pixel
            M = M0(:,p);
        else
            M = M1(:,p);
        end
        
        for k = 1:3  % fill in sips for each share
            share((i-1)*Q+iIdxSip(1), (j-1)*Q+jIdxSip(1), k) = M(k,1)*255;
            share((i-1)*Q+iIdxSip(2), (j-1)*Q+jIdxSip(2),k ) = M(k,2)*255;    
            share((i-1)*Q+iIdxSip(3), (j-1)*Q+jIdxSip(3),k ) = M(k,3)*255;    
            share((i-1)*Q+iIdxSip(4), (j-1)*Q+jIdxSip(4),k ) = M(k,4)*255;    
        end
        
    end
end

% Fill in the ABPs for each share
for k = 1:3
    for i = 1:nRowsCover
        for j = 1:nColsCover
            if abps(i,j,k) == 1
                share(i,j,k) = 0;
            end
        end
    end
end


% Error Diffusion for the information carrying pixels in each share
for k = 1:noShares
    for i = 4:nRowsCover-2
        for j = 4:nColsCover-2
            
            blkIdxI = i:i+2;
            blkIdxJ = j-2:j+2;
            
            diffMaskWeighted = diffusionMask .* ...
                (1 - sips(blkIdxI,blkIdxJ)) .* (1- double(abps(blkIdxI, blkIdxJ, k)));
            %diffMaskWeighted = diffMaskWeighted ./ (sum(diffMaskWeighted(:))); 不能归一化
            
            if (sips(i,j) == 1) | (abps(i,j,k) == 1)   % SIPs 和 ABPs处保持填过的值不变，
                % 将误差扩散到周围
                
                err1 = share(i,j,k) - cover(i,j,k);
                cover(blkIdxI, blkIdxJ, k) = cover(blkIdxI, blkIdxJ, k) ...
                    - err1 .* diffMaskWeighted;
                
                continue;
            end
            
            % 对于其他Pixel正常误差扩散
            % % using threshold modulation
            %threshold = 255*(0.25+0.33*0.25*(share(i,j-1,k)+share(i,j-2,k)+share(i,j-3,k)));
            threshold = 127; % no threshold modulation
            share(i,j,k) = 255 * (cover(i,j,k)>threshold);
            err1 = share(i,j,k) - cover(i,j,k);
            cover(blkIdxI, blkIdxJ, k) = cover(blkIdxI, blkIdxJ, k) ...
                - err1 .* diffMaskWeighted;
            
        end
    end
end


% show shares
figure;
for k = 1:noShares
 subplot(2,2,k)
 imshow(share(:,:,k),[]); title(['share' num2str(k)]);
end
figure; imshow(share(:,:,1),[]); 
figure; imshow(share(:,:,2),[]); 
figure; imshow(share(:,:,3),[]); 
% print and stacking
% 模拟打印和叠加前，需要将黑色用1表示
share = ~share;           % for printing&stacking, interpret '1' as black
stacked = (share(:,:,1) | share(:,:,2)) | share(:,:,3);  % simulate stacking ( OR operation on printed dots)
figure; imshow(stacked,[]); title('printed');

% 在电脑屏幕上显示解码结果，需要将黑色用0表示
stacked = ~stacked;             % to display on screen, change '1' to '0' for black
figure; imshow(stacked,[]); title('displayed on screen');
