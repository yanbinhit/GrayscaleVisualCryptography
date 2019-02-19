function HalftoneVcViaErrorDiff2_2outof2()
% Implementation of a 2-out-of-2 VC scheme, with meanful shaddows
% The VC schemem is based on construction 2 in "Halftone Visual Cryptography Via Error
% Diffusion" 2009.
close all; clear all; 
Q = 4; % each block is Q-by-Q
gamma = 4; % number of SIPs
nAbps = 6; % number of ABPs in each share
nShares = 2;
nAbpsPerShare = 6;
nRowsCover = 512;
nColsCover = 512;


%diffusionMask = [0 0 0 8 4; 2 4 8 4 2; 1 2 4 2 1]/42;  % Stucki
diffusionMask = [0 0 0 7 0; 0 3 5 1 0; 0 0 0 0 0]/16;  % Floyd-Steinberg
%diffusionMask = [0 0 0 7 5; 3 5 7 5 3; 1 3 5 3 1]/48;  % Jarvis

% Generate locations for SIPs, and ABPs for the 2 shares, using random method.
% (not optimal)
[sips, abps] = AbpSipLocatoins2outof2Random(nRowsCover, nColsCover);
save 'temp2outof2.mat' sips abps;
load 'temp2outof2.mat';

shares = zeros(nRowsCover, nColsCover, 2); % each slice is a share

% Cover images
cover = zeros(nRowsCover, nColsCover, 2); % 3 cover images
cover(:,:,1) = double(rgb2gray(imread('../Images/lena.tiff','tiff')));  % cover image 1
cover(:,:,2) = double(rgb2gray(imread('../Images/peppers.tiff','tiff'))); % cover image 2
coverOriginal = cover;

% secretImgIn = rgb2gray(imread('abc.bmp','bmp'));
% secretImg = secretImgIn > 0; % The secret image to be shared
% figure; imshow(secretImg,[]);
secretImgIn = imread('../Images/geometry.bmp','bmp');
secretImgIn = imresize(secretImgIn,[nRowsCover/4 nColsCover/4],'nearest');
secretImg = secretImgIn > 0;

% set of matices in visual cryptography
M0 = zeros(2,2);          % set of permutated elementary matrix for WHITE pixel
M1 = zeros(2,2);          % set of permutated elementary matrix for BLACK pixel
M0 = [1 1 0 0; 1 1 0 0];
M1 = [1 1 0 0; 0 0 1 1];

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
        
        for k = 1:2  % fill in sips for each share
            share((i-1)*Q+iIdxSip(1), (j-1)*Q+jIdxSip(1), k) = M(k,1)*255;
            share((i-1)*Q+iIdxSip(2), (j-1)*Q+jIdxSip(2),k ) = M(k,2)*255;    
            share((i-1)*Q+iIdxSip(3), (j-1)*Q+jIdxSip(3),k ) = M(k,3)*255;    
            share((i-1)*Q+iIdxSip(4), (j-1)*Q+jIdxSip(4),k ) = M(k,4)*255;    
        end
        
    end
end

% Fill in the ABPs for each share
for k = 1:2
    for i = 1:nRowsCover
        for j = 1:nColsCover
            if abps(i,j,k) == 1
                share(i,j,k) = 0;  % 填入黑色
            end
        end
    end
end


% Error Diffusion for the information carrying pixels in each share
for k = 1:nShares
    for i = 3:nRowsCover-2
        for j = 3:nColsCover-2
            
            blkIdxI = i:i+2;
            blkIdxJ = j-2:j+2;
            
            diffMaskWeighted = diffusionMask .* ...
                (1 - sips(blkIdxI,blkIdxJ)) .* (1- abps(blkIdxI, blkIdxJ, k));
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
for k = 1:nShares
 subplot(1,2,k)
 imshow(share(:,:,k),[]); title(['share' num2str(k)]);
end
figure; imshow(share(:,:,1),[]); 
figure; imshow(share(:,:,2),[]); 

% print and stacking
% 模拟打印和叠加前，需要将黑色用1表示
share = ~share;           % for printing&stacking, interpret '1' as black
stacked = (share(:,:,1) | share(:,:,2));  % simulate stacking ( OR operation on printed dots)
figure; imshow(stacked,[]); title('printed');

% 在电脑屏幕上显示解码结果，需要将黑色用0表示
stacked = ~stacked;             % to display on screen, change '1' to '0' for black
figure; imshow(stacked,[]); 
%title('displayed on screen');

