function HalftoneVcViaErrorDiff2Auth2sips()
% Implementation of a 2-out-of-2 VC scheme, with meanful shaddows, with
% watermark authentication
% The VC schemem is based on construction 2 in "Halftone Visual Cryptography Via Error
% Diffusion" 2009.
% consider 2sips
close all; clear all; 
Q = 4; % each block is Q-by-Q
gamma = 2; % number of SIPs
nShares = 2;
nAbpsPerShare = 7; % number of ABPs in each share
nRowsCover = 768;
nColsCover = 768;


%diffusionMask = [0 0 0 8 4; 2 4 8 4 2; 1 2 4 2 1]/42;  % Stucki
diffusionMask = [0 0 0 7 0; 0 3 5 1 0; 0 0 0 0 0]/16;  % Floyd-Steinberg
%diffusionMask = [0 0 0 7 5; 3 5 7 5 3; 1 3 5 3 1]/48;  % Jarvis

% Generate locations for SIPs, and ABPs for the 2 shares, using random method.
% (not optimal)
[sips, abps] = AbpSipLocatoins2outof2Random2sips(nRowsCover, nColsCover);
sips(:,(nColsCover/2+1) : nColsCover) = sips(:,1:nColsCover/2);      % copy right half to left half
abps(:,(nColsCover/2+1) : nColsCover,1) = abps(:, 1:(nColsCover/2), 1) ;
abps(:,(nColsCover/2+1) : nColsCover,2) = abps(:, 1:(nColsCover/2), 2) ;

save 'temp2outof2.mat' sips abps;
load 'temp2outof2.mat';

share = zeros(nRowsCover, nColsCover, 2); % each slice is a share

% Cover images
cover = zeros(nRowsCover, nColsCover, 2); % 2 cover images
img1 = double(rgb2gray(imread('../Images/lena.tiff','tiff')));  % cover image 1
cover(:,:,1) = imresize(img1,[nRowsCover nColsCover]);
img2 = double(rgb2gray(imread('../Images/peppers.tiff','tiff'))); % cover image 2
cover(:,:,2) = imresize(img2,[nRowsCover nColsCover]);
coverOriginal = cover;
figure; imshow(cover(:,:,1),[]); title('Cover 1');
figure; imshow(cover(:,:,2),[]); title('Cover 2');
imwrite(uint8(cover(:,:,1)),'cover1.bmp','bmp');
imwrite(uint8(cover(:,:,2)),'cover2.bmp','bmp');



secretImgIn = imread('../Images/geometry.bmp','bmp');
secretImgIn = imresize(secretImgIn, [nRowsCover/4 nColsCover/4], 'nearest');
secretImg = secretImgIn > 0; % The secret image to be shared
watermarkImgIn = (imread('../Images/wen.bmp','bmp'));
watermarkImgIn = imresize(watermarkImgIn, [[nRowsCover/4 nColsCover/8]],'nearest');
watermarkImg = watermarkImgIn > 0; 
figure; imshow(secretImg,[]);
figure; imshow(watermarkImg,[]);
imwrite(secretImg, 'secretImg.bmp','bmp');
imwrite(watermarkImg, 'watermarkImg.bmp','bmp');


% set of matices in visual cryptography
M0 = zeros(2,2);          % set of permutated elementary matrix for WHITE pixel
M1 = zeros(2,2);          % set of permutated elementary matrix for BLACK pixel
M0 = [1 0; 1 0];
M1 = [1 0; 0 1];

% Fill in the SIPs in the halftone share images
permMt = zeros(nRowsCover/Q, nColsCover/Q/2, 2); % each position is a permutation vector
% For the left half of both shares
for i = 1:ceil(nRowsCover/Q)
    for j = 1:ceil(nColsCover/Q/2)
        cellIdxI = ((i-1)*Q+1):(i*Q);
        cellIdxJ = ((j-1)*Q+1):(j*Q);
        [iIdxSip,jIdxSip] = find(sips(cellIdxI, cellIdxJ)==1);
        
        permMt(i,j,:) = randperm(2); % a random permuation of columns
        if secretImg(i,j) == 1   % the secret pixel is a white pixel
            M = M0(:,permMt(i,j,:));
        else
            M = M1(:,permMt(i,j,:));
        end
        
        for k = 1:2  % fill in sips for each share
            share((i-1)*Q+iIdxSip(1), (j-1)*Q+jIdxSip(1),k ) = M(k,1)*255;
            share((i-1)*Q+iIdxSip(2), (j-1)*Q+jIdxSip(2),k ) = M(k,2)*255;    
        end
        
    end
end

% For the right half of share 2
for i = 1:ceil(nRowsCover/Q)
    for j = 1:ceil(nColsCover/Q/2)
        cellIdxI = ((i-1)*Q+1):(i*Q);
        cellIdxJ = ((j-1)*Q+1):(j*Q);
        [iIdxSip,jIdxSip] = find(sips(cellIdxI, cellIdxJ)==1);
        
        shift = nColsCover/2;
        if watermarkImg(i,j) == 1 % the secret pixel is a white pixel, 填相同
            share((i-1)*Q+iIdxSip(1), (j-1)*Q+jIdxSip(1)+shift, 2) = ...
                share((i-1)*Q+iIdxSip(1), (j-1)*Q+jIdxSip(1), 1);
            share((i-1)*Q+iIdxSip(2), (j-1)*Q+jIdxSip(2)+shift, 2 ) = ...
                share((i-1)*Q+iIdxSip(2), (j-1)*Q+jIdxSip(2), 1 );
        else                                                           % 填互补
            share((i-1)*Q+iIdxSip(1), (j-1)*Q+jIdxSip(1)+shift, 2) = ...
                255 - share((i-1)*Q+iIdxSip(1), (j-1)*Q+jIdxSip(1), 1);
            share((i-1)*Q+iIdxSip(2), (j-1)*Q+jIdxSip(2)+shift, 2 ) = ...
                255 - share((i-1)*Q+iIdxSip(2), (j-1)*Q+jIdxSip(2), 1 );
        end
    end
end

% for the right half of share 1
for i = 1:ceil(nRowsCover/Q)
    for j = 1:ceil(nColsCover/Q/2)
        cellIdxI = ((i-1)*Q+1):(i*Q);
        cellIdxJ = ((j-1)*Q+1):(j*Q);
        [iIdxSip,jIdxSip] = find(sips(cellIdxI, cellIdxJ)==1);     
        shift = nColsCover/2;
        if secretImg(i,j+(ceil(nColsCover/Q/2))) == 1   % the secret pixel is a white pixel
            share((i-1)*Q+iIdxSip(1), (j-1)*Q+jIdxSip(1)+shift, 1) = ...
                share((i-1)*Q+iIdxSip(1), (j-1)*Q+jIdxSip(1)+shift, 2);
            share((i-1)*Q+iIdxSip(2), (j-1)*Q+jIdxSip(2)+shift, 1 ) = ...
                share((i-1)*Q+iIdxSip(2), (j-1)*Q+jIdxSip(2)+shift, 2 );
        else
            share((i-1)*Q+iIdxSip(1), (j-1)*Q+jIdxSip(1)+shift, 1) = ...
                255 - share((i-1)*Q+iIdxSip(1), (j-1)*Q+jIdxSip(1)+shift, 2);
            share((i-1)*Q+iIdxSip(2), (j-1)*Q+jIdxSip(2)+shift, 1 ) = ...
                255 - share((i-1)*Q+iIdxSip(2), (j-1)*Q+jIdxSip(2)+shift, 2 );
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
    threshold = (max(max(cover(:,:,k))) - min(min(cover(:,:,k))))/2
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
            %threshold = 127 - *coverOriginal(i,j,k); % input threshold modulation
            %threshold = 127; % no threshold modulation          
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
figure; imshow(share(:,:,1),[]); title('share 1');
figure; imshow(share(:,:,2),[]); title('share 2');

%--------------------------------------------
% attack share 1
%--------------------------------------------
% % attack #1: known sips positions
% share1Att = share(:,:,1);
% share1Att(sips==1) = 0;
% figure; imshow(share1Att,[]); title('attacked share 1');
% share(:,:,1) = share1Att;

% attack #2: unknown sips positions, randomly choose \gamma pixels and set them 0
rndPositions = zeros(nRowsCover, nColsCover);
rndPositions = rand([nRowsCover nColsCover])<(gamma/Q^2);
share1Att = share(:,:,1);
share1Att(rndPositions==1) = 0;
figure; imshow(share1Att,[]); title('attacked share 1');
share(:,:,1) = share1Att;

% 保存文件画图用
imwrite(share(:,:,1),'share1.bmp','bmp');
imwrite(share(:,:,2),'share2.bmp','bmp');


%----------------------------------------------------
% 对掩盖图像做普通半色调化，用于对比视觉质量
%----------------------------------------------------
cover = coverOriginal;
htImg = zeros(nRowsCover, nColsCover, 2);
for k = 1:2
    for i = 3:nRowsCover-2
        for j = 3: nColsCover-2
            blkIdxI = i:i+2;
            blkIdxJ = j-2:j+2;
            htImg(i,j,k) = 255 * (cover(i,j,k)>threshold);
            err1 = htImg(i,j,k) - cover(i,j,k);
            cover(blkIdxI, blkIdxJ, k) = cover(blkIdxI, blkIdxJ, k) ...
                - err1 .* diffusionMask;
        end
    end
end
save 'htImgs2sips.mat' htImg share;

figure;imshow(htImg(:,:,1),[]); title('Halftone image');
figure; imshow(htImg(:,:,2),[]); title('Halftone image');
%--------------------------------------------------
% print and stacking
% 模拟打印和叠加前，需要将黑色用1表示
%--------------------------------------------------
share = ~share;           % for printing&stacking, interpret '1' as black
stacked = (share(:,:,1) | share(:,:,2));  % simulate stacking ( OR operation on printed dots)
figure; imshow(stacked,[]); title('printed');
shiftStack = zeros(nRowsCover, nColsCover+nColsCover/2);
shiftStack = ([share(:,:,2) zeros(nRowsCover, nColsCover/2)]) | ([zeros(nRowsCover,nColsCover/2) share(:,:,1)]);
figure; imshow(shiftStack,[]); title('shifted and stacked');


% 在电脑屏幕上显示解码结果，需要将黑色用0表示
stacked = ~stacked;             % to display on screen, change '1' to '0' for black
figure; imshow(stacked,[]); title('displayed on screen');
shiftStack = ~shiftStack;
figure; imshow(shiftStack,[]); title('shifted and stacked');

imwrite(stacked,'stackedImg.bmp','bmp');
imwrite(shiftStack,'shiftStack.bmp','bmp');

% 保存文件以便于制作动画
save 'shares.mat' share;



