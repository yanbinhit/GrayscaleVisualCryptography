function HalftoneVcViaErrorDiff2Auth2sipsBatchTest()
% 对本文算法做批量测试
close all; clear all;
DEBUG = 0;
nRowsCover = 768;
nColsCover = 768;
Q = 4;

nTests = 100; % total number of tests
contrastSecret = zeros(nTests,1);   % record the batch test result
mssimSecret = zeros(nTests,1);
mssimAvgShare = zeros(nTests,1);
psnrAvgShare = zeros(nTests,1);

prevDir = pwd;
[dir, dummy, dummy2] = fileparts(mfilename('fullpath'));
cd(dir);

fh = waitbar(0,'Processing ...');
for ii = 1:nTests
    waitbar((ii-1)/nTests,fh);
    %==================================================================
    % Randomly select 2 cover images from the set of cover images
    %==================================================================
    nCoverImgs = 10; % total number of cover images used in the test
    cover = zeros(nRowsCover, nColsCover, 2); % 2 cover images
    img1Name = [dir '\coverImages\' num2str(randi(nCoverImgs)) '.bmp'];
    img1 = (imread(img1Name,'bmp'));  % cover image 1
    if size(img1,3)>1
        img1 = rgb2gray(img1);
    end
    cover(:,:,1) = imresize(img1,[nRowsCover nColsCover]);
    img2Name = [dir '\coverImages\' num2str(randi(nCoverImgs)) '.bmp'];
    img2 = (imread(img2Name,'bmp')); % cover image 2
    if size(img2,3)>1
        img2 = rgb2gray(img2);
    end
    cover(:,:,2) = imresize(img2,[nRowsCover nColsCover]);
    
%     figure; imshow(cover(:,:,1),[]);
%     figure; imshow(cover(:,:,2),[]);
    
    %====================================================================
    % Randomly choose a secret image and a watermark
    %====================================================================
    nSecretImgs = 6;
    nWatermarkImgs = 4;
    secretImgName = [dir '\secretImagesLarge\' num2str(randi(nSecretImgs)) '.bmp'];
    secretImgIn = (imread(secretImgName,'bmp'));  % cover image 1
    if size(secretImgIn,3)>1
        secretImgIn = rgb2gray(secretImgIn);
    end
    secretImgIn = imresize(secretImgIn, [nRowsCover nColsCover]);
    secretImgIn = (secretImgIn>0);
    
    watermarkImgName = [dir '\watermarkImages\' num2str(randi(nWatermarkImgs)) '.bmp'];
    watermarkImgIn = (imread(watermarkImgName,'bmp')); % cover image 2
    if size(watermarkImgIn,3)>1
        watermarkImgIn = rgb2gray(watermarkImgIn);
    end
    watermarkImgIn = imresize(watermarkImgIn, [nRowsCover nColsCover/2],'nearest');
    watermarkImgIn = (watermarkImgIn>0);
    
%     figure; imshow(secretImgIn,[]);
%     figure; imshow(watermarkImgIn,[]);
    
    %======================================================================
    % generate shares and simulate stacking
    %======================================================================
    % generate shares and simulate the stacking operations
    tic
    [shares, stacked, wmRecovered] = ...
       HalftoneVcViaErrorDiff2Auth2sipsFun(cover, secretImgIn,...
                                            watermarkImgIn ,...
                                            nRowsCover,...
                                            nColsCover, DEBUG);
    toc
    %secretImgInUp = imresize(secretImgIn, [nRowsCover nColsCover]);
    %secretImgInUp = (secretImgInUp > 0);
    contrastSecret(ii) = SecretContrast(secretImgIn, stacked);
    [mssim, ssim_map] = ssim(secretImgIn*255, stacked*255); % input dynamic range 0-255
    mssimSecret(ii) = mssim;
    
    if DEBUG
        mssim
        figure;
        subplot(1,2,1); imshow(secretImgIn,[]);
        subplot(1,2,2); imshow(stacked,[]);
        figure; imshow(secretImgIn,[]);
        figure; imshow(stacked,[]);
    end
    %----------------------------------------
    % Compare quality
    %----------------------------------------
    % Compare the PSNR of two filtered images
    sigma = 2;
    hSize = 11;
    h = fspecial('gaussian', hSize, sigma);
    for k = 1:2
        sharesFilt(:,:,k) = imfilter(shares(:,:,k), h, 'replicate');
    end
    
    psnr = zeros(1,2);
    for k = 1:2
        psnr(k) = calPSNR(cover(:,:,k), sharesFilt(:,:,k));
    end
    psnrAvg = sum(psnr)/numel(psnr);
    psnrAvgShare(ii) = psnrAvg;
    
%     figure;subplot(1,2,1); imshow(sharesFilt(:,:,1),[]);
%     subplot(1,2,2); imshow(sharesFilt(:,:,2),[]);
    
    % calculate SSIM
    for k = 1:2
        [mssim(k), dummy] = ssim(cover(:,:,k), shares(:,:,k));
    end
    mssimAvgShare(ii) = (mssim(1) + mssim(2))/2;
    
end % end of ii

close(fh);

save('HalftoneVcViaErrorDiff2Auth2sipsBatchTest.mat', 'contrastSecret',...
    'mssimSecret', 'mssimAvgShare', 'psnrAvgShare');
disp('Done!');
