function [shares, stacked] = ...
    FrgvssChen2011(cover, secretImgIn, nRowsCover, nColsCover)
% Implementation of friendly random grid (RG) based visual secret sharing (VSS):
% FRGVSS, proposed by T. H. Chen(2011) in "User-Friendly Random-Grid-Based
% Visual Secret Sharing".
% 
% Assumption of the inputs: cover image is gray scale, secret image is binary
%<Inputs>
% cover: cover image cube, each slice is a cover image
% secretImgIn: input secret image 
% nRowsCover: number of rows
% nColsCover: number of columns 
% 
%<outpus>
% shares: share image cube, each slice is a share
% stacked: stacked shares (recovered secret image)
% 
% Bin Yan, yanbinhit@hotmail.com 
% 2015.5.11
%


% Important parameters
nShares = 2;
%gamma = 1/3; % value taken by Chen's 2011 paper
gamma = 2/16; % fraction of sips over all pixels, for comparison with our algorithm
%----------------------------
% Pre-process cover image
%----------------------------
% Cover images
coverHt = zeros(nRowsCover, nColsCover, 2); % halftone images
coverOriginal = cover;
% figure; imshow(cover(:,:,1),[]); title('Cover 1');
% figure; imshow(cover(:,:,2),[]); title('Cover 2');
% halftone the cover
for i = 1:nShares
    coverHt(:,:,i) = HalftoningED(cover(:,:,i));
end
% figure; imshow(coverHt(:,:,1),[]); title('Cover 1: Halftone');
% figure; imshow(coverHt(:,:,2),[]); title('Cover 2: Halftone');

%-----------------------------
% Pre-process the secret image
%-----------------------------
% Make the size of the secret image the same as that of the cover images
secretImg = secretImgIn > 0; % The secret image to be shared
%figure; imshow(secretImg,[]); 

%---------------------------------------------------
% seclect positions to embed secret or cover image
%----------------------------------------------------
% 产生从1到Mc-by-Nc的置乱，然后取前beta作为secret pixel的位置，剩下的一半作为
% cover1的位置，剩下的另外一半用于承载cover2.最后将线性位置转化为二维下标，
% 并分别用三个指示矩阵指示这些位置
nPixels = nRowsCover * nColsCover;
sips = zeros(nRowsCover, nColsCover); % indicator matrix for sips
cips1 = zeros(nRowsCover, nColsCover);% indicator matrix for cover information pixels 1
cips2 = zeros(nRowsCover, nColsCover); % indicator matrix for cover information pixels 2
rndVec = randperm( nPixels );
[idxI, idxJ] = ind2sub([nRowsCover nColsCover], rndVec);
pim = zeros(nRowsCover, nColsCover); % pixel indicator matrix, 0: sip, 1:cip1, 2:cip2
                                     % PIM is another way of indicating sips and cips
for k = 1: floor(nPixels*gamma)
    sips(idxI(k),idxJ(k)) = 1;
    pim(idxI(k),idxJ(k)) = 0;
end
for k = (floor(nPixels*gamma)+1):(floor(nPixels*(gamma+(1-gamma)/2)))
    cips1(idxI(k), idxJ(k)) = 1;
    pim(idxI(k),idxJ(k)) = 1;
end
for k = (floor(nPixels*(gamma+(1-gamma)/2))+1) : nPixels
   cips2(idxI(k), idxJ(k)) = 1; 
   pim(idxI(k),idxJ(k)) = 2;
end
%figure; imshow(pim,[]);

%=================================
% Encode to shares
%=================================
shares = zeros(nRowsCover, nColsCover, nShares);
for i = 1:nRowsCover
    for j = 1:nColsCover
        if pim(i,j) == 0 % sips, use random grid VSS
            shares(i,j,1) = normrnd(0,1,1,1)>0;
            if secretImg(i,j) == 1 % white 
                shares(i,j,2) = shares(i,j,1);
            else
                shares(i,j,2) = 1 - shares(i,j,1);
            end
        elseif pim(i,j) == 1 % white cover
            shares(i,j,1) = coverHt(i,j,1);
            shares(i,j,2) = 0;
        elseif pim(i,j) == 2
            shares(i,j,1) = 0;
            shares(i,j,2) = coverHt(i,j,2);
        end      
    end
end

% figure;imshow(shares(:,:,1),[]);
% figure;imshow(shares(:,:,2),[]);

%--------------------------------
% Stack to recover secret
%--------------------------------
%--------------------------------------------------------------
% print and stacking
% 模拟打印和叠加前，需要将黑色用1表示
sharesInv = ~shares;           % for printing&stacking, interpret '1' as black
stacked = (sharesInv(:,:,1) | sharesInv(:,:,2));  % simulate stacking ( OR operation on printed dots)
% figure; imshow(stacked,[]); title('printed');


% 在电脑屏幕上显示解码结果，需要将黑色用0表示
stacked = ~stacked;             % to display on screen, change '1' to '0' for black
% figure; imshow(stacked,[]); title('displayed on screen');

%sum(sum(stacked(secretImg==1)))/sum(sum(secretImg==0))



