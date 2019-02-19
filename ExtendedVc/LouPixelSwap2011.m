function  [shares, stacked] = ...
    LouPixelSwap2011(cover, secretImgIn, nRowsCover, nColsCover)
% Implement the pixel swapping based size-invariant extended visual cryptography
% No authentication capability.
% 
% <Inputs>
% cover: cover images, each slice of cover is a gray-level cover image
% secretImgIn: secret image, binary
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
%
% <References>
% Section 3.1 in "2011 A novel authenticatable color visual secret sharing 
%               scheme using non-expanded meaningful shares "


% get halftone cover images
ha = HalftoningED(cover(:,:,1));
hb = HalftoningED(cover(:,:,2));

shares = zeros(nRowsCover, nColsCover, 2);
shares(:,:,1) = ha;
shares(:,:,2) = hb;

% secret image
secretImgIn = imresize(secretImgIn, [nRowsCover nColsCover], 'nearest');
secretImg = secretImgIn > 0;

% make sure that the size of ha, hb and s should be the same
% skip
nBlkRow = floor(nRowsCover/2);
nBlkCol = floor(nColsCover/2);

h = waitbar(0,'Please wait...');
for i = 1: nBlkRow
    waitbar(i/nBlkRow,h)
    for j = 1:nBlkCol
        blkIdxR =((i-1)*2+1): (i*2);
        blkIdxC =((j-1)*2+1): (j*2);
        sBlk = secretImg(blkIdxR, blkIdxC);
        haBlk = ha(blkIdxR, blkIdxC);
        hbBlk = hb(blkIdxR, blkIdxC);
        
        if sum(haBlk(:)) < sum(hbBlk(:)) % modify the block with more black pixels
            mBlk = haBlk;
            fBlk = hbBlk;
        else
            mBlk = hbBlk;
            fBlk = haBlk;
        end
        
        mVec = mBlk(:); %转成一维处理，应该更便于处理
        fVec = fBlk(:);  
        
        if (sum(sBlk(:)) == 0 ) % for all black secret block,modify M to make the
                                % stacked resut having more blacks
            mVec = FindOptimalStackBlack(~mVec, ~fVec);
            mVec = ~mVec * 255;
        elseif (sum(sBlk(:)) == 4) % for all white secret block,modify M to make 
                                    % the stacked resut having more white
            mVec = FindOptimalStackWhite(~mVec, ~fVec);
            mVec = ~mVec * 255;
        end
        
        % 填回
        if sum(haBlk(:)) < sum(hbBlk(:)) 
            shares(blkIdxR, blkIdxC, 1) = reshape(mVec, 2,2);
        else 
            shares(blkIdxR, blkIdxC, 2) = reshape(mVec, 2,2);
        end
        
    end
end
close(h);
% figure; imshow(shares(:,:,1),[]);
% figure; imshow(shares(:,:,2),[]);
%--------------------------------------------------------------
% print and stacking
% 模拟打印和叠加前，需要将黑色用1表示
sharesInv = ~shares;           % for printing&stacking, interpret '1' as black
stacked = (sharesInv(:,:,1) | sharesInv(:,:,2));  % simulate stacking ( OR operation on printed dots)
% figure; imshow(stacked,[]); title('printed');


% 在电脑屏幕上显示解码结果，需要将黑色用0表示
stacked = ~stacked;             % to display on screen, change '1' to '0' for black
% figure; imshow(stacked,[]); title('displayed on screen');


% Find the area where both the two cover images are white
whiteAreas = (ha & hb);
% figure; imshow(whiteAreas,[]);
%   

