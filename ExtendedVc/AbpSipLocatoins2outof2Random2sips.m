function[sips, abps] = AbpSipLocatoins2outof2Random2sips(nRows, nCols)
 %clear all; close all; 
 % 2¸ösips
Q = 4;
gamma = 2;
noAbpsPerShare = 7;
noShares = 2;

nRowsCover = nRows;
nColsCover = nCols;


sips = zeros(nRowsCover, nColsCover); % SIPs positions are marked with 1
abps = zeros(nRowsCover, nColsCover, noShares);

for i = 1:ceil(nRowsCover/Q)
    for j = 1:ceil(nColsCover/Q)
        rndVec = randperm(Q*Q);
        % sip
        sipVec = rndVec(1:gamma);
        [sipI, sipJ] = ind2sub([Q Q], sipVec);
        for k = 1:gamma
            sips((i-1)*Q+sipI(k), (j-1)*Q+sipJ(k)) = 1;
        end
        % abp1
        abp1Vec = rndVec(gamma+1:gamma+7);
        [abp1I, abp1J] = ind2sub([Q Q], abp1Vec);
        for k = 1:noAbpsPerShare
            abps((i-1)*Q+abp1I(k), (j-1)*Q+abp1J(k), 1) = 1;
        end
        % abp2
        abp2Vec = rndVec(gamma+8:gamma+14);
        [abp2I, abp2J] = ind2sub([Q Q], abp2Vec);
        for k = 1:noAbpsPerShare
            abps((i-1)*Q+abp2I(k), (j-1)*Q+abp2J(k), 2) = 1;
        end
           
        
    end
end

%  figure; imshow(sips,[]); title('sips');
%  figure; imshow(abps(:,:,1)|sips,[]); title('abps1 and sips');
%  figure;imshow((abps(:,:,1)| abps(:,:,2))  , []); title('all abps stacked');
