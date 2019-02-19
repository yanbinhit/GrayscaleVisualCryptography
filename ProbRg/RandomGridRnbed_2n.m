%function [shares, stacked] = RandomGridRnbed_2n(imageIn, n)
function RandomGridRnbed_2n()
% Implementation of the improved random grid using random noise balanced
% error diffusion(RNBED), the (2,n) scheme
%<Reference>
% Xiao-Tian Wu, et.al, Improving the visual quality of random grid-based
%   visual secret sharing via error diffusion. J. Vis. Commun. Image R.
%   2013.
% change working directory
prevDir = pwd;
[dir, dummy, dummy2] = fileparts(mfilename('fullpath'));
cd(dir);

close all; clear;
imageIn = imread('../Images/lena.tiff','tif');

%imageIn = ones(512,512)*255;
%imageIn = imread('../images/ruler512.tiff','tif');
n = 2;

% Some importan parameters used by the paper
u = 0.1;
N = [0.26 0.25 0.25 0.25 0.26];
%N = [0.27, 0.2, 0, 0.2, 0.27];


if size(imageIn,3)>1
    imageIn = rgb2gray(imageIn);
end
%imageIn = (imageIn>127).*255;
imageIn = HalftoningED(imageIn);

imageIn = padarray(imageIn, [1 1], 'replicate','both');
[nR, nC] = size(imageIn);
shares = zeros(nR, nC, n);
stacked = zeros(nR, nC);
G = zeros(nR, nC, n);

% halftoning to obtain binary secret image
sImg = HalftoningED(imageIn);
%sImg = imageIn>128;
sImg = sImg>0;

for k = 1:n
    G(:,:,k) = unifrnd(0.5-u, 0.5+u, [nR nC]);
end
gH = zeros(nR, nC, n); % halftoned G

% for i = 1:n
%     gH(:,:,i) = BlueNoiseRnbed(G(:,:,i), nR, nC, N(1));
% end

for i=2:2:nR-1
    for j=2:2:nC-1
        bt =  4 - sum(sum(sImg(i:i+1,j:j+1)));
        %=======================================================
        for jj = j:j+1;
            for ii = i:i+1;
                if sImg(i,j) == 0 % black, random
                    for kk = 1:n  % Halftone the images G with RNBED
                        sai = (2*(randi(2)-1)-1);
                        G(ii,jj,kk) = G(ii,jj,kk) + N(bt+1) * sai;
                        if G(ii,jj,kk) > 1/2
                            gH(ii,jj,kk) = 1;
                        else
                            gH(ii,jj,kk) = 0;
                        end
                        err = G(ii,jj,kk) - gH(ii,jj,kk) - N(bt+1)*sai;
                        % error diffusion
                        G(ii,jj+1,kk) = G(ii,jj+1,kk) + err * (7/16);
                        G(ii+1,jj,kk) = G(ii+1,jj,kk) + err * (5/16);
                        G(ii+1,jj-1,kk) = G(ii+1,jj-1,kk) + err * (3/16);
                        G(ii+1,jj+1,kk) = G(ii+1,jj+1,kk) + err * (1/16);
                        shares(ii,jj,kk) = gH(ii,jj,kk);
                    end
                else % sImg(i,j) = 1, white pixel, copy
                    %[d, vec] = RandomSelect(n);
                    d=1; vec=2;
                    % random generate for share d
                    sai = (2*(randi(2)-1)-1);
                    G(ii,jj,d) = G(ii,jj,d) + N(bt+1) * sai;
                    if G(ii,jj,d) > 1/2
                        gH(ii,jj,d) = 1;
                    else
                        gH(ii,jj,d) = 0;
                    end
                    err = G(ii,jj,d) - gH(ii,jj,d) - N(bt+1)*sai;
                    % error diffusion
                    G(ii,jj+1,d) = G(ii,jj+1,d) + err * (7/16);
                    G(ii+1,jj,d) = G(ii+1,jj,d) + err * (5/16);
                    G(ii+1,jj-1,d) = G(ii+1,jj-1,d) + err * (3/16);
                    G(ii+1,jj+1,d) = G(ii+1,jj+1,d) + err * (1/16);
                    shares(ii,jj,d) = gH(ii,jj,d);
                    
                    % other shares copy from d
                    for kk=vec
                        G(ii,jj,kk) = G(ii,jj,kk) + N(bt+1) * sai;
                        gH(ii,jj,kk) = gH(ii,jj,d);   % copy
                        err = G(ii,jj,kk) - gH(ii,jj,kk) - N(bt+1)*sai;
                        % error diffusion
                        G(ii,jj+1,kk) = G(ii,jj+1,kk) + err * (7/16);
                        G(ii+1,jj,kk) = G(ii+1,jj,kk) + err * (5/16);
                        G(ii+1,jj-1,kk) = G(ii+1,jj-1,kk) + err * (3/16);
                        G(ii+1,jj+1,kk) = G(ii+1,jj+1,kk) + err * (1/16);
                        shares(ii,jj,kk) = gH(ii,jj,kk);
                    end
                end
            end
        end
        
    end
end


stacked = ones(nR, nC);
for i = 1:n
    stacked = ~(~shares(:,:,i) | ~stacked);
end


figure; imshow(gH(:,:,1),[]); title('gH1');
figure; imshow(gH(:,:,2),[]); title('gH2');
figure; imshow(sImg,[]);
figure; imshow(shares(:,:,1),[]); title('share 1');
figure; imshow(shares(:,:,2),[]); title('share 2');
figure; imshow(stacked,[]);

sImgRgRnbed = sImg;
stackedRgRnbed = stacked;
save dataRgRnbed.mat sImgRgRnbed stackedRgRnbed;

% % anisotropy
% [D, a]=ddf(uint8(sImg), 1.5, 2.5, 32);
% figure;
% secplot(a, D);
% [D, a]=ddf(uint8(stacked), 1.5, 2.5, 32);
% figure;
% secplot(a, D);
% [P,fr]=rapsd(double(sImg>0), 128, [], []);
% figure;plot(fr,P,'k');
% [P,fr]=rapsd(double(stacked>0), 128, [], []);
% figure;plot(fr,P,'k');
%=======================================
%=======================================
% utility function
%=======================================
%=======================================

% randomly choose a index d from 1...n, then return d, and the remainnig
% indeces in another vector
function [d, vec] = RandomSelect(n)
d = randi(n);
vec = 1:n;
vec(d) = [];
