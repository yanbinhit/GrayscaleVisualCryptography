function ExtendedVcNaorShamir1994()
close all; clear all;
s = imread('../Images/geometry.bmp','bmp');
s = imresize(s,[512 512]);
imshow(s);
c1 = imread('../Images/cover1.bmp','bmp');
c2 = imread('../Images/cover2.bmp','bmp');

c1 = HalftoningED(imresize(c1, [512 512]));
c2 = HalftoningED(imresize(c2, [512 512]));

% % Use binary cover. Comment this our if you use halftone cover images
% c1 = imresize(imread('../Images/sdkd2.bmp','bmp'),[512 512]);
% c2 = imresize(rgb2gray(imread('../Images/binaryRice','png'))>0, [512 512]);


figure; imshow(c1,[]);
 figure; imshow(c2,[]);

s =1- s>0; % normalize and use 1 for black
c1 =1- c1>0;
c2 =1- c2>0;

[M, N] = size(s);
s1 = zeros(2*M, 2*N);
s2 = zeros(2*M, 2*N);

M0(:,:,1) = [0 0 1 1; 1 0 1 0];
M0(:,:,2) = [0 0 1 1; 1 0 1 1];
M0(:,:,3) = [1 0 1 1; 0 0 1 1];
M0(:,:,4) = [1 0 1 1; 1 0 1 1];

M1(:,:,1) = [0 0 1 1; 1 1 0 0];
M1(:,:,2) = [0 0 1 1;1 1 1 0];
M1(:,:,3) = [1 1 1 0; 0 0 1 1];
M1(:,:,4) = [0 1 1 1; 1 1 1 0];

for i = 1:M
    for j = 1:N
        if s(i,j) == 0
            switch c1(i,j)
                case 1
                    switch c2(i,j)
                        case 1
                            M = M0(:,:,4);
                        case 0
                            M = M0(:,:,3);
                    end
                case 0
                    switch c2(i,j)
                        case 1
                            M = M0(:,:,2);
                        case 0
                            M = M0(:,:,1);
                    end
            end
        elseif s(i,j) == 1
            switch c1(i,j)
                case 1
                    switch c2(i,j)
                        case 1
                            M = M1(:,:,4);
                        case 0
                            M = M1(:,:,3);
                    end
                case 0
                    switch c2(i,j)
                        case 1
                            M = M1(:,:,2);
                        case 0
                            M = M1(:,:,1);
                    end
            end
        end
        % use M to fill s1 and s2
        permVec = randperm(4);
        M = M(:,permVec);
        s1((i-1)*2+1:i*2, (j-1)*2+1:j*2) = reshape(M(1,:),2,2);
        s2((i-1)*2+1:i*2, (j-1)*2+1:j*2) = reshape(M(2,:),2,2);
    end
end

% stacking
sHat = s1 | s2;

% show shares
s1 = 1-s1;
s2 = 1-s2;
figure; 
imshow(s1,[]);
figure;
imshow(s2,[]);

% show target image
figure; 
imshow(1-sHat,[]);

