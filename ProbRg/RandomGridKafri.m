function [shares, stacked] = RandomGridKafri(img, type)
% Implmentation of the basic random grid VC proposed by Kafri, 3 types are
% implemented
%<Inputs>

img = img>0; % convert to binary {0,1}
switch type
    case 1
        [shares, stacked] = Kafri1(img);
    case 2
        [shares, stacked] = Kafri2(img);
    case 3
        [shares, stacked] = Kafri3(img);
    otherwise
        error('Wrong type! Type should be 1, 2 or 3');
end

%================================================
% Local functions
%================================================
%================================================
% Algorithm 1
%================================================
function [shares, stacked] = Kafri1(img)
[M, N] = size(img);
shares = zeros(M, N, 2);
stacked = zeros(M, N, 1);
shares(:,:, 1) = randi(2, M, N) -1;
for i = 1:M
    for j = 1:N
        if img(i,j) == 1 % white
            shares(i,j,2) = shares(i,j,1); % copy
        elseif img(i,j) == 0 % black
            shares(i,j,2) = 1 - shares(i,j,1); % complement
        else
            error('Secrect image should be binary with range {0,1}');
        end
    end
end
stacked = ~((~shares(:,:,1))|(~shares(:,:,2)));

%======================================================
% Algorithm 2 : 
%======================================================
function [shares, stacked] =  Kafri2(img)
[M, N] = size(img);
shares = zeros(M, N, 2);
stacked = zeros(M, N, 1);
shares(:,:, 1) = randi(2, M, N) -1;
for i = 1:M
    for j = 1:N
        if img(i,j) == 1 % white
            shares(i,j,2) = shares(i,j,1); % copy
        elseif img(i,j) == 0 % black
            shares(i,j,2) = randi(2)-1; % random
        else
            error('Secrect image should be binary with range {0,1}');
        end
    end
end
stacked = ~((~shares(:,:,1))|(~shares(:,:,2)));

%======================================================
% Algorithm 3 
%======================================================
function [shares, stacked] = Kafri3(img)
[M, N] = size(img);
shares = zeros(M, N, 2);
stacked = zeros(M, N, 1);
shares(:,:, 1) = randi(2, M, N) -1;
for i = 1:M
    for j = 1:N
        if img(i,j) == 1 % white
            shares(i,j,2) = randi(2)-1; % random
        elseif img(i,j) == 0 % black
            shares(i,j,2) = 1 - shares(i,j,1); % complement
        else
            error('Secrect image should be binary with range {0,1}');
        end
    end
end
stacked = ~((~shares(:,:,1))|(~shares(:,:,2)));