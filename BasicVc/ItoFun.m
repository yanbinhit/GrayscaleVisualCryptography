function [shares, stacked, imageHT] = ItoFun(imageIn)
if size(imageIn, 3)>1
    imageIn = rgb2gray(imageIn);
end

M = 1;  % the size of the block is M-by-N
N = 1;
M1 = [1 0; 1 0]; % subscript indicate birghtness
M0 = [1 0; 0 1];
[nR, nC] = size(imageIn); 
imageHT = zeros(nR+2*M, nC+2*N);
stacked = zeros(nR+2*M, nC+2*N);
imageHT = HalftoningED(imageIn);
for i = (M+1):M:(nR)
    for j = (N+1):N:(nC)     
        p = unidrnd(2);
        if imageHT(i,j) == 255
            Mp = M1(:,p);
            Ms = ~(~Mp(1,1) | ~Mp(2,1) );
        
        else
            Mp = M0(:,p);
            Ms = ~(~Mp(1,1) | ~Mp(2,1) );
        
        end  
        shares(i,j,1) = Mp(1,1);
        shares(i,j,2) = Mp(2,1);
        stacked(i,j) = Ms * 255;
            
    end
end
stacked = uint8(stacked(M+1:M+nR, N+1:N+nC));
