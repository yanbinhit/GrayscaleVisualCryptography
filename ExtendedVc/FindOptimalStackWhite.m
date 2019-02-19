function ms = FindOptimalStackWhite(m, f)
% 输入两个二值矢量，同长度，其中可以改变向量m中1的位置，而向量f不可改变。改变m，
% 从而保证叠加（或）后白色像素（0）最多。
%
% <inputs>
% m: length N binary vector, can be modified, 0: white, 1:black
% f: length N binary vector, fixed,
%
% <outputs>
% ms: m after swapping of 1's
%
N = length(m);

m = (m>0);
f = (f>0);


if sum(m | f)==0 % 叠加后全白，则不修改m
    ms = m;
elseif sum(m) == 4 % 没有白像素，无法修改
    ms = m;
else
    
    ww = ((m|f) == 0); % find white white positions
    wwIdx = find(~(m|f) == 0);
    bb = ((m&f) == 1); % find black black positions
    bbIdx = find(~(m&f) == 1);
    bw = 1-(ww|bb);
    bwIdx = find(bw == 1);
    ms = m;
    nBw = sum(bw); % number of black-white positions
    if nBw == 1
        ms = ms;
    else
        for i = 1:nBw-1
            for j = i+1:nBw
                if ms(bwIdx(i))*ms(bwIdx(j)) == 0
                    temp = ms(bwIdx(j));
                    ms(bwIdx(j)) = ms(bwIdx(i));
                    ms(bwIdx(i)) = temp;
                    i = j+1;
                    break;
                end
            end
        end
    end
end