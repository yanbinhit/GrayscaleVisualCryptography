function ms = FindOptimalStackBlack(m, f)
% 输入两个二值矢量，同长度，其中可以改变向量m中1的位置，而向量f不可改变。改变m，
% 从而保证叠加（或）后黑色像素（1）最多。
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

if sum(m | f)>=N % 叠加后全黑，则不修改m
    ms = m;
elseif sum( m | f) == 0 % 叠加后全白，则没有可修改的黑像素，放弃修改
    ms = m;
else
    
    ww = ((m|f) == 0); % find white white positions
    wwIdx = find((m|f) == 0);
    bb = ((m&f) == 1); % find black black positions
    bbIdx = find((m&f) == 1);
    ms = m;
    if sum(ww)<=sum(bb)
        ms(wwIdx) = 1;
        ms(bbIdx(1:sum(ww))) = 0;
    else
        ms(bbIdx) = 0;
        ms(wwIdx(1:sum(bb))) = 1;
    end
end