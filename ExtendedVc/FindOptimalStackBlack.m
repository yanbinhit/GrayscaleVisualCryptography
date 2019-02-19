function ms = FindOptimalStackBlack(m, f)
% ����������ֵʸ����ͬ���ȣ����п��Ըı�����m��1��λ�ã�������f���ɸı䡣�ı�m��
% �Ӷ���֤���ӣ��򣩺��ɫ���أ�1����ࡣ
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

if sum(m | f)>=N % ���Ӻ�ȫ�ڣ����޸�m
    ms = m;
elseif sum( m | f) == 0 % ���Ӻ�ȫ�ף���û�п��޸ĵĺ����أ������޸�
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