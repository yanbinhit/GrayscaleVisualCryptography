function DisplayDitherMatrix()

clear all; close all;
L = 4;
I = [14 10 11 15; ...
    9, 3, 0, 4;...
    8, 2, 1, 5;...
    13, 7, 6, 12];
figure;
imagesc(I);
colormap gray;
axis off;
axis square;


%% L=8

I = [62  57  48  36  37  49  58  63; ...
   56  47  35  21  22  38  50  59;  ...
   46  34  20  10  11  23  39  51;  ...
   33  19  9  3  0  4  12  24;  ...
   32  18  8  2  1  5  13  25;  ...
   45  31  17  7  6  14  26  40;  ...
   45  44  30  16  15  27  41  52;  ...
   61  54  43  29  28  42  53  60]
figure;
imagesc(I);
colormap gray;
axis off;
axis square;




