function BlueNoiseRnbedTest()
M = 512;
N = M;
nR = 0.25;
mask = BlueNoiseRnbed(M, N, nR);
figure; imshow(mask,[]);