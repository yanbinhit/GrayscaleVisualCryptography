function ErrorDiffusionExample()
close all; clear all; 

img = zeros(128, 256);
for i = 1:256
    img(:,i) = 256-i;
end

% Using Floyd-Steinberg kernel
imageHT = HalftoningED(img);
figure;
imshow(imageHT,[]);

% Using Jarvis kernel
imageHT = HalftoningEdJarvis(img);
figure;
imshow(imageHT, []);

